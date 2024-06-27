include("graph/hypergraph.jl")
include("graph/graphification.jl")
include("graph/matrices.jl")
include("cmg/CombinatorialMultigrid.jl")
include("utils.jl")

function LinearAlgebra.ldiv!(c::AbstractVecOrMat{T}, 
        P::CombinatorialMultigrid.lPreconditioner, 
        b::AbstractVecOrMat{T}) where {T}
    (n, m) = size(b)
    for k in 1:m
        c[:, k] = P.f(b[:, k])
    end
end

function make_a_op(hgr)
    return LinearMap((config_lapOpVariant == "paper" ?
            x -> hgr_laplacian(hgr, x) + config_lapOpShift * x
            : x -> hgr_laplacian(hgr, x, config_lapOpVariant)),
        hgr.num_vertices, issymmetric=true, isposdef=true)
end

function make_b_op(hgr, hint, acc)
    return LinearMap(x -> config_weightVsHint * 
                (config_weightOpVariant == "paper" ? weight_balance_laplacian(hgr, x)
                    : weight_balance_laplacian(x, hgr.vwts, acc))
            + (1.0 - config_weightVsHint) *
                (config_hintOpVariant == "paper" ? hint_laplacian(hint, x)
                    : hint_laplacian_kspecpart(hint, x))
            + config_balanceShift * x,
        hgr.num_vertices, issymmetric=true, isposdef=true)
end

function solve_lobpcg(hgr::__hypergraph__,
        hint_partition::AbstractArray{Int},
        deflation_space::AbstractArray{Float64, 2},
        nev::Int,
        adjacency_matrix::AbstractArray{Float64, 2})
    n = hgr.num_vertices
    m = hgr.num_hyperedges
    inform(n, false, () -> string(hgr))
    is_graph = check_hypergraph_is_graph(hgr)

    inform("received " * (is_graph ? "" : "hyper") * "graph with n=$n, m=$m")
    inform_dbg(n, false, () -> "received hint partiton: $hint_partition\nreceived deflation space: $deflation_space")
    
    amap = make_a_op(hgr)
    bacc = ones(size(hgr.vwts, 2))
    bmap = make_b_op(hgr, hint_partition, bacc)

    inform_dbg(n, true, "building laplacian...")
    lap_matrix = laplacianize_adj_mat(adjacency_matrix, is_graph ? hgr : nothing)
    inform(n, false, () -> pretty_print(lap_matrix))
    
    inform_dbg(n, true, "preconditioning...")
    preconditioner = nothing
    try
        t = @elapsed begin
            (pfunc, hierarchy) = CombinatorialMultigrid.cmg_preconditioner_lap(
                spdiagm(ones(n) .* config_preCondShift) + lap_matrix,
                timeout = config_stepTimeout)
            preconditioner = CombinatorialMultigrid.lPreconditioner(pfunc)
        end
        inform("preconditioning took $(t)s")
    catch e
        inform("didnt use preconditioner due to " * sprint(showerror, e))
        @print_backtrace
    end
    
    inform("launching LOBPCG...")
    t = @elapsed begin
        results = lobpcg(amap, 
            bmap, 
            false, 
            nev, 
            tol = 1e-40,
            maxiter = config_lobpcgMaxIters, 
            P = preconditioner,
            C = deflation_space,
            log = true)
    end
    evecs = results.X
    inform("LOBPCG successful in $(t)s")
    inform_dbg(n, false, () -> "result: " * string(map(x -> round(x, sigdigits = 2), evecs)))

    return reshape(evecs, n, config_numEvecs)
end
