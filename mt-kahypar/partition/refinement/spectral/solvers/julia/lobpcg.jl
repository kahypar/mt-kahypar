using LinearAlgebra
using LinearMaps
using IterativeSolvers
using SparseArrays
using Random

include("hypergraph.jl")
include("cmg/CombinatorialMultigrid.jl")
include("graphification.jl")

include("config.jl")

function LinearAlgebra.ldiv!(c::AbstractVecOrMat{T}, 
        P::CombinatorialMultigrid.lPreconditioner, 
        b::AbstractVecOrMat{T}) where {T}
    (n, m) = size(b)
    for k in 1:m
        c[:, k] = P.f(b[:, k])
    end
end

function hgr_laplacian(hg::__hypergraph__, x::AbstractArray)
    n = hg.num_vertices
    m = hg.num_hyperedges
    eind = hg.eind
    eptr = hg.eptr
    
    factor = zeros(Float64, n)
    subtrahend = zeros(Float64, n)
    
    for e in 1:m
        first_valid_entry = eptr[e]
        first_invalid_entry = eptr[e + 1]
        k = first_invalid_entry - first_valid_entry
        
        operand_dot_e = 0.0
        for pin in first_valid_entry : first_invalid_entry - 1
            operand_dot_e += x[eind[pin]]
        end
        
        clique_edge_weight = convert(Float64, hg.hwts[e]) / (convert(Float64, k) - 1.0)
        
        for pin in first_valid_entry : first_invalid_entry - 1
            factor[eind[pin]] += clique_edge_weight * convert(Float64, k)
            subtrahend[eind[pin]] += clique_edge_weight * operand_dot_e
        end
    end

    return [x[v] * factor[v] - subtrahend[v] for v in 1 : n]
end

function weight_balance_laplacian(hg::__hypergraph__, x::AbstractArray)
    n = hg.num_vertices
    m = hg.num_hyperedges
    eind = hg.eind
    eptr = hg.eptr
    
    w_dot_x = sum([hg.vwts[v] * x[v] for v in 1 : n])
    total_weight = sum(hg.vwts)

    return [hg.vwts[v] * (total_weight * x[v] - w_dot_x) for v in 1 : n]
end

function hint_laplacian(hint::AbstractArray, x::AbstractArray)
    n = length(hint)
    
    # partition counts
    n_0 = 0

    # partial sums of x
    x_0 = 0.0
    x_1 = 0.0
    
    for i in 1 : n
        if hint[i] == 0
            n_0 += 1
            x_0 += x[i]
        else
            x_1 += x[i]
        end
    end

    n_1 = n - n_0

    y = zeros(Float64, n)
    for i in 1 : n
        p0 = hint[i] == 0
        y[i] = x[i] * (p0 ? n_1 : n_0) - (p0 ? x_1 : x_0)
    end

    return y
end

function import_hypergraph(hgr_data::AbstractArray)
    data = convert(AbstractArray{Int64, 1}, hgr_data)
    n = data[1]
    m = data[2]
    pin_list_indices = data[(2 + n + m + 1) : (2 + n + m + (m + 1))]
    return build_hypergraph(n,
        m,
        pin_list_indices .+ 1,
        data[(2 + n + m + (m + 1) + 1) : length(data)] .+ 1,
        -ones(Int, n),
        data[(2 + 1) : (2 + n)],
        data[(2 + n + 1) : (2 + n + m)])
end

function calc_degrees(hgr::__hypergraph__)
    degs = zeros(hgr.num_vertices)
    for v in hgr.eind
        degs[v] += 1
    end
    return degs
end

function laplacianize_adj_mat(adj::SparseMatrixCSC)
    res = deepcopy(adj)
    inform("deepcopy created")
    @sync Threads.@threads for i in 1 : adj.n
        res[i, i] = -sum(adj[i, 1 : adj.n])
    end
    return -res
end

function inform(message::String)
    if (config_verbose)
        @info message
    end
end

function inform(graph_size::Integer, big_graphs::Bool, message::String)
    if (graph_size < 25 && !big_graphs) || (graph_size > 25000 && big_graphs)
        inform(message)
    end
end

function make_a_op(hgr)
    return LinearMap(x -> hgr_laplacian(hgr, x) + config_lapOpShift * x,
        hgr.num_vertices, issymmetric=true, isposdef=true)
end

function make_b_op(hgr, hint)
    return LinearMap(x -> config_weightVsHint * weight_balance_laplacian(hgr, x)
            + (1.0 - config_weightVsHint) *hint_laplacian(hint, x) + config_lapOpShift * x,
        hgr.num_vertices, issymmetric=true, isposdef=true)
end

# TODO: set number of evecs
function solve_lobpcg(hgr_data::AbstractArray, hint::AbstractArray, deflation_evecs::AbstractArray)
    hgr = import_hypergraph(hgr_data)
    n = hgr.num_vertices
    m = hgr.num_hyperedges
    inform("received hypergraph with n=$n, m=$m")
    inform(n, false, string(convert(AbstractArray{Int64}, hgr_data)))

    hint_partition = convert(AbstractArray{Int64, 1}, hint)
    deflation_space = reshape(convert(AbstractArray{Float64, 1}, deflation_evecs), n, convert(Int64, length(deflation_evecs) / n))
    inform(n, false, "received hint partiton: $hint_partition\nreceived deflation space: $deflation_space")
    
    amap = make_a_op(hgr)
    bmap = make_b_op(hgr, hint_partition)

    inform(n, true, "dehyperize to adjaciency matrix...")
    rand_adj_matrix = hypergraph2graph(hgr, config_randLapCycles)
    inform(n, true, "building laplacian...")
    rand_lap_matrix = laplacianize_adj_mat(rand_adj_matrix)
    inform(n, true, "preparing preconditioning...")
    (pfunc, hierarchy) = CombinatorialMultigrid.cmg_preconditioner_lap(spdiagm(ones(n) ./ 1e06) + rand_lap_matrix)
    inform(n, true, "preconditioning...")
    preconditioner = CombinatorialMultigrid.lPreconditioner(pfunc)
    
    evecs = Float64[]
    try
        inform("launching LOBPCG...")
        results = lobpcg(amap, 
            bmap, 
            false, 
            config_numEvecs, 
            tol = 1e-40,
            maxiter = config_lobpcgMaxIters, 
            P = preconditioner,
            C = deflation_space,
            log = true)
        evecs = results.X
        inform("successful")
    catch e
        inform("failed due to $e")
        evecs = ones(Float64, n)
    end

    rounded_evecs = map(x -> round(x, sigdigits = 2), evecs)
    inform(n, false, "result: $rounded_evecs")

    return convert(AbstractArray{Float64}, evecs)
end
