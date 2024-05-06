using LinearAlgebra
using LinearMaps
using IterativeSolvers
using SparseArrays
using Random
using GraphSignals

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

function check_hypergraph_is_graph(hgr::__hypergraph__)
    for i in 1 : hgr.num_hyperedges
        if hgr.eptr[i + 1] - hgr.eptr[i] > 2
            return false
        end
    end
    return true
end

function calc_degrees(hgr::__hypergraph__)
    degs = zeros(hgr.num_vertices)
    for v in hgr.eind
        degs[v] += 1
    end
    return degs
end

function laplacianize_adj_mat(adj::SparseMatrixCSC, graph::Union{__hypergraph__, Nothing} = nothing)
    degree(v) = sum(adj.nzval[adj.colptr[v] : adj.colptr[v + 1] - 1])#sum(adj[v, 1 : adj.n])#(isnothing(graph) ? -sum(adj[v, 1 : adj.n]) : (graph.vptr[v + 1] - graph.vptr[v]) .* weights TODO)
    degs = zeros(adj.n)
    Threads.@threads for i in 1 : adj.n
        degs[i] = degree(i)
    end
    
    # degs = adj * ones(adj.n)
    
    # degs = GraphSignals.degrees(adj)

    inform(adj.n, true, "degrees calculated")

    (is, js, vs) = findnz(adj)
    vs_ = -vs

    inform(adj.n, true, "matrix prepared")

    return sparse(vcat(is, 1 : adj.n), vcat(js, 1 : adj.n), vcat(vs_, degs), adj.n, adj.n)
end

function graph_adj_matrix(g::__hypergraph__)
    is = g.eind[1 : 2 : end]
    js = g.eind[2 : 2 : end]
    vs = g.hwts
    res = sparse(vcat(is, js), vcat(js, is), convert(AbstractArray{Float64, 1}, vcat(vs, vs)), g.num_vertices, g.num_vertices)

    return res
end

function pretty_print(A)
    str = IOBuffer()
    show(IOContext(str, :compact => false), "text/plain", A)
    return String(take!(str))
end

function inform(message::String)
    if (config_verbose)
        @info message
    end
end

function inform(graph_size::Integer, big_graphs::Bool, message::String)
    if (graph_size < config_verbose_limits[1] && !big_graphs) || (graph_size > config_verbose_limits[2] && big_graphs)
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
    inform(hgr_data[1], false, "transmitted (hyper)graph data: " * string(convert(AbstractArray{Int64}, hgr_data)))

    hgr = import_hypergraph(hgr_data)
    n = hgr.num_vertices
    m = hgr.num_hyperedges
    is_graph = check_hypergraph_is_graph(hgr)
    inform("received " * (is_graph ? "" : "hyper") * "graph with n=$n, m=$m")
    
    hint_partition = convert(AbstractArray{Int64, 1}, hint)
    deflation_space = reshape(convert(AbstractArray{Float64, 1}, deflation_evecs), n, convert(Int64, length(deflation_evecs) / n))
    inform(n, false, "received hint partiton: $hint_partition\nreceived deflation space: $deflation_space")
    
    amap = make_a_op(hgr)
    bmap = make_b_op(hgr, hint_partition)
    
    try
        inform(n, true, "building adjaciency matrix...")
        lap_matrix = spdiagm([])
        if is_graph
            lap_matrix = graph_adj_matrix(hgr)
        else
            lap_matrix = hypergraph2graph(hgr, config_randLapCycles)
        end

        inform(n, true, "building laplacian...")
        lap_matrix = laplacianize_adj_mat(lap_matrix, is_graph ? hgr : nothing)
        inform(n, false, pretty_print(lap_matrix))

        inform(n, true, "preconditioning...")
        (pfunc, hierarchy) = CombinatorialMultigrid.cmg_preconditioner_lap(spdiagm(ones(n) ./ 1e06) + lap_matrix)
        preconditioner = CombinatorialMultigrid.lPreconditioner(pfunc)
       
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
        inform("LOBPCG successful")
        inform(n, false, "result: " * string(map(x -> round(x, sigdigits = 2), evecs)))
        return convert(AbstractArray{Float64}, evecs)
    catch e
        inform("failed due to " * sprint(showerror, e))
        return ones(Float64, n)
    end
end
