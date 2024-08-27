# adapted from [K-SpecPart](https://github.com/TILOS-AI-Institute/HypergraphPartitioning/tree/main/K_SpecPart) under [BSD license](https://github.com/TILOS-AI-Institute/HypergraphPartitioning/blob/main/LICENSE)

include("../utils.jl")

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

# kspecpart variant
function hgr_laplacian(hypergraph::__hypergraph__, x::AbstractArray, epsilon::Int)
    eind = hypergraph.eind
    eptr = hypergraph.eptr
    n = length(x)
    m = hypergraph.num_hyperedges
    y = zeros(Float64, n)
    w = hypergraph.hwts

    for j in 1:m
        first_valid_entry = eptr[j]
        first_invalid_entry = eptr[j+1]
        k = first_invalid_entry - first_valid_entry
        scale = (floor(k/2) * ceil(k/2))/(k-1)
        sm = 0.0
        for t in first_valid_entry:first_invalid_entry-1
            sm += x[eind[t]]
        end
        sm /= k
        for t in first_valid_entry:first_invalid_entry-1
            idx = eind[t]
            #y[idx] += w[j] * (x[idx] - sm)/scale
            y[idx] += w[j] * (x[idx] - sm)/(scale*epsilon)
        end
    end
    return y
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

#kspecpart variant
function weight_balance_laplacian(x::AbstractArray, 
    vwts::Vector{Int}, 
    multiplier::AbstractArray)
twt = sum(vwts)
n = size(x, 1)
y = zeros(Float64, n)
s = multiplier[1]/twt
kvec = vwts'x
Threads.@sync Threads.@threads for j in 1:n
y[j] += twt * ((vwts[j] * x[j]) - ((kvec * vwts[j])/twt)) * s
end
return y
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

function hint_laplacian_kspecpart(hint::AbstractArray, x::AbstractArray)
    (p1, p2) = [map(((i, p),) -> i, Iterators.filter(((i, p),) -> p == pid, enumerate(hint))) for pid in 0 : 1]
    n = length(x)
    y = zeros(n)
    d1 = ones(length(p1))
    d2 = ones(length(p2))
    t1 = Threads.@spawn (sum(d2) .* d1 .* x[p1] - d1 * (d2' * x[p2]))
    t2 = Threads.@spawn (sum(d1) .* d2 .* x[p2] - d2 * (d1' * x[p1])) 
    t1 = fetch(t1)
    t2 = fetch(t2)
    y[p1] = t1
    y[p2] = t2
    return y
end

function laplacianize_adj_mat(adj::SparseMatrixCSC, graph::Union{__hypergraph__, Nothing} = nothing)
    degree(v) = sum(adj.nzval[adj.colptr[v] : adj.colptr[v + 1] - 1])
    degs = zeros(adj.n)
    #=Threads.@sync Threads.@threads=# for i in 1 : adj.n
        degs[i] = degree(i)
    end
    
    inform_dbg(adj.n, true, "degrees calculated")

    (is, js, vs) = findnz(adj)
    vs_ = -vs

    inform_dbg(adj.n, true, "matrix prepared")

    return sparse(vcat(is, 1 : adj.n), vcat(js, 1 : adj.n), vcat(vs_, degs), adj.n, adj.n)
end

function graph_adj_matrix(g::__hypergraph__)
    is = g.eind[1 : 2 : end]
    js = g.eind[2 : 2 : end]
    vs = g.hwts
    res = sparse(vcat(is, js), vcat(js, is), convert(AbstractArray{Float64, 1}, vcat(vs, vs)), g.num_vertices, g.num_vertices)

    return res
end

function adjacency_matrix(g::__hypergraph__)
    inform_dbg(g.num_vertices, true, "building adjaciency matrix...")
    if check_hypergraph_is_graph(g) && !config_approximateGraphs
        return graph_adj_matrix(g)
    else
        return hypergraph2graph(g, config_randLapCycles)
    end
end