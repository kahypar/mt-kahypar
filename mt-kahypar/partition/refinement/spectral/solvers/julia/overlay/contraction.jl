include("optimal_attempt_partitioner.jl")

function contract_vtx_weights(vwts::Vector{Int}, clusters::Vector{Int})
    cmax = maximum(clusters)
    vwts_contracted = zeros(Int, cmax) 
    for i in 1:length(clusters)
        vwts_contracted[clusters[i]] += vwts[i]
    end
    return vwts_contracted
end

function contract_hypergraph(hgraph::__hypergraph__, clusters::Vector{Int})
    hedges = hgraph.eind
    eptr = hgraph.eptr
    e = hgraph.num_hyperedges
    e_ = zeros(Int, length(hedges)+1)
    e_[eptr] .= 1
    e_ = cumsum(e_)
    e_ = e_[1:end-1]
    v_ = clusters[hedges]
    E = sparse(v_, e_, ones(Int, length(e_)))
    (v_, e_, ~) = findnz(E)
    E = sparse(v_, e_, ones(Int, length(e_)))
    E_sum = sum(E, dims=1)[1, :]
    ndx = findall(x-> x>1, E_sum)
    E = E[:,ndx] 
    w_ = hgraph.hwts
    w_ = w_[ndx]
    (n, ~) = size(E)
    r = rand(1, n)
    h = (r*E)[1, :]
    prm = sortperm(h)
    hs = h[prm] 
    ind = findall(hs[1:end-1] .!= hs[2:end])
    push!(ind, length(hs))
    w_ = w_[prm]
    w1 = cumsum(w_)
    w1[ind[2:end]] -= w1[ind[1:end-1]]
    ind = prm[ind]
    w1[prm] = w1
    E = E[:,ind]      
    (v_, e_, ~) = findnz(E)  
    w_ = w1[ind]
    loc_ = findall(e_[1:end-1] .!= e_[2:end])
    push!(loc_, length(e_))
    push!(loc_, 0)
    loc_[2:end] = loc_[1:end-1] .+ 1
    loc_[1] = 1
    num_vertices_contracted = maximum(v_)
    num_hyperedges_contracted = length(loc_)-1
    eptr_contracted = loc_
    eind_contracted = v_
    hwts_contracted = w_
    vwts_contracted = contract_vtx_weights(hgraph.vwts, clusters)
    fixed = hgraph.fixed
    fixed_contracted = -ones(Int, num_vertices_contracted)
    for i in 1:length(fixed)
        if fixed[i] > -1
            fixed_contracted[clusters[i]] - fixed[i]
        end
    end

    return build_hypergraph(num_vertices_contracted,
                        num_hyperedges_contracted,
                        eptr_contracted,
                        eind_contracted,
                        fixed_contracted,
                        vwts_contracted,
                        hwts_contracted)
end

function overlay(partitions::Vector, 
                hgraph::__hypergraph__)
    (clusters, clusters_sizes) = hyperedges_removal(partitions, hgraph)
    return contract_hypergraph(hgraph, clusters), clusters
end

function iterative_overlay(partitions::Vector, 
                hgraph::__hypergraph__)
    (clusters, clusters_sizes) = iterative_hyperedges_removal(partitions, hgraph)
    return contract_hypergraph(hgraph, clusters), clusters
end

function gen_partition_profile(hgraph::__hypergraph__, partition::Vector{Int})
    n = hgraph.num_vertices
    e = hgraph.num_hyperedges
    eptr = hgraph.eptr
    w_ = hgraph.hwts
    hedges = hgraph.eind
    nIncidentEdgesCut = zeros(Int, n)
    cutEdges = zeros(Int, e)
    cutStatistics = zeros(Int, e)
    edgePart = zeros(Int, e)
    fp = zeros(Int, n)
    S1 = zeros(Int, n)
    nodesInCutOnly = ones(Int, n)
    b_hedges = hgraph.vind
    b_eptr = hgraph.vptr

    for j in 1:e
        cutFlag = false
        
        for k in eptr[j]:eptr[j+1]-2
            if partition[hedges[k]] != partition[hedges[k+1]] 
                cutFlag = true
                break
            end
        end

        if cutFlag == true    
            edge = hedges[eptr[j]:(eptr[j+1]-1)]
            hye_parts = partition[edge]
            cutEdges[j] = 1
            edgePart[j] = -1
            s = sum(hye_parts)
        
            if s == length(edge)-1
                id = findall(iszero, hye_parts)
                fp[edge[id]] .+= w_[j]
            end
            
            if s == 1
                id = findall(x->x==1, hye_parts)
                fp[edge[id]] .+= w_[j]
            end
        else
            edgePart[j] = partition[hedges[eptr[j]]] 
        end
    end
    
    for j in 1:n
        for k in b_eptr[j]:b_eptr[j+1]-1
            if cutEdges[b_hedges[k]] == 0
                nodesInCutOnly[j] = 0
                break
            end
        end
    end

    ce = findall(!iszero, cutEdges)
    ww = 0

    for j in 1:length(ce)
        edge = ce[j]
        ww = w_[edge]
        
        for k in eptr[edge]:(eptr[edge+1]-1)
            nIncidentEdgesCut[hedges[k]] = nIncidentEdgesCut[hedges[k]] + ww
        end
    end

    return (cutEdges, nodesInCutOnly, nIncidentEdgesCut, fp, edgePart)
end

function hyperedges_removal(partitions::Vector, 
                            hgraph::__hypergraph__)
    m = length(partitions)
    union_cut = Int[]
    intersect_cut = Vector{Int}(1:hgraph.num_hyperedges)
    hyperedges = Vector{Int}(1:hgraph.num_hyperedges)
    for i in 1:m
        partition = partitions[i]
        (cut_edges, ~, ~, ~, ~) = gen_partition_profile(hgraph, partition)
        cut_edges_indices = findall(!iszero, cut_edges)
        union!(union_cut, cut_edges_indices)
        intersect!(intersect_cut, cut_edges_indices)
    end

    return island_removal(hgraph, union_cut)
end

function iterative_hyperedges_removal(partitions::Vector, 
                                    hgraph::__hypergraph__)
    m = length(partitions)
    union_cut = Int[]
    intersect_cut = Vector{Int}(1:hgraph.num_hyperedges)
    hyperedges = Vector{Int}(1:hgraph.num_hyperedges)
    partition = partitions[1]
    (cut_edges, ~, ~, ~, ~) = gen_partition_profile(hgraph, partition)
    cut_edges_indices = findall(!iszero, cut_edges)
    union!(union_cut, cut_edges_indices)
    (clusters, clusters_sizes)  = island_removal(hgraph, union_cut)
    for i in 2:m
        if length(clusters_sizes) > 500
            break
        end
        partition = partitions[i]
        (cut_edges, ~, ~, ~, ~) = gen_partition_profile(hgraph, partition)
        cut_edges_indices = findall(!iszero, cut_edges)
        union!(union_cut, cut_edges_indices)
        (clusters, clusters_sizes)  = island_removal(hgraph, union_cut)
    end
    return (clusters, clusters_sizes)
end