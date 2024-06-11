function isolate_islands(hgraph::__hypergraph__)
    cluster_id = 0
    num_vertices = hgraph.num_vertices
    num_hyperedges = hgraph.num_hyperedges
    eptr = hgraph.eptr
    eind = hgraph.eind
    vptr = hgraph.vptr
    vind = hgraph.vind
    fixed = hgraph.fixed
    vwts = hgraph.vwts
    hwts = hgraph.hwts
    vwts_processed = similar(vwts)
    hwts_processed = similar(hwts)
    original_indices = zeros(Int, num_vertices)
    unused_indices = zeros(Int, num_vertices)
    new_indices = zeros(Int, num_vertices)
    eptr_processed = zeros(Int, length(eptr))
    eind_processed = zeros(Int, length(eind))
    (clusters, clusters_sizes) = island_removal(hgraph, Vector{Int}())
    max_cluster = findmax(clusters_sizes)[2]
    new_index_ptr = 0
    unused_index_ptr = 0

    for i in 1:length(clusters)
        if clusters[i] == max_cluster
            new_index_ptr += 1
            original_indices[new_index_ptr] = i
            new_indices[i] = new_index_ptr
            vwts_processed[new_index_ptr] = vwts[i]
        else 
            unused_index_ptr += 1
            unused_indices[unused_index_ptr] = i
        end
    end

    original_indices = original_indices[1:new_index_ptr]
    vwts_processed = vwts_processed[1:new_index_ptr]
    unused_indices = unused_indices[1:unused_index_ptr]
    cluster_ptr = 0

    for i in 1:length(eind)
        vtx = eind[i]
        vtx_cluster = clusters[vtx]
        if vtx_cluster == max_cluster
            cluster_ptr += 1
            eind_processed[cluster_ptr] = new_indices[vtx]
        end
    end

    eind_processed = eind_processed[1:cluster_ptr]
    cluster_ptr = 0
    buffer = 1
    num_vertices_processed = maximum(new_indices)
    num_hyperedges_processed = 0
    size = 0

    for k in 1:num_hyperedges
        loc = eptr[k]
        if clusters[eind[loc]] == max_cluster
            num_hyperedges_processed += 1
            cluster_ptr += 1
            size = eptr[k+1] - eptr[k]
            eptr_processed[cluster_ptr] = buffer
            hwts_processed[cluster_ptr] = hwts[k]
            buffer += size
        end
    end
    
    cluster_ptr += 1
    eptr_processed[cluster_ptr] = buffer
    eptr_processed = eptr_processed[1:cluster_ptr]
    hwts_processed = hwts_processed[1:cluster_ptr-1]
    vind_processed = zeros(length(eind_processed))
    vind_processed[eptr_processed[1:end-1]] .= 1
    vind_processed = cumsum(vind_processed)
    he_prm = sortperm(eind_processed)
    v = eind_processed[he_prm]
    vind_processed = vind_processed[he_prm] 
    vptr_processed = zeros(Int, num_vertices_processed+1)
    vptr_processed[2:end-1] = findall(v[1:end-1] .!= v[2:end]) .+ 1
    vptr_processed[1] = 1
    vptr_processed[end] = eptr_processed[end]
    fixed_processed = -ones(Int, num_vertices_processed)

    for k in 1:length(fixed)
        if fixed[clusters[k]] > -1
            fixed_processed[clusters[k]] = fixed[k]
        end
    end

    return __hypergraph__(num_vertices_processed, 
                        num_hyperedges_processed, 
                        eptr_processed, 
                        eind_processed, 
                        vptr_processed, 
                        vind_processed, 
                        fixed_processed, 
                        vwts_processed, 
                        hwts_processed), 
            original_indices, 
            new_indices,
            unused_indices
end

function island_removal(hgraph::__hypergraph__,
                        excluded_hyperedges::Vector{Int})
    eptr = hgraph.eptr
    eind = hgraph.eind
    num_hyperedges = hgraph.num_hyperedges
    num_vertices = hgraph.num_vertices
    cids = Vector{Int}(1:num_vertices)
    sz = ones(Int, num_vertices)
    included_hyperedges = setdiff(Vector{Int}(1:num_hyperedges), 
                        excluded_hyperedges)
    for k in 1:length(included_hyperedges)
        he = included_hyperedges[k]
        first_valid_entry = eptr[he]
        first_invalid_entry = eptr[he+1]
        len_span = first_invalid_entry - first_valid_entry
        for j in 0:len_span-2
            u = eind[eptr[he] + j]
            v = eind[eptr[he] + j + 1]
            while u != cids[u]
                cids[u] = cids[cids[u]]
                u = cids[u]
            end
            while v != cids[v]
                cids[v] = cids[cids[v]]
                v = cids[v]
            end
            if u != v
                if sz[u] < sz[v]
                    cids[u] = v
                    sz[v] += sz[u]
                else
                    cids[v] = u
                    sz[u] += sz[v]
                end
            end
        end
    end
    clusters = zeros(Int, num_vertices)
    clusters_sizes = zeros(Int, num_vertices)
    for k in 1:num_vertices
        u = k
        while u != cids[u]
            u = cids[u]
        end
        clusters[k] = u
        clusters_sizes[u] += 1
    end

    indices = findall(x-> x > 0, clusters_sizes)
    ordered_indices = zeros(Int, maximum(indices))
    ordered_indices[indices] = Vector{Int}(1:length(indices))

    for k in 1:num_vertices
        clusters[k] = ordered_indices[clusters[k]]
    end

    clusters_sizes = clusters_sizes[indices]
    return (clusters, clusters_sizes) 
end