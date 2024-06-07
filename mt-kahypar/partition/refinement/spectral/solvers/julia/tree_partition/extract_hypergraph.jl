function extract_hypergraph(hgraph::__hypergraph__, 
                            partition::Vector{Int}, 
                            side::Int)
    num_hyperedges = hgraph.num_hyperedges
    num_vertices = hgraph.num_vertices
    eptr = hgraph.eptr
    eind = hgraph.eind
    fixed = hgraph.fixed
    vwts = hgraph.vwts
    hwts = hgraph.hwts
    cluster_labels = zeros(Int, num_vertices)
    vwts_clustered = Int[]
    cluster_id = 0
    for i in eachindex(partition)
        if partition[i] == side
            cluster_id += 1
            cluster_labels[i] = cluster_id
            push!(vwts_clustered, vwts[i])
        end
    end

    hwts_clustered = Int[]
    hash_map = Dict{Int, Int}()
    he_cluster = 0
    hyperedges_vec = [Vector{Int}() for i in 1:num_hyperedges]
    for i in 1:hgraph.num_hyperedges
        first_valid_entry = eptr[i]
        first_invalid_entry = eptr[i+1]
        hash_val = 0
        for j in first_valid_entry:first_invalid_entry-1
            vtx = eind[j]
            clabel = cluster_labels[vtx]
            if clabel > 0
                push!(hyperedges_vec[i], clabel)
                hash_val += clabel*clabel
            end
        end
        if isempty(hyperedges_vec[i]) == true
            continue
        end
        if haskey(hash_map, hash_val)
            hashed_he = hash_map[hash_val]
            hwts_clustered[hashed_he] += hwts[i]
        else
            he_cluster += 1
            push!(hash_map, hash_val => he_cluster)
            if length(hyperedges_vec[i]) > 1
                push!(hwts_clustered, hwts[i])
            else
                push!(hwts_clustered, 0)
            end
        end
    end

    eind_clustered = Int[]
    eptr_clustered = [1]
    num_vertices_clustered = cluster_id
    num_hyperedges_clustered = he_cluster
    for i in eachindex(hyperedges_vec)
        if isempty(hyperedges_vec[i]) == true
            continue
        end
        append!(eind_clustered, hyperedges_vec[i])
        push!(eptr_clustered, length(eind_clustered)+1)
    end

    vind_clustered = zeros(length(eind_clustered))
    vind_clustered[eptr_clustered[1:end-1]] .= 1
    vind_clustered = cumsum(vind_clustered)
    he_prm = sortperm(eind_clustered)
    v = eind_clustered[he_prm]
    vind_clustered = vind_clustered[he_prm] 
    vptr_clustered = zeros(Int, num_vertices_clustered+1)
    vptr_clustered[2:end-1] = findall(v[1:end-1] .!= v[2:end]) .+ 1
    vptr_clustered[1] = 1
    vptr_clustered[end] = eptr_clustered[end]
    fixed = hgraph.fixed
    fixed_clustered = -ones(Int, num_vertices_clustered)
    for i in eachindex(fixed)
        if fixed[i] > 0
            fixed_clustered[cluster_labels[i]] = fixed[i]
        end
    end

    return __hypergraph__(num_vertices_clustered, 
                        num_hyperedges_clustered, 
                        eptr_clustered, 
                        eind_clustered, 
                        vptr_clustered, 
                        vind_clustered, 
                        fixed_clustered, 
                        vwts_clustered, 
                        hwts_clustered), cluster_labels
end 