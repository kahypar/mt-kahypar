function knuthshuffle!(v::AbstractVector)
    for i in length(v):-1:2
        j = rand(1:i)
        v[i], v[j] = v[j], v[i]
    end
end

function hypergraph2graph(hgraph::__hypergraph__,
                        cycles::Int)
    num_vertices = hgraph.num_vertices
    num_hyperedges = hgraph.num_hyperedges
    eptr = hgraph.eptr
    eind = hgraph.eind
    eptr_last_val = eptr[end]
    grph_i = zeros(Int, 2*cycles*eptr_last_val)
    grph_j = zeros(Int, 2*cycles*eptr_last_val)
    grph_w = zeros(Float64, 2*cycles*eptr_last_val)
    ptr = 1
    he_sizes = eptr[2:end] - eptr[1:end-1]
    max_pin_size = maximum(he_sizes)
    hwts = hgraph.hwts

    for k in 1:num_hyperedges
        he_wt = hwts[k]
        first_valid_entry = eptr[k]
        first_invalid_entry = eptr[k+1]
        he_size = first_invalid_entry - first_valid_entry
        if he_size == 2
            grph_i[ptr] = eind[first_valid_entry]
            grph_j[ptr] = eind[first_valid_entry+1]
            grph_w[ptr] = he_wt
            ptr += 1
        elseif he_size == 3
            grph_i[ptr] = eind[first_valid_entry]
            grph_j[ptr] = eind[first_valid_entry+1]
            grph_w[ptr] = he_wt/2.0
            ptr += 1
            grph_i[ptr] = eind[first_valid_entry+1]
            grph_j[ptr] = eind[first_valid_entry+2]
            grph_w[ptr] = he_wt/2.0   
            ptr += 1
            grph_i[ptr] = eind[first_valid_entry+2]
            grph_j[ptr] = eind[first_valid_entry]
            grph_w[ptr] = he_wt/2.0  
            ptr += 1 
        else 
            scale = (floor(he_size/2)*ceil(he_size/2))/(he_size-1)
            cw = 1/(cycles*2*scale)
            hedge = eind[first_valid_entry:first_invalid_entry-1]
            p = Vector{Int}(1:he_size)
            for t in 1:cycles
                p = randperm(he_size)
                grph_i[ptr:ptr+he_size-2] = hedge[p[1:end-1]]
                grph_j[ptr:ptr+he_size-2] = hedge[p[2:end]]
                grph_w[ptr:ptr+he_size-2] .= cw * he_wt
                ptr += he_size-1
                grph_i[ptr] = eind[p[end]]
                grph_j[ptr] = eind[p[1]]
                grph_w[ptr] = cw*he_wt
                ptr += 1
            end
        end
    end
    ptr -= 1
    grph_i = grph_i[1:ptr]
    grph_j = grph_j[1:ptr]
    grph_w = grph_w[1:ptr]
    n = max(maximum(grph_i), maximum(grph_j))
    adj_matrix = sparse(grph_i, 
                        grph_j, 
                        grph_w, 
                        n, 
                        n)
    return adj_matrix + adj_matrix'
end