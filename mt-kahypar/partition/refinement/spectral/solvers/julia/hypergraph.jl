struct __hypergraph__
    num_vertices::Int
    num_hyperedges::Int
    eptr::Vector{Int} # pin list indices
    eind::Vector{Int} # pin lists
    vptr::Vector{Int} # edge membership list indices
    vind::Vector{Int} # edge membership lists
    fixed::Vector{Int}
    vwts::Vector{Int}
    hwts::Vector{Int}
end

function build_hypergraph(num_vertices::Int,
                        num_hyperedges::Int,
                        eptr::Vector{Int},
                        eind::Vector{Int},
                        fixed_vertices::Vector{Int},
                        vwts::Vector{Int},
                        hwts::Vector{Int})
    vertices_list = [Vector{Int}() for i in 1:num_vertices]
    for i in 1:num_hyperedges
        first_valid_entry = eptr[i]
        first_invalid_entry = eptr[i+1]
        for j in first_valid_entry:first_invalid_entry-1
            v = eind[j]
            push!(vertices_list[v], i)
        end
    end
    vind = Int[]
    vptr = Int[]
    push!(vptr, 1)
    for vertices in vertices_list
        append!(vind, vertices)
        push!(vptr, length(vind)+1)
    end

    return __hypergraph__(num_vertices, 
                        num_hyperedges, 
                        eptr, 
                        eind, 
                        vptr, 
                        vind, 
                        fixed_vertices, 
                        vwts, 
                        hwts)
end