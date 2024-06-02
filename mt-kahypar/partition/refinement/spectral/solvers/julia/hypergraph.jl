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
    vertices_lists = [[Vector{Int}() for i in 1:num_vertices] for _ in 1:Threads.nthreads()]
    Threads.@sync Threads.@threads for i in 1:num_hyperedges
        first_valid_entry = eptr[i]
        first_invalid_entry = eptr[i+1]
        for j in first_valid_entry:first_invalid_entry-1
            v = eind[j]
            push!(vertices_lists[findfirst(x->x==Threads.threadid(), Threads.threadpooltids(Threads.threadpool()))][v], i)
        end
    end
    vertices_list = [Vector{Int}() for _ in 1:num_vertices]
    vertices_lists .|> enumerate .|> (vlist -> (vlist .|> (i_vvec -> append!(vertices_list[i_vvec[1]], i_vvec[2]))))
    Threads.@threads for i in 1:num_vertices
        sort!(vertices_list[i])
    end

    vind = Int[]
    vptr = Int[1]
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
