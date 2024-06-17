include("../utils.jl")

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
    inform_dbg(num_vertices, true, "building hypergraph object...")
    vertices_list = [Vector{Int}() for _ in 1:num_vertices] # vertices_lists = [[Vector{Int}() for i in 1:num_vertices] for _ in 1:Threads.nthreads()]
    #=Threads.@sync Threads.@threads=# for i in 1:num_hyperedges
        first_valid_entry = eptr[i]
        first_invalid_entry = eptr[i+1]
        for j in first_valid_entry:first_invalid_entry-1
            v = eind[j]
            push!(vertices_list[v], i)#s[findfirst(x->x==Threads.threadid(), Threads.threadpooltids(Threads.threadpool()))][v], i)
        end
    end
    inform_dbg(num_vertices, true, "vertices lists computed")
    # vertices_list = [Vector{Int}() for _ in 1:num_vertices]
    # vertices_lists .|> enumerate .|> (vlist -> (vlist .|> (i_vvec -> append!(vertices_list[i_vvec[1]], i_vvec[2]))))
    # inform(num_vertices, true, "vertices lists assembled")
    # Threads.@sync Threads.@threads for i in 1:num_vertices
        # sort!(vertices_list[i])
    # end
    # inform(num_vertices, true, "vertices lists sorted")

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

function import_hypergraph(hgr_data::AbstractArray)
    data = convert(AbstractArray{Int64, 1}, hgr_data)
    n = data[1]
    m = data[2]
    pin_list_indices = data[(2 + n + m + 1) : (2 + n + m + (m + 1))]
    pin_lists = data[(2 + n + m + (m + 1) + 1) : length(data)]

    return build_hypergraph(n,
        m,
        pin_list_indices .+ 1,
        pin_lists .+ 1,
        -ones(Int, n),
        data[(2 + 1) : (2 + n)],
        data[(2 + n + 1) : (2 + n + m)])
end

function write_hypergraph(hgr::__hypergraph__, fname::Union{String, Nothing} = nothing)
    n = hgr.num_vertices
    e = hgr.num_hyperedges
    hedges = hgr.eind
    eptr = hgr.eptr
    hwts = hgr.hwts
    vwts = hgr.vwts
    fixed = hgr.fixed
    wt_flag = maximum(vwts) > 1 ? true : false
    hgr_file = isnothing(fname) ? "$config_tmpDir/$(hgr.num_vertices)_$(hgr.num_hyperedges)_$(time()).hgr" : fname
    f = open(hgr_file, "w")
    println(f, e, " ", n, " 11")

    for i in 1:e
        start_idx = eptr[i]
        end_idx = eptr[i+1]-1
        print(f, hwts[i])
        for j in start_idx:end_idx
            print(f, " ", hedges[j])
        end
        print(f, "\n")
    end
    for i in 1:n
        if wt_flag == 0
            println(f, vwts[i])
        else
            println(f, vwts[i]+1)
        end
    end

    close(f)

    if maximum(fixed) > -1
        f = open(hgr_file*".fixed", "w")
        for i in 1:n
            println(f, fixed_vtxs[i])
        end
        close(f)
    end
    
    return hgr_file
end

function check_hypergraph_is_graph(hgr::__hypergraph__)
    for i in 1 : hgr.num_hyperedges
        if hgr.eptr[i + 1] - hgr.eptr[i] > 2
            return false
        end
    end
    return true
end

function write_partition(partition::Vector{Int},
    partition_file_name::String)
    f = open(partition_file_name, "w")
    for i in 1:length(partition)
        println(f, partition[i])
    end
    close(f)
end

function golden_evaluator(hypergraph::__hypergraph__,
    num_parts::Int,
    partition::AbstractArray{Int})
    cutsize = 0
    for i in 1:hypergraph.num_hyperedges
        first_valid_entry = hypergraph.eptr[i]
        first_invalid_entry = hypergraph.eptr[i+1]
        for j in first_valid_entry+1 : first_invalid_entry-1
            v = hypergraph.eind[j]
            if (partition[hypergraph.eind[j]] != partition[hypergraph.eind[j-1]])
                cutsize += hypergraph.hwts[i]
                break
            end
        end
    end
    balance = zeros(Int, num_parts)
    for i in 1:hypergraph.num_vertices
        p = partition[i]
        balance[p+1] += hypergraph.vwts[i]
    end
    return (cutsize, balance)
end

function golden_evaluator_glob(g, k, p)
    result = golden_evaluator(g, k, p)
    if isnothing(global_best) || (result[1] < global_best[1] && check_balance(g, result[2], floor(Int, config_e * 50.)))
        global global_best = (result[1], p)
    end
    return result
end
