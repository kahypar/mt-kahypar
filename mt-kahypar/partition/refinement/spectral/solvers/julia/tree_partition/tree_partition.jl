# adapted from [K-SpecPart](https://github.com/TILOS-AI-Institute/HypergraphPartitioning/tree/main/K_SpecPart) under [BSD license](https://github.com/TILOS-AI-Institute/HypergraphPartitioning/blob/main/LICENSE)

using Combinatorics
using Laplacians 
using SimpleGraphs
include("cut_distillation.jl")
include("metis.jl")
include("degree_aware_prims.jl")
include("extract_hypergraph.jl")

function find_labels(clusters::AbstractArray, n::Int)
    labels = zeros(Int, n)
    for i in 1:length(clusters)
        for j in 1:length(clusters[i])
            labels[clusters[i][j]] = i-1
        end
    end
    return labels
end

function reweigh_graph(adj::SparseMatrixCSC, 
                    X::AbstractArray, 
                    lst::Bool)
    n = size(X, 1)
    nev = size(X, 2)
    g = SimpleWeightedGraph(adj)
    ewts = g.weights
    (is, js, vs) = findnz(g.weights)
    is_r = Int[]
    js_r = Int[]
    vs_r = Int[]
    # lower left
    offset = 0
    for index in 1:length(vs)
        i = is[index]
        j = js[index]
        if j <= i
            continue
        end
        distance = 0.0
        for d in 1:nev
            span = X[i, d] - X[j, d]
            if lst == true
                if span <= 1e-09
                    distance += 1e09
                else
                    distance += 1/(span*span) 
                end
            else
                distance += abs(span)
            end
        end
        push!(is_r, i)
        push!(js_r, j)
        push!(vs_r, convert(Int, isnan(distance) ? 1e09 : round(distance)))
    end
    return SimpleWeightedGraph(sparse(vcat(is_r, js_r), vcat(js_r, is_r), vcat(vs_r, vs_r)))
end

function reweigh_graph_with_cuts(adj::SparseMatrixCSC, 
                                hgraph::__hypergraph__,
                                X::AbstractArray, 
                                lst::Bool)
    n = size(X, 1)
    nev = size(X, 2)
    offset = 0
    g = SimpleWeightedGraph(adj)
    g_copy = deepcopy(g)
    ewts = g_copy.weights
    (~, ~, he_wts) = findnz(ewts)
    (min_wt, max_wt) = extrema(he_wts)
    norm_factor = max_wt - min_wt
    he_wts ./= norm_factor
    for i in 1:length(ewts.colptr)-1
        row_len = ewts.colptr[i+1] - ewts.colptr[i]
        for row in 1:row_len
            row_x = ewts.rowval[row+offset]
            distance = 0.0
            for d in 1:nev
                span = X[row_x, d] - X[i, d]
                if lst == true
                    if span == 0.0
                        distance += 1e09
                    else
                        distance += 1/(span*span) 
                    end
                else
                    distance += abs(span)
                end
            end
            ewts.nzval[row+offset] = distance
        end
        offset += row_len
    end
    (ii, jj, algebraic_wts) = findnz(ewts)
    (min_wt, max_wt) = extrema(algebraic_wts)
    norm_factor = max_wt - min_wt
    algebraic_wts ./= norm_factor
    he_factor = 100.0 
    algebraic_factor = 1000.0
    total_wts = he_factor .* he_wts + algebraic_factor .* algebraic_wts
    println("he wts ", he_wts[1:10])
    println("algebraic_wts ", algebraic_wts[1:10])
    println("total_wts ", total_wts[1:10])
    g_matrix = sparse(ii, jj, total_wts)
    return SimpleWeightedGraph(g_matrix)
end

function construct_tree(g::SimpleWeightedGraph, 
                        X::AbstractArray, 
                        tree_type::Int)
    tree = SimpleWeightedGraphs.SimpleGraph(g)
    n = SimpleWeightedGraphs.nv(g)
    tree_matrix = spzeros(n, n)
    if tree_type == 1   #lsst construction
        vtx_ids = Vector{Int}(1:n)
        fiedler = X[:,1]
        sorted_vtx_ids = sortperm(fiedler)
        reverse_map = zeros(Int, n)
        for i in 1:n
            reverse_map[sorted_vtx_ids[i]] = i
        end
        (i, j, w) = findnz(g.weights)
        tree_reordered_mat = sparse(sorted_vtx_ids[i], sorted_vtx_ids[j], w)
        tree_reordered_graph = SimpleWeightedGraphs.SimpleWeightedGraph(tree_reordered_mat)
        lsst = akpw(tree_reordered_graph.weights)
        (i, j, w) = findnz(lsst)
        lsst = sparse(reverse_map[i], reverse_map[j], w)
        tree = SimpleWeightedGraphs.SimpleGraph(lsst)
        tree_matrix = lsst
    elseif tree_type == 2   #mst construction
        #mst = SimpleWeightedGraphs.prim_mst(g)
        mst = degrees_aware_prim_mst(g, 10)
        i = zeros(Int, length(mst))
        j = zeros(Int, length(mst))
        w = zeros(length(mst))
        for k in 1:length(mst)
            vpair = mst[k]
            vsrc = vpair.src
            vdst = vpair.dst
            i[k] = vsrc
            j[k] = vdst
            w[k] = g.weights[vsrc, vdst]
        end
        tree = SimpleWeightedGraphs.SimpleGraph(mst)
        tree_matrix = sparse(i, j, w, n, n)
        tree_matrix += tree_matrix'
    elseif tree_type == 3   #path construction
        vtxs = sortperm(X)
        i = vtxs[1:end-1]
        j = vtxs[2:end]
        w = abs.(X[j] - X[i])
        w_z = findall(w .== 0.0)
        w[w_z] .= 1e-6
        tree_matrix = sparse(i, j, w, n, n)
        tree_matrix += tree_matrix'
        tree = SimpleWeightedGraphs.SimpleGraph(tree_matrix)
    else
        @warn "Please select correct tree type!"
    end
    return tree, tree_matrix
end

function two_way_linear_tree_sweep(T::SimpleWeightedGraphs.SimpleGraph, 
                                distilled_cuts::__cut_profile__,
                                hgraph::__hypergraph__, 
                                capacities::Vector{Int}, 
                                solns::Int)
    hwts = hgraph.hwts
    num_vertices = hgraph.num_vertices
    num_hyperedges = hgraph.num_hyperedges
    vwts = hgraph.vwts
    partition_matrix = zeros(Int, solns, num_vertices)
    vtx_cuts = distilled_cuts.vtx_cuts
    edge_cuts = distilled_cuts.edge_cuts
    edge_diff = distilled_cuts.edge_diff
    pred = distilled_cuts.pred
    edge_terminators = distilled_cuts.edge_terminators
    p = distilled_cuts.p
    forced_type = distilled_cuts.forced_type
    forced_0 = distilled_cuts.forced_0
    forced_1 = distilled_cuts.forced_1
    forced_01 = distilled_cuts.forced_01
    FB0 = distilled_cuts.FB0
    FB1 = distilled_cuts.FB1
    edge_cuts_0 = distilled_cuts.edge_cuts_0
    edge_cuts_1 = distilled_cuts.edge_cuts_1
    nforced_0 = sum(hwts[forced_0])
    nforced_1 = sum(hwts[forced_1])
    nforced_01 = sum(hwts[forced_01])
    twt = sum(vwts)
    edge_cuts[1] = num_hyperedges
    vtx_cuts[1] = 0
    exc_0 = zeros(Int, num_vertices)
    exc_1 = zeros(Int, num_vertices)
    area_part = zeros(Int, 2)
    cut_cost_0 = zeros(num_vertices)
    cut_cost_1 = zeros(num_vertices)
    ratio_cost = zeros(num_vertices)
    cut_cost = zeros(num_vertices)
    area_cost = zeros(num_vertices)
    total_cost = zeros(num_vertices)
    polarity = zeros(Int, num_vertices)
    status_flag = zeros(Int, num_vertices)
    area_util_0 = zeros(Int, num_vertices)
    area_util_1 = zeros(Int, num_vertices)

    for i in 1:num_vertices
        exc_0[i] = edge_cuts[i] + nforced_0 - FB0[i] + edge_cuts_1[i] + nforced_01
        exc_1[i] = edge_cuts[i] + nforced_1 - FB1[i] + edge_cuts_0[i] + nforced_01
    end

    exc_0[1] = exc_1[1] = num_hyperedges

    for i in 1:num_vertices
        cut_cost_0[i] = exc_0[i]
        cut_cost_1[i] = exc_1[i]
        (cut_cost[i], pol) = findmin([cut_cost_0[i], cut_cost_1[i]])
        polarity[i] = pol-1

        if pol == 0
            area_util_0[i] = vtx_cuts[i]
            area_util_1[i] = twt - vtx_cuts[i]
        else
            area_util_1[i] = vtx_cuts[i]
            area_util_0[i] = twt - vtx_cuts[i]
        end

        #=SimpleWeightedGraphs.rem_edge!(T, i, pred[i])
        comps = SimpleWeightedGraphs.connected_components(T)
        cc = find_labels(comps, num_vertices)
        (cutsize, bal) = golden_evaluator(hgraph, 2, cc)
        println("[debug] bal comparison ", bal, ", ", area_util_0[i], " and ", area_util_1[i])
        SimpleWeightedGraphs.add_edge!(T, i, pred[i])=#

        if area_util_0[i] > capacities[2] || area_util_1[i] > capacities[2]
            area_cost[i] = 1e09
            #=area_cost_0 = area_util_0[i] > capacities[2] ? area_util_0[i] - capacities[2] : 0
            area_cost_1 = area_util_1[i] > capacities[2] ? area_util_1[i] - capacities[2] : 0
            area_cost[i] = (area_cost_0 + area_cost_1) * 1000
            if (i < 10) 
                println("Iteration ", i, " cutcost ", cut_cost[i], " areacost ", area_cost[i])
            end=#
        end

        ratio_cost[i] = cut_cost[i]/(area_util_0[i] * area_util_1[i])
        total_cost[i] = cut_cost[i] + area_cost[i]
    end

    cut_idxs = sortperm(total_cost)
    cut_point = cut_idxs[1]
    
    if total_cost[cut_idxs[1]] >= 1e09
        cut_idxs = sortperm(ratio_cost)
        cut_point = cut_idxs[1]
        overflow = false
        if area_util_0[cut_point] > capacities[2] || 
            area_util_1[cut_point] > capacities[2]
            overflow = true
        end
        i = 2
        while overflow == true
            cut_point = cut_idxs[i]
            overflow = false
            if area_util_0[cut_point] > capacities[2] || 
                area_util_1[cut_point] > capacities[2]
                overflow = true
            end
            i += 1
            if i > num_vertices
                cut_point= -1
                break
            end
        end
    end
    
    #println("[debug] cutpoint ", cut_point)
    partition = -ones(Int, hgraph.num_vertices)
    cutsize = 1e09
    if cut_point > -1
        SimpleWeightedGraphs.rem_edge!(T, cut_point, pred[cut_point])
        comps = SimpleWeightedGraphs.connected_components(T)
        partition = find_labels(comps, num_vertices)
        inform("Cutsize from tree sweep $(cut_cost[cut_point])")
        (cutsize, ~) = golden_evaluator(hgraph, config_k, partition)
    else
        inform("Tree sweep failed to return a valid cut")
    end
    return (partition, cutsize, cut_point)
end

function METIS_tree_partition(T::SimpleWeightedGraphs.SimpleGraph, 
                            distilled_cuts::__cut_profile__,
                            hgraph::__hypergraph__, 
                            seed::Int,
                            metis_opts::Int,
                            num_parts::Int,
                            ub_factor::Int)
    hwts = hgraph.hwts
    num_vertices = hgraph.num_vertices
    num_hyperedges = hgraph.num_hyperedges
    vwts = hgraph.vwts
    vtx_cuts = distilled_cuts.vtx_cuts
    edge_cuts = distilled_cuts.edge_cuts
    edge_diff = distilled_cuts.edge_diff
    pred = distilled_cuts.pred
    edge_terminators = distilled_cuts.edge_terminators
    p = distilled_cuts.p
    forced_type = distilled_cuts.forced_type
    forced_0 = distilled_cuts.forced_0
    forced_1 = distilled_cuts.forced_1
    forced_01 = distilled_cuts.forced_01
    FB0 = distilled_cuts.FB0
    FB1 = distilled_cuts.FB1
    edge_cuts_0 = distilled_cuts.edge_cuts_0
    edge_cuts_1 = distilled_cuts.edge_cuts_1
    nforced_0 = sum(hwts[forced_0])
    nforced_1 = sum(hwts[forced_1])
    nforced_01 = sum(hwts[forced_01])
    exc_0 = zeros(Int, num_vertices)
    exc_1 = zeros(Int, num_vertices)
    cut_cost = zeros(num_vertices)
    T_matrix = sparse(T)

    for i in 1:hgraph.num_vertices
        if pred[i] == i
            cut_cost[i] = 1e09
            continue
        end
        exc_0[i] = edge_cuts[i] + nforced_0 - FB0[i] + edge_cuts_1[i] + nforced_01
        exc_1[i] = edge_cuts[i] + nforced_1 - FB1[i] + edge_cuts_0[i] + nforced_01
        cut_cost[i] = min(exc_0[i], exc_1[i])
        parent = pred[i]
        T_matrix[i, parent] = cut_cost[i]
        T_matrix[parent, i] = cut_cost[i]
    end
    g = SimpleWeightedGraph(T_matrix)
    gname = build_metis_graph(g, metis_opts)
    metis(gname, config_k, seed, ub_factor, metis_opts)
    run(`rm -f $gname`)
    pname = gname * ".part.$config_k"
    pfile = open(pname, "r")
    partition = zeros(Int, hgraph.num_vertices)
    partition_i = 0
    for ln in eachline(pname) 
        partition_i += 1
        partition[partition_i] = parse(Int, ln)
    end
    close(pfile)
    run(`rm -f $pname`)
    (cutsize, ~) = golden_evaluator(hgraph, config_k, partition)
    inform("Cutsize from metis  $cutsize")
    return (partition, cutsize)
end

function k_way_linear_tree_sweep(T::SimpleWeightedGraphs.SimpleGraph, 
                                distilled_cuts::__cut_profile__,
                                fixed_vertices::__pindex__,
                                hgraph::__hypergraph__, 
                                capacities::Vector{Int}, 
                                num_parts::Int,
                                solns::Int)
    recursion_levels = num_parts-1
    hypergraph_recursive = deepcopy(hgraph)
    tree_recursive = deepcopy(T)
    total_wts = sum(hgraph.vwts) 
    distilled_cuts_recursive = distilled_cuts
    min_capacity = capacities[1];
    new_max_capacity = total_wts - min_capacity 
    memento = Int[]
    no_soln = false
    capacities_recursive = [min_capacity, new_max_capacity]
    #println("Capacities starting ", capacities_recursive)
    for i in 1:recursion_levels
        (partition, cutsize, cutpoint) = two_way_linear_tree_sweep(tree_recursive, 
                                distilled_cuts_recursive,
                                hypergraph_recursive, 
                                capacities_recursive, 
                                solns)
        #println("Recursive level ", i)
        #println("Capacities recursive  ", capacities_recursive)
        #println("Cutsize recorded ", cutsize)
        if cutsize >= 1e09
            no_soln = true
            break
        end
        blocks = zeros(Int, 2)
        for j in 1:length(partition) 
            blocks[partition[j]+1] += hypergraph_recursive.vwts[j]
        end
        (~, smaller_side) = findmin(blocks)
        recursive_total_vwts = 0
        for j in 1:hypergraph_recursive.num_vertices
            if partition[j] == smaller_side-1
                hypergraph_recursive.vwts[j] = 0
            else 
                recursive_total_vwts += hypergraph_recursive.vwts[j]
            end
        end
        distilled_cuts_recursive = distill_cuts_on_tree(hypergraph_recursive, 
                                                        fixed_vertices, 
                                                        T)
        capacities_recursive = [min_capacity, recursive_total_vwts-min_capacity]
        push!(memento, cutpoint)
        tree_recursive = deepcopy(T)
    end
    
    recursive_partition = -ones(Int, hgraph.num_vertices)
    cutsize = 1e09
    if no_soln == false
        for i in 1:length(memento)
            SimpleWeightedGraphs.rem_edge!(tree_recursive, memento[i], 
                                        distilled_cuts.pred[memento[i]])
        end
        comps = SimpleWeightedGraphs.connected_components(tree_recursive)
        recursive_partition = find_labels(comps, hgraph.num_vertices)
        (cutsize, balance) = golden_evaluator(hgraph, num_parts, recursive_partition)
        inform("Cutsize from tree sweep $cutsize with balance $balance")
    else 
        inform("Tree sweep failed to return a valid cut")
    end
    return (recursive_partition, cutsize)
end

function generate_next_level(partition::Vector{Int},
                            hgraph::__hypergraph__,
                            T::SimpleWeightedGraphs.SimpleGraph,
                            original_capacities::Vector{Int},
                            capacities::Vector{Int},
                            num_parts::Int)
    balance = zeros(Int, num_parts)
    for i in eachindex(partition)
        balance[partition[i]+1] += hgraph.vwts[i]
    end
    (~, max_part) = findmax(balance)
    induced_hypergraph, clusters = extract_hypergraph(hgraph, partition, max_part-1)
    cluster_labels = findall(!iszero, clusters)
    new_min_capacity = capacities[1]
    new_max_capacity = sum(induced_hypergraph.vwts) - new_min_capacity
    new_capacities = [new_min_capacity, new_max_capacity]
    T_clustered = T[cluster_labels]
    comps = SimpleWeightedGraphs.connected_components(T_clustered)
    fixed = induced_hypergraph.fixed
    p1 = findall(x-> x == 0, fixed)
    p2 = findall(x-> x == 1, fixed)
    fixed_vertices = __pindex__(p1, p2)
    distilled_cuts_clustered = distill_cuts_on_tree(induced_hypergraph, 
                                                fixed_vertices, 
                                                T_clustered)
    return __recursive_parts__(induced_hypergraph,
                            T_clustered,
                            distilled_cuts_clustered,
                            new_capacities,
                            clusters)
end


# main function ------------------------------------------------------------------------------------
function tree_partition(adj::SparseMatrixCSC,
                        X::Array{Float64},
                        hgraph::__hypergraph__,
                        fixed_vertices::__pindex__,
                        ub_factor::Int,
                        capacities::Vector{Int})
    dims = Vector{Int}(1:size(X, 2))
    types = 2
    seed = config_seed
    partitions = []
    best_partition = zeros(Int, hgraph.num_vertices)
    best_cutsize = -1

    for type in 1:types
        lst = type == 1 ? true : false
        for i in 1:length(dims)
            evecs = collect(Combinatorics.combinations(dims, i))
            for j in eachindex(evecs)
                X_thr = X[:, evecs[j]]
                if (size(X_thr, 2) > 1 && type == 3)
                    continue
                end
                if type == 3
                    X_thr = X_thr[:,1]
                end
                inform_dbg("Using eigenvectors $(evecs[j])")
                #clique_expansion = reweigh_graph_with_cuts(adj, hgraph, X_thr, lst)
                clique_expansion = reweigh_graph(adj, X_thr, lst)
                (tree, tree_matrix) = construct_tree(clique_expansion, X_thr, type)
                #=degs = [SimpleWeightedGraphs.degree(tree, i) for i in 1:SimpleWeightedGraphs.nv(tree)]
                (max_deg, node) = findmax(degs)
                println("Node ", node, " has max degree of ", max_deg)=#
                distilled_cuts = distill_cuts_on_tree(hgraph, fixed_vertices, tree)
                cutsize_tree = 1e09
                cutsize_metis = 1e09
                kway = config_k > 2
                if kway == false
                    (partition_tree, cutsize_tree, ~) = 
                                                two_way_linear_tree_sweep(
                                                    tree, 
                                                    distilled_cuts, 
                                                    hgraph, 
                                                    capacities, 
                                                    2)
                    (partition_metis, cutsize_metis) = 
                                                METIS_tree_partition(
                                                    tree, 
                                                    distilled_cuts,
                                                    hgraph, 
                                                    seed, 
                                                    type, 
                                                    config_k, 
                                                    ub_factor)
                else 
                    partition_tree, cutsize_tree = k_way_linear_tree_sweep(tree, 
                                                    distilled_cuts,
                                                    fixed_vertices,
                                                    hgraph, 
                                                    capacities, 
                                                    config_k,
                                                    2)
                    (partition_metis, cutsize_metis) = 
                                                METIS_tree_partition(
                                                    tree, 
                                                    distilled_cuts, 
                                                    hgraph, 
                                                    seed, 
                                                    type, 
                                                    config_k, 
                                                    ub_factor)
                end
                #=if cutsize_metis < cutsize_tree
                    best_partition = partition_metis
                    best_cutsize = cutsize_metis
                else
                    best_partition = partition_tree
                    best_cutsize = cutsize_tree
                end
                push!(partitions, (best_partition, best_cutsize))=#
                push!(partitions, [partition_metis, cutsize_metis])
                if (cutsize_tree < 1e09)
                    push!(partitions, [partition_tree, cutsize_tree])
                end
            end
        end
    end
    return partitions
end