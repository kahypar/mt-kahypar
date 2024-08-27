include("main.jl")
include("/home/julian/Dokumente/Studium/BA/kspecpart/HypergraphPartitioning/K_SpecPart/io.jl")

# hgr_data = [7,4,1,1,1,1,1,1,1,1,1,1,1,0,2,6,9,12,0,2,0,1,3,4,3,4,6,2,5,6]
# hgr_data = [7, 4, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 2, 6, 9, 12, 0, 2, 0, 1, 3, 6, 3, 5, 6, 2, 4, 5]
# hgr_data = [8, 22, 4, 2, 5, 3, 1, 6, 2, 1, 1, 2, 1, 1, 2, 1, 2, 2, 3, 2, 1, 2, 2, 5, 1, 3, 2, 2, 2, 6, 5, 6, 0, 2, 4, 6, 8, 10, 12, 14, 16, 18, 20, 22, 24, 26, 28, 30, 32, 34, 36, 38, 40, 42, 44, 0, 4, 0, 2, 0, 1, 1, 0, 1, 2, 1, 3, 2, 0, 2, 1, 2, 4, 2, 3, 3, 1, 3, 2, 3, 5, 3, 6, 4, 0, 4, 2, 4, 5, 5, 3, 5, 4, 5, 6, 6, 3, 6, 5]
# hgr_data = [8, 11, 4, 2, 5, 3, 1, 6, 2, 1, 1, 2, 1, 2, 1, 3, 2, 2, 5, 2, 6, 0, 2, 4, 6, 8, 10, 12, 14, 16, 18, 20, 22, 0, 4, 0, 2, 0, 1, 1, 2, 1, 3, 2, 4, 2, 3, 3, 5, 3, 6, 4, 5, 5, 6]
#hint = [1, 1, 1, 0, 0, 1, 0]
# hint = [0, 0, 0, 1, 0, 1, 1, 1]
#defl = ones(hgr_data[1])
# defl = [1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 1]

hgr_data = [10, 6, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 2, 4, 6, 10, 13, 16, 0, 1, 1, 2, 3, 5, 3, 4, 6, 9, 6, 8, 9, 5, 7, 8]
# defl = ones(10)
# hint = [1, 0, 0, 1, 1, 1, 0, 0, 1, 0]

# main_lobpcg_tree_distill(hgr_data, hint, defl)
#main_lobpcg(hgr_data, hint, defl)

# hgr = import_hypergraph(hgr_data)
# n = hgr.num_vertices
# lap = sparse(reduce(hcat, [hgr_laplacian(hgr, spdiagm(ones(n))[i, 1 : n]) for i in 1 : n]))
# r = lobpcg(lap, false, 2, C = ones(Float64, n, 1))
# @info "$r"

# include("/home/julian/Dokumente/Studium/BA/test/amazon.jl")

# solve_lobpcg(amazon_hgr_data, ones(amazon_hgr_data[1]), ones(amazon_hgr_data[1])) 

#hgr_data = readlines("/home/julian/Dokumente/Studium/BA/test/sat14_atco_enc1_opt2_10_16_log.txt") .|> (x->parse(UInt64, x))

# hgr=import_hypergraph(hgr_data)
# print(hgr)

#hgr_data = readlines("/home/julian/Dokumente/Studium/BA/bt/experiments/export/hgr_data_1.718032130387584e9") .|> (x->parse(UInt64, x))
# hgr_data = readlines("/home/julian/Dokumente/Studium/BA/bt/experiments/export/hgr_data_1.718116513103965e9") .|> (x->parse(UInt64, x))
# hgr_data = readlines("/home/julian/Dokumente/Studium/BA/bt/experiments/export/1360625_data_1.718140816519633e9") .|> (x->parse(UInt64, x))
defl = ones(hgr_data[1])
hint = ones(hgr_data[1])
# hgr = import_hypergraph(hgr_data)
# queue=[1]
# visited=[]
# while length(queue)>0
#     current = pop!(queue)
#     if current in visited
#         continue
#     end
#     push!(visited, current)
#     for i_e in hgr.vptr[current]:hgr.vptr[current+1]-1
#         e = hgr.vind[i_e]
#         append!(queue,hgr.eind[hgr.eptr[e]:hgr.eptr[e+1]-1])
#     end
# end
# count = length(visited)
#hgr = read_hypergraph_file("/home/julian/Dokumente/Studium/BA/bt/experiments/export/578_4466_1.718020195936232e9.hgr")
# hint = ones(Int, hgr.num_vertices)
# defl = ones(hgr.num_vertices, 1)
config_weightVsHint = 1.0
# adj_matrix = adjacency_matrix(hgr)
# adj_matrix_deterministic = sparse([hgr_laplacian(hgr, sparsevec([i], [1], hgr.num_vertices)) for i in 1 : hgr.num_vertices])
# tree_distill(solve_lobpcg(hgr, hint, defl, config_numEvecs, adj_matrix), hgr, adj_matrix, hint)
# main_lobpcg_tree_distill(hgr_data, hint, defl)
# main_kspecpart(hgr_data, hint, defl)
main_nothing(hgr_data, hint, defl)

print("done\n")