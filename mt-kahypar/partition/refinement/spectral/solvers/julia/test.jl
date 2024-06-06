include("main.jl")

# hgr_data = [7,4,1,1,1,1,1,1,1,1,1,1,1,0,2,6,9,12,0,2,0,1,3,4,3,4,6,2,5,6]
hgr_data = [7, 4, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 2, 6, 9, 12, 0, 2, 0, 1, 3, 6, 3, 5, 6, 2, 4, 5]
# hgr_data = [8, 22, 4, 2, 5, 3, 1, 6, 2, 1, 1, 2, 1, 1, 2, 1, 2, 2, 3, 2, 1, 2, 2, 5, 1, 3, 2, 2, 2, 6, 5, 6, 0, 2, 4, 6, 8, 10, 12, 14, 16, 18, 20, 22, 24, 26, 28, 30, 32, 34, 36, 38, 40, 42, 44, 0, 4, 0, 2, 0, 1, 1, 0, 1, 2, 1, 3, 2, 0, 2, 1, 2, 4, 2, 3, 3, 1, 3, 2, 3, 5, 3, 6, 4, 0, 4, 2, 4, 5, 5, 3, 5, 4, 5, 6, 6, 3, 6, 5]
# hgr_data = [8, 11, 4, 2, 5, 3, 1, 6, 2, 1, 1, 2, 1, 2, 1, 3, 2, 2, 5, 2, 6, 0, 2, 4, 6, 8, 10, 12, 14, 16, 18, 20, 22, 0, 4, 0, 2, 0, 1, 1, 2, 1, 3, 2, 4, 2, 3, 3, 5, 3, 6, 4, 5, 5, 6]
hint = [1, 1, 1, 0, 0, 1, 0]
# hint = [0, 0, 0, 1, 0, 1, 1, 1]
defl = ones(hgr_data[1])
# defl = [1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 1]

main_lobpcg_tree_distill(hgr_data, hint, defl)

# hgr = import_hypergraph(hgr_data)
# n = hgr.num_vertices
# lap = sparse(reduce(hcat, [hgr_laplacian(hgr, spdiagm(ones(n))[i, 1 : n]) for i in 1 : n]))
# r = lobpcg(lap, false, 2, C = ones(Float64, n, 1))
# @info "$r"

# include("/home/julian/Dokumente/Studium/BA/test/amazon.jl")

# solve_lobpcg(amazon_hgr_data, ones(amazon_hgr_data[1]), ones(amazon_hgr_data[1])) 

#hgr_data =  readlines("/home/julian/Dokumente/Studium/BA/test/sat14_atco_enc1_opt2_10_16_log.txt") .|> (x->parse(UInt64, x))

# hgr=import_hypergraph(hgr_data)
# print(hgr)
