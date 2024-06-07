function lca2rmq(g::SimpleWeightedGraphs.SimpleGraph, cM::Int)
    eulerVecs = euler(g, cM)
    eulerLevel = zeros(Int, length(eulerVecs[1]))
    levelVec = genNodeLevel(eulerVecs[1], SimpleWeightedGraphs.nv(g))
    for idx in 1:length(eulerVecs[1])
        eulerLevel[idx] = levelVec[eulerVecs[1][idx]]
    end
	rmqSparseTable = rmq_solve(eulerLevel)
    return (rmqSparseTable, eulerLevel, eulerVecs, levelVec)
end