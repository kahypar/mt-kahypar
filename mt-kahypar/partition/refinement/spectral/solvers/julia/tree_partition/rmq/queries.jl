function lcaQuery(u::Int, v::Int, eulerLevel::Vector{Int}, eulerTour::Vector{Int}, eulerIndex::Vector{Int}, rmqSparseTable::Array{Int, 2})
    u_idx = eulerIndex[u]
	v_idx = eulerIndex[v]
    lca_idx = rmqQuery(u_idx, v_idx, eulerLevel, rmqSparseTable)
	lca = eulerTour[lca_idx]
    return lca
end

function rmqQuery(i::Int, j::Int, eulerLevel::Vector{Int}, rmqSparseTable::Array{Int, 2})
    tmp = 0
    idx = 0
    k = 0

    if i > j
        k = i
        i = j
        j = k
	end

	k = convert(Int, floor(log2(j-i+1)))

 	if eulerLevel[rmqSparseTable[i, k+1]] <= eulerLevel[rmqSparseTable[j-2^k+1, k+1]]
        idx = rmqSparseTable[i, k+1]
    else
        idx = rmqSparseTable[j-2^k+1, k+1]
	end
    return idx
end