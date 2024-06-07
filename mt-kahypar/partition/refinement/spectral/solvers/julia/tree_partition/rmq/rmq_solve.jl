function rmq_solve(eulerLevel::Vector{Int}, T=Array{Int, 2})
	n_ = length(eulerLevel)
    l_ = convert(Int, ceil(log2(n_)))+1	
    st_ = zeros(Int, n_, l_)
    k = 0
    @simd for idx in 1:n_
        st_[idx, 1] = idx
    end
    for j in 2:l_
        k = 1 << (j-2)
        for i in 1:n_
            if eulerLevel[st_[i, j-1]] <= eulerLevel[st_[min(i+k, n_), j-1]]
            	st_[i, j] = st_[i, j-1]
            else
                st_[i, j] = st_[min(i+k, n_), j-1]
            end
        end
	end
    return st_
end