@inline euler(g::SimpleWeightedGraphs.SimpleGraph, s::Int; dir=:out) =
	(dir == :out) ? dfsEuler(g, s, SimpleWeightedGraphs.outneighbors) : dfsEuler(g, s, inneighbors)

function dfsEuler(g::SimpleWeightedGraphs.SimpleGraph{T}, s::Integer, neighborfn::Function) where T<:Int
	seen = Vector{Bool}()
	resize!(seen, SimpleWeightedGraphs.nv(g))
	fill!(seen, false)
	finishOrder = zeros(T, SimpleWeightedGraphs.nv(g))
	parentOf = zeros(T, SimpleWeightedGraphs.nv(g))
	eulerIndex = zeros(T, SimpleWeightedGraphs.nv(g))
	S = Vector{T}([s])
	eulerTour = Vector{T}()
	push!(eulerTour, s)
	eulerIndex[s] = 1
	seen[s] = true
	parentOf[s] = s
	idx = 0
	
	while !isempty(S)
		v = S[end]
		u = 0
		for n in neighborfn(g, v)
			if !seen[n]
				u = n
				break
			end
		end
		if u == 0
			idx += 1
			finishOrder[idx] = S[end]
			push!(eulerTour, parentOf[v])
			pop!(S)
		else
			seen[u] = true
			push!(eulerTour, u)
			eulerIndex[u] = length(eulerTour)
			append!(S, u)
			parentOf[u] = v
		end
	end

	eulerTour = eulerTour[1:length(eulerTour)-1]

	return (eulerTour, eulerIndex, finishOrder, parentOf)
end