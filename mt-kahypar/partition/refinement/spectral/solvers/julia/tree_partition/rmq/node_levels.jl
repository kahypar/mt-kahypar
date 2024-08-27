# adapted from [K-SpecPart](https://github.com/TILOS-AI-Institute/HypergraphPartitioning/tree/main/K_SpecPart) under [BSD license](https://github.com/TILOS-AI-Institute/HypergraphPartitioning/blob/main/LICENSE)

function genNodeLevel(eulerTour::Vector{Int}, n::Int)
	level = -1
	levelVec = zeros(Int, n)
	seen = zeros(Int, n)

	@inbounds for vtx in eulerTour
		@inbounds if seen[vtx] == 0
			level += 1
			@inbounds levelVec[vtx] = level
			@inbounds seen[vtx] = 1
		else
			level -= 1
		end
	end
	return levelVec
end