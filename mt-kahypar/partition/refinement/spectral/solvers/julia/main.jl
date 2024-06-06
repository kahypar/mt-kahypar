include("utils.jl")
include("graph/hypergraph.jl")
include("lobpcg.jl")
include("tree_distill.jl")

function main(input, method)
    inform("julia launched")
    n = convert(Int, input[1][1])
    inform(n, false, () -> "transmitted (hyper)graph data: " * string(convert(AbstractArray{Int64}, input[1])))

    try
        hgr = import_hypergraph(input[1])
        hint_partition = convert(AbstractArray{Int64}, input[2])
        deflation_space = reshape(convert(AbstractArray{Float64}, input[3]), n, convert(Int, length(input[3]) / n))

        method(hgr, hint_partition, deflation_space)
    catch e
        inform("failed due to " * sprint(showerror, e))
        @print_backtrace
        
        return ones(Float64, n)
    end
end

function main_lobpcg(hgr_data::AbstractArray, hint::AbstractArray, deflation_evecs::AbstractArray)
    main((hgr_data, hint, deflation_evecs), (g, h, e) -> solve_lobpcg(g, h, e, 1))
end

function main_lobpcg_tree_distill(hgr_data::AbstractArray, hint::AbstractArray, deflation_evecs::AbstractArray)
    main((hgr_data, hint, deflation_evecs), function(g::__hypergraph__, h::AbstractArray{Int}, e::AbstractArray{Float64, 2})
        return tree_distill(solve_lobpcg(g, h, e, config_numEvecs), g) 
    end)
end
