using LinearAlgebra
using LinearMaps
using IterativeSolvers
# using GraphSignals
using SparseArrays

include("hypergraph.jl")
include("cmg/CombinatorialMultigrid.jl")
include("graphification.jl")

function LinearAlgebra.ldiv!(c::AbstractVecOrMat{T}, 
        P::CombinatorialMultigrid.lPreconditioner, 
        b::AbstractVecOrMat{T}) where {T}
    (n, m) = size(b)
    for k in 1:m
        c[:, k] = P.f(b[:, k])
    end
end

function hgr_laplacian(hg::__hypergraph__, x::AbstractArray)
    n = hg.num_vertices
    m = hg.num_hyperedges
    eind = hg.eind
    eptr = hg.eptr
    
    factor = zeros(Float64, n)
    subtrahend = zeros(Float64, n)
    
    for e in 1:m
        first_valid_entry = eptr[e]
        first_invalid_entry = eptr[e + 1]
        k = first_invalid_entry - first_valid_entry
        
        operand_dot_e = 0.0
        for pin in first_valid_entry : first_invalid_entry - 1
            operand_dot_e += x[eind[pin]]
        end
        
        clique_edge_weight = convert(Float64, hg.hwts[e]) / (convert(Float64, k) - 1.0)
        
        for pin in first_valid_entry : first_invalid_entry - 1
            factor[eind[pin]] += clique_edge_weight * convert(Float64, k)
            subtrahend[eind[pin]] += clique_edge_weight * operand_dot_e
        end
    end

    return [x[v] * factor[v] - subtrahend[v] for v in 1 : n]
end

function weight_balance_laplacian(hg::__hypergraph__, x::AbstractArray)
    n = hg.num_vertices
    m = hg.num_hyperedges
    eind = hg.eind
    eptr = hg.eptr
    
    w_dot_x = sum([hg.vwts[v] * x[v] for v in 1 : n])
    total_weight = sum(hg.vwts)

    return [hg.vwts[v] * (total_weight * x[v] - w_dot_x) for v in 1 : n]
end

function hint_laplacian(hint::AbstractArray, x::AbstractArray)
    n = length(hint)
    
    # partition counts
    n_0 = 0

    # partial sums of x
    x_0 = 0.0
    x_1 = 0.0
    
    for i in 1 : n
        if hint[i] == 0
            n_0 += 1
            x_0 += x[i]
        else
            x_1 += x[i]
        end
    end

    n_1 = n - n_0

    y = zeros(Float64, n)
    for i in 1 : n
        p0 = hint[i] == 0
        y[i] = x[i] * (p0 ? n_1 : n_0) - (p0 ? x_1 : x_0)
    end

    return y
end

function import_hypergraph(hgr_data::AbstractArray)
    data = [convert(Int64, x) for x in hgr_data]
    n = data[1]
    m = data[2]
    pin_list_indices = data[(2 + n + m + 1) : (2 + n + m + (m + 1))]
    return build_hypergraph(n,
        m,
        pin_list_indices .+ 1,
        data[(2 + n + m + (m + 1) + 1) : length(data)] .+ 1,
        -ones(Int, n),
        data[(2 + 1) : (2 + n)],
        data[(2 + n + 1) : (2 + n + m)])
end

function calc_degrees(hgr::__hypergraph__)
    degs = zeros(hgr.num_vertices)
    for v in hgr.eind
        degs[v] += 1
    end
    return degs
end

function solve_lobpcg(hgr_data::AbstractArray, hint::AbstractArray)
    hgr = import_hypergraph(hgr_data)
    lap_matrix = sparse(reduce(hcat, [hgr_laplacian(hgr, [i == j ? 1. : 0. for j in 1 : hgr.num_vertices]) for i in 1 : hgr.num_vertices]))
    
    amap = LinearMap(x -> hgr_laplacian(hgr, x), issymmetric=true, hgr.num_vertices)
    bmap = LinearMap(x -> weight_balance_laplacian(hgr, x) + hint_laplacian(hint, x), hgr.num_vertices)

    evecs = Float64[]
    (pfunc, hierarchy) = CombinatorialMultigrid.cmg_preconditioner_lap(spdiagm(ones(hgr.num_vertices) ./ 1e06) + lap_matrix)
    results = lobpcg(amap, 
            bmap, 
            false, 
            1, 
            tol=1e-40,
            # maxiter=20, 
            P = CombinatorialMultigrid.lPreconditioner(pfunc), 
            log = true)
    evecs = results.X
    
    n = length(evecs)
    @info "$results"
    @info "$n:$evecs"

    return convert(AbstractArray{Float64}, evecs)
end

function test_julia_from_c(args...)
    print("hello world!")
    return convert(AbstractArray{Float64}, [42.])
end
