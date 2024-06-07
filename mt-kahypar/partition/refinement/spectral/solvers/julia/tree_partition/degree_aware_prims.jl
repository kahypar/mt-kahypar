using Graphs
using SimpleTraits
using DataStructures: IntDisjointSets, PriorityQueue, dequeue!, dequeue_pair!, enqueue!, heappop!, heappush!, in_same_set, peek, union!, find_root!
using LinearAlgebra: I, Symmetric, diagm, eigen, eigvals, norm, rmul!, tril, triu
import LinearAlgebra: Diagonal, issymmetric, mul!
using Random: AbstractRNG, GLOBAL_RNG, MersenneTwister, randperm, randsubseq!, seed!, shuffle, shuffle!
using SparseArrays: SparseMatrixCSC, nonzeros, nzrange, rowvals
import SparseArrays: blockdiag, sparse

import Base: adjoint, write, ==, <, *, ≈, convert, isless, issubset, union, intersect,
            reverse, reverse!, isassigned, getindex, setindex!, show,
            print, copy, in, sum, size, eltype, length, ndims, transpose,
            join, iterate, eltype, get, Pair, Tuple, zero
function degrees_aware_prim_mst end

@traitfn function degrees_aware_prim_mst(g::AG::(!Graphs.IsDirected), degree_threshold::Int,
    distmx::AbstractMatrix{T}=Graphs.weights(g)) where {T <: Real, U, AG <: Graphs.AbstractGraph{U}}
    nvg = Graphs.nv(g)
    pq = PriorityQueue{U, T}()
    finished = zeros(Bool, nvg)
    wt = fill(typemax(T), nvg) #Faster access time
    parents = zeros(U, nvg)
    degrees = zeros(U, nvg)
    pq[1] = typemin(T)
    wt[1] = typemin(T)

    while !isempty(pq)
        v = dequeue!(pq)
        finished[v] = true
        connection_flag = false
        for u in Graphs.neighbors(g, v)
            if degrees[u]+1 > degree_threshold && degrees[v]+1 > degree_threshold
                continue
            end
            finished[u] && continue
            
            if wt[u] > distmx[u, v]
                wt[u] = distmx[u, v] 
                pq[u] = wt[u]
                parents[u] = v
                connection_flag = true
            end
        end
        if connection_flag == false
            for u in Graphs.neighbors(g, v)
                finished[u] && continue
            
                if wt[u] > distmx[u, v]
                    wt[u] = distmx[u, v] 
                    pq[u] = wt[u]
                    parents[u] = v
                    connection_flag = true
                end
            end
        end 
    end

    return [Graphs.Edge{U}(parents[v], v) for v in Graphs.vertices(g) if parents[v] != 0]
end