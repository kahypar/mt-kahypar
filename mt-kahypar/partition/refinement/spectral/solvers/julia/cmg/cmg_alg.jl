## define structures
struct lPreconditioner{T}
    f::T
    function lPreconditioner(pfunc)
        f = pfunc
        return new{typeof(f)}(f)
    end
end

struct CholT
    ld::SparseMatrixCSC{Float64, Int64}
    ldT::SparseMatrixCSC{Float64, Int64}
    d::Vector{Float64}
    p::Vector{Int64}
    invp::Vector{Int64}
end

struct HierarchyLevel
    sd::Bool
    islast::Bool
    iterative::Bool
    A::SparseMatrixCSC{Float64, Int64}
    invD::Vector{Float64}
    cI::Vector{Int64}
    nc::Int64
    n::Int64
    nnz::Int64
    #chol::CholT
    chol::LDLFactorizations.LDLFactorization{Float64, Int64, Int64, Int64}
end

struct Hierarchy
    A::SparseMatrixCSC{Float64, Int64}
    invD::Vector{Float64}
    cI::Vector{Int64}
    #chol::CholT
    chol::LDLFactorizations.LDLFactorization{Float64, Int64, Int64, Int64}
end

mutable struct LevelAux
    fwd::Bool
    rc::Int32
    repeat::Int32
    #
    sd::Bool
    islast::Bool
    iterative::Bool
    n::Int64
    nc::Int64
end

struct Workspace
    x::Vector{Float64}
    b::Vector{Float64}
    tmp::Vector{Float64}
end

## Hierarchy Structures Initialization

function init_LevelAux(H::Vector{HierarchyLevel})
    n_levels = length(H)
    LevelAux_ = Vector{LevelAux}(undef, n_levels)

    @inbounds for j in 1:n_levels
        if j==1
            repeat = Int32(1)
        elseif j==n_levels
            repeat = Int32(0)
        else
            repeat = Int32(max(floor(nnz(H[j-1].A)/nnz(H[j].A) - 1), 1))
        end

        LevelAux_[j] = LevelAux(true,1, repeat, H[j].sd, H[j].islast, H[j].iterative, H[j].n, H[j].nc)
    end

    return LevelAux_
end

function init_Workspace(H::Vector{HierarchyLevel})
    n_levels = length(H)
    Workspace_ = Vector{Workspace}(undef, n_levels)

    @inbounds for j in 1:n_levels
        n_ = size(H[j].A, 1)
        if j==n_levels && !H[j].iterative
            n_ = n_ + 1
        end
        Workspace_[j] = Workspace(zeros(n_), zeros(n_), zeros(n_))
    end

    return Workspace_
end


function init_Hierarchy(H::Vector{HierarchyLevel})
    n_levels = length(H)
    Hierarchy_ = Vector{Hierarchy}(undef, n_levels)

    @inbounds for j in 1:n_levels
        Hierarchy_[j] = Hierarchy(H[j].A, H[j].invD, H[j].cI, H[j].chol)
    end

    return Hierarchy_
end

##

function cmg_preconditioner_lap(A::SparseMatrixCSC)
    cmg_!(A, A)
end

function cmg_!(A::T, A_::T) where {T<:SparseMatrixCSC}
    A_o = A
    flag = Int64(0)
    sflag = Int64(0)
    loop = true
    sd = true
    iterative = true
    j = Int64(0)
    h_nnz = Int64(0)
    n = Int64(0)
    H = Vector{HierarchyLevel}()

    if size(A_, 1) > size(A, 1)
        sflag = 1
        sd = true
    else
        sflag = 1
        sd = false
    end

    # build up H
    while loop
        iterative = true
        n = size(A_, 1)
        # direct method for small size
        if (n < 500)
            iterative = false
            break
        end
        dA_ = Array(diag(A_))
        (cI, ~) = steiner_group(A_, dA_)
        nc = maximum(cI)
        islast = false
        A = A_ # !
        invD = 1 ./ (2 * dA_) # !
        R = sparse(cI, 1:n, ones(n), nc, n) # added for efficiency

        if nc == 1
            islast = true
            iterative = true
            flag = 1
        end

        # check for hierarchy stagnation for potentially bad reasons
        h_nnz = h_nnz + nnz(A_)
        if (nc >= n - 1) || (h_nnz > 5 * nnz(A_o))
            islast = true
            iterative = true
            flag = 3 # indicates stagnation
            @warn "CMG convergence may be slow due to matrix density. Future versions of CMG will eliminate this problem."
            break
        end

        Rt = sparse(cI, 1:n, 1, nc, n) # ! take out double
        A_ = (Rt * A) * Rt'
        push!(H, HierarchyLevel(sd, islast, iterative, A, invD, cI, nc, n, nnz(A), ldl([1.0 0; 0 1.0])))

        if sflag == 1
            sd = true
            sflag = 0
        end

        if nc == 1
            break
        end
    end

    # code for last hierarchy level
    if flag == 0
        j += 1
        B = A_[1:n - 1, 1:n - 1]
        ldlt = ldl(B)
        push!(H, HierarchyLevel(true, true, false,B, Float64[], Int64[], 0, n, nnz(ldlt.L), ldlt))
    end

    X = init_LevelAux(H)
    W = init_Workspace(H)
    M = init_Hierarchy(H)


    # create precondition function
    pfunc = make_preconditioner(M, W, X)

    return (pfunc, H)
end

function findRowSumAndDominance(A::SparseMatrixCSC)
    colptr::Vector{Int64} = A.colptr
    nzval::Vector = A.nzval
    rowval::Vector = A.rowval
    n = length(nzval)
    s = Int64(1)
    flag = Int64(0)
    k = Int64(0)
    il = Vector{Int64}(undef, n)
    jl = Vector{Int64}(undef, n)
    vl = Vector{Float64}(undef, n)
    sumR = zeros(Float64, length(colptr)-1)

    @inbounds for i in 1:length(colptr)-1
        l = colptr[i+1] - colptr[i]
        t = s + l - 1
        sum = 0
        for j in s:t
            val = nzval[j]
            sum += val
            row = rowval[j]
            k += 1
            il[k] = row
            jl[k] = i
            vl[k] = val
            if row != i && flag == 0
                if val > 0
                    flag = 2
                end
            end
        end
        sumR[i] = sum
        s += l
    end

    return sumR, flag, il, jl, vl
end

function validateInput!(A::SparseMatrixCSC)
    flag = Int64(0)
    # check symmetry
    if !issymmetric(A)
        return 1, A
    end
    # detect strict dominance && positive off diagonals
    n = size(A, 1)
    sAp = Vector{Float64}(undef, n)
    sd = Vector{Int64}(undef, n)
    dA = diag(A)
    (sA, flag, i, j, v) = findRowSumAndDominance(A)

    if flag == 2
        return 2, A
    end

    @inbounds @simd for i in 1:length(sA)
        sAp[i] = (sA[i] + abs(sA[i]))/2
        sd[i] = (sAp[i]/dA[i]) > 1e-13 ? 1 : 0
    end

    # augment by extra coordinate if strictly dominant
    if maximum(sd) > 0.0
        ex_v = -sAp[sd]
        ex_v_sum = -sum(ex_v)
        exd = length(ex_v)
        exd_c = findall(!iszero, sd) #nonzeros(sd) # get coordinates
        i_ = ones(Int64, exd) * (n+1)
        i = vcat(i, i_, exd_c, n+1)
        j = vcat(j, exd_c, i_, n+1)
        v = vcat(v, ex_v, ex_v, ex_v_sum)
        A = sparse(i, j, v, n + 1, n + 1)
    end
    return 3, A
end

"""
    function(cI, nc) = steiner_group(A, dA_)
    Steiner groups
    # Arguments
    -`A`: Laplacian
"""

function steiner_group(A::SparseMatrixCSC, dA_::Vector{Float64})
    (C, M) = findMinSparse(A)
    split_forest_!(C)
    efd = abs.(M ./ dA_)
    if minimum(efd) < 1 / 8 # low effective degree nodes found
        C = update_groups_(A, C, dA_)
    end
    #return C, efd
    (cI, nc, ~) = forest_components_(C)
end

"""
    function C1 = update_groups_(A, C, dA_)
    update groups based on nodes with low effective degree
    # Arguments
    -`A`: Laplacian
"""
function update_groups_(A::SparseMatrixCSC, C::Vector{Int64}, dA_::Vector{Float64})
    n = length(C)
    B = zeros(Float64, n)
    # B[j] is the total tree weight incident to node j
    @inbounds for i in 1:n
        if C[i] != i
            B[i] = A[i, C[i]] + B[i]
            B[C[i]] = A[i,C[i]] + B[C[i]]
        end
    end
    ndx = findall(x->x > -0.125, B ./ dA_)
    @inbounds @simd for i in 1:length(ndx)
        C[ndx[i]] = Int32(ndx[i])
    end
    return C
end # update_groups_

"""
    function[cI, nc, csizes] = forest_components_(C)
    forest components, connected components in unimodal forest
    # Arguments
    -`C`: unimodal tree
"""

function forest_components_(C::Vector{Int64})
    n = length(C)
    cI = zeros(Int64, n)
    cSizes = zeros(Int64, n)
    buffer = zeros(Int64, 100)
    ccI = 1

    @inbounds for j in 1:n
        bufferI = 1
        jwalk = j
        # tree walk
        while cI[jwalk] == 0
            cI[jwalk] = ccI
            buffer[bufferI] = jwalk
            bufferI = bufferI + 1
            if bufferI == size(buffer, 1) # for memory efficiency
                buffer_ = zeros(Int64, min(2 * size(buffer, 1), n))
                buffer_[1:size(buffer, 1)] = buffer
                buffer = buffer_
            end
            jwalk = C[jwalk]
        end # while

        bufferI = bufferI - 1
        en = C[jwalk] # end node
        if cI[en] != ccI
            @simd for i in 1:bufferI
                cI[buffer[i]] = cI[en]
            end
        else
            ccI = ccI + 1
        end
        cSizes[en] = cSizes[en] + bufferI
    end # for
    if cSizes[ccI] == 0
        ccI = ccI - 1
    end

    cSizes = cSizes[1:ccI]
    return cI, ccI, cSizes
end

"""
   function pfun = make_preconditioner(H)
Make preconditioner
"""
function make_preconditioner(H::Vector{Hierarchy}, W::Vector{Workspace}, X::Vector{LevelAux})
    pfun = b->preconditioner_i(H, W, X, b)
    return pfun

    #=
    if !X[1].sd
        pfun = b->preconditioner_i(H, W, X, b)
        return pfun
    end
    if X[1].sd
        pfun = b->preconditioner_sd(b, H, X,  W)
        return pfun
    end
    =#
end

"""
   function x = preconditioner_sd(b,H)
preconditioner sd
"""
function preconditioner_sd(b::Array{Float64}, H::Array{Hierarchy}, X::Vector{LevelAux}, W::Vector{Workspace})
    n = length(b)
    push!(b, -sum(b))
    x = b->preconditioner_i(H, W, X, b)
    x = x[1:n] + x[n + 1]
    return x
end

"""
   function C = split_forest_(C1::Array{Int})
decompose unimodal forest into low conductance components
"""

function split_forest_!(C1::Vector{Int64})
    n = length(C1)
    C = C1
    new_front = Int64(0)
    removed_ancestors = Int64(0)
    k = Int64(0)
    jwalk = Int64(0)
    jwalka = Int64(0)
    jwalkb = Int64(0)
    cut_mode = false
    startwalk = false
    ancestors_in_path = Int64(0)
    ancestors = zeros(Int64, n)
    indegree = zeros(Int64, n + 2)
    visited = falses(n) # logical sparse array
    walkbuffer = zeros(Int64, 20)
    newancestorbuff = zeros(Int64, 20)

   # compute indegrees
    for i in 1:n
        indegree[C[i]] = indegree[C[i]] + 1
    end

   # partition into clusters of small diameter

    for j in 1:n
        jwalk = j
        startwalk = true

        while (startwalk && (indegree[jwalk] == 0 ) && !visited[jwalk])
            startwalk = false
            ancestors_in_path = 0 # change over C-CMG
            k = 1
            walkbuffer[k] = jwalk
            newancestorbuff[k] = 0
            while (k <= 6 || visited[jwalk]) # +1 for c-indexing adjust
                jwalk = C[jwalk]
                walkterminated = (jwalk == walkbuffer[k]) || ((k > 1) && (jwalk == walkbuffer[k - 1]))
                if walkterminated
                    break # while
                end
                k += 1
                walkbuffer[k] = jwalk
                if visited[jwalk]
                    newancestorbuff[k] = ancestors_in_path
                else
                    ancestors_in_path = ancestors_in_path + 1
                    newancestorbuff[k] = ancestors_in_path
                end
            end

            if k > 6 # large diameter - cut
                middlek = Int64(ceil(k / 2))
                C[walkbuffer[middlek]] = walkbuffer[middlek] # cut middle edge
                indegree[walkbuffer[middlek + 1]] = indegree[walkbuffer[middlek + 1]] - 1 # update indegree

                for ik in (middlek + 1):k
                    ancestors[walkbuffer[ik]] = ancestors[walkbuffer[ik]] - ancestors[walkbuffer[middlek]]
                end

            # update ancestors and visited flag
                for ik in 1:middlek
                    visited[walkbuffer[ik]] = true
                    ancestors[walkbuffer[ik]] = ancestors[walkbuffer[ik]] + newancestorbuff[ik]
                end

            # set first vertex in new walk
                jwalk = walkbuffer[middlek + 1]
                startwalk = true
            end # end cut procedure

         # commit walk changes
            if !startwalk
                for ik in 1:k
                    ancestors[walkbuffer[ik]] = ancestors[walkbuffer[ik]] + newancestorbuff[ik]
                    visited[walkbuffer[ik]] = true
                end
            end
        end # outer while
    end

   # tree partition into clusters of high conductance
    for j in 1:n
        jwalk = j
        startwalk = true

        while startwalk && (indegree[jwalk] == 0)
            startwalk = false
            jwalkb = jwalk
            cut_mode = false
            # initialize new_front
            new_front = 0
            removed_ancestors = 0

            while true
                jwalka = C[jwalk]
                walkterminated = (jwalka == jwalk) || (jwalka == jwalkb)
                if walkterminated
                    break; # while
                end

                if (!cut_mode && (ancestors[jwalk] > 2) && (ancestors[jwalka] - ancestors[jwalk] > 2)) # possibly low conductance - make cut
                    C[jwalk] = jwalk # cut edge
                    indegree[jwalka] = indegree[jwalka] - 1
                    removed_ancestors = ancestors[jwalk]
                    new_front = jwalka
                    cut_mode = true
                end # end making cut

                jwalkb = jwalk
                jwalk = jwalka
                if cut_mode
                    ancestors[jwalk] = ancestors[jwalk] - removed_ancestors
                end
            end
            if cut_mode
                startwalk = true
                jwalk = new_front
            end
        end
    end
end # split_forest_

function findMinSparse(A::SparseMatrixCSC)
    minval = Inf
    row_ = Int64(0)
    colptr = A.colptr
    nzval = A.nzval
    rowval = A.rowval
    (~, m) = size(A)
    n = length(nzval)
    min_cols = zeros(Int64, m)
    min_vals = zeros(Float64, m)

    for i in 1:m
        r_loc_1 = colptr[i]
        r_loc_2 = colptr[i+1]
        size_l = r_loc_2-r_loc_1

        if size_l == 0
            continue
        end
        minval = Inf
        row_ = 0

        for j in r_loc_1:r_loc_2-1
            val = nzval[j]
            row = rowval[j]
            if val < minval
                minval = val
                row_ = row
            end
        end

        min_cols[i] = row_
        min_vals[i] = minval
    end
    return (min_cols, min_vals)
end

function findSumSparse(A::SparseMatrixCSC)
    colptr = A.colptr
    nzval = A.nzval
    s = 1
    m = length(colptr)-1
    k = 0
    sumVec = Vector{Float64}(undef, m)

    @inbounds for i in 1:m
        l = colptr[i+1] - colptr[i]
        t = s + l - 1
        if l > 0
            k += 1
        end
        sum = 0.0
        @simd for j in s:t
            val = nzval[j]
            sum += val
        end
        sumVec[k] = sum
        s += l
    end
    return sumVec
end


@inline function interpolate!(x::Vector{Float64}, cI::Vector{Int64}, z::Vector{Float64})
    @inbounds @simd for i in 1:length(x)
        x[i] = 0
    end

    @inbounds @simd for i in 1:length(z)
        idx = cI[i]
        x[idx] += z[i]
    end
end

function preconditioner_i(H::Vector{Hierarchy},  W::Vector{Workspace}, X::Vector{LevelAux}, b::Vector{Float64})
    level = Int64(1)
    #@inbounds W[1].b = b
    BLAS.blascopy!(length(b),b,1,W[1].b,1)

    while level > 0
        x = W[level].x
        b = W[level].b
        tmp = W[level].tmp

        invD = H[level].invD
        A = H[level].A
        cI = H[level].cI

        if X[level].islast && !X[level].iterative
            @inbounds ldiv!(x, H[level].chol,b)
            level -= 1
        elseif X[level].islast && X[level].iterative
            W[level].x .= b .* invD
            level -= 1
        elseif !X[level].islast && X[level].fwd
            repeat = X[level].repeat   #number of 'recursive' calls

            if X[level].rc > repeat
                X[level].rc = 1
                level -= 1
            else
                if X[level].rc == 1
                    W[level].x .= b.* invD
                else
                    mul!(tmp, A, x)
                    tmp .-= b
                    tmp .*= -invD
                    W[level].x .+= tmp
                end

                mul!(tmp, A, x)
                tmp .-= b
                X[level].fwd = false
                interpolate!(W[level+1].b, cI, -tmp)
                level += 1
            end
        elseif !X[level].islast && !X[level].fwd
            z = W[level+1].x
            W[level].x .+= z[cI]
            mul!(tmp, A, x)
            tmp .-= b
            tmp .*= -invD
            W[level].x .+= tmp
            X[level].rc += 1
            X[level].fwd = true
        end
    end
    
    return W[1].x
end