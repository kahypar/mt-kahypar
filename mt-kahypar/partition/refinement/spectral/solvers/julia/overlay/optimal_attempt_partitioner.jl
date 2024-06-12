function write_hypergraph_ol(hgraph::__hypergraph__, fname::String)
    n = hgraph.num_vertices
    e = hgraph.num_hyperedges
    hedges = hgraph.eind
    eptr = hgraph.eptr
    hwts = hgraph.hwts
    vwts = hgraph.vwts
    fixed = hgraph.fixed
    wt_flag = maximum(vwts) > 1 ? true : false
    f = open(fname, "w")
    println(f, e, " ", n, " 11")

    for i in 1:e
        start_idx = eptr[i]
        end_idx = eptr[i+1]-1
        print(f, hwts[i])
        for j in start_idx:end_idx
            print(f, " ", hedges[j])
        end
        print(f, "\n")
    end
    for i in 1:n
        if wt_flag == 0
            println(f, vwts[i])
        else
            println(f, vwts[i]+1)
        end
    end

    close(f)

    if maximum(fixed) > -1
        f = open(fname*".fixed", "w")
        for i in 1:n
            println(f, fixed_vtxs[i])
        end
        close(f)
    end
end

function check_balance(hgraph::__hypergraph__, partition::Vector{Int}, num_parts::Int, ub_factor::Int)
    blocks = zeros(Int, num_parts)
    for i in 1:length(partition)
        blocks[partition[i]+1] += hgraph.vwts[i]
    end

    max_balance = Int(ceil((50 + ub_factor) * sum(hgraph.vwts)/100))
    for i in 1:num_parts
        if blocks[i] > max_balance
            return false
        end
    end
    return true
end

function optimal_partitioner(hgraph::__hypergraph__, num_parts::Int, ub_factor::Int)
    partition = zeros(Int, hgraph.num_vertices)
    hgr_file_name = config_tmpDir * "/" * "coarse.hgr"
    write_hypergraph_ol(hgraph, hgr_file_name)
    if (hgraph.num_hyperedges < 1500 && num_parts == 2)
        ilp_part(hgr_file_name = hgr_file_name, ub_factor = ub_factor)
        pfile = hgr_file_name * ".part." * string(num_parts)
        f = open(pfile, "r")
        itr = 0
        for ln in eachline(f)
            itr += 1
            p = parse(Int, ln)
            partition[itr] = p
        end
        close(f)
        (cutsize, ~) = golden_evaluator(hgraph, num_parts, partition)
        if check_balance(hgraph, partition, num_parts, ub_factor) == false || cutsize == 0
            hmetis(
                hgr_file_name = hgr_file_name,
                runs = 10,
                ctype = 1,
                rtype = 1,
                vcycle = 1,
                reconst = 0,
                dbglvl = 0
            )            
        end

        pfile = hgr_file_name * ".part." * string(num_parts)
        f = open(pfile, "r")
        itr = 0
        for ln in eachline(f)
            itr += 1
            p = parse(Int, ln)
            partition[itr] = p
        end
        close(f)
        
        rm_cmd = `rm $hgr_file_name`
        run(rm_cmd, wait=true)
        rm_cmd = `rm $pfile`
        run(rm_cmd, wait=true)
    elseif (hgraph.num_hyperedges < 300 && num_parts > 2)
        ilp_part(hgr_file_name = hgr_file_name, ub_factor = ub_factor)
        pfile = hgr_file_name * ".part." * string(num_parts)
        f = open(pfile, "r")
        itr = 0
        for ln in eachline(f)
            itr += 1
            p = parse(Int, ln)
            partition[itr] = p
        end
        close(f)
        (cutsize, ~) = golden_evaluator(hgraph, num_parts, partition)
        if check_balance(hgraph, partition, num_parts, ub_factor) == cutsize == 0
            hmetis(
                hgr_file_name = hgr_file_name,
                runs = 10,
                ctype = 1,
                rtype = 1,
                vcycle = 1,
                reconst = 0,
                dbglvl = 0
            )
        end

        pfile = hgr_file_name * ".part." * string(num_parts)
        f = open(pfile, "r")
        itr = 0
        for ln in eachline(f)
            itr += 1
            p = parse(Int, ln)
            partition[itr] = p
        end
        close(f)
        
        rm_cmd = `rm $hgr_file_name`
        run(rm_cmd, wait=true)
        rm_cmd = `rm $pfile`
        run(rm_cmd, wait=true)
    else
        # 10 parallel runs of hMETIS
        parallel_runs = 10
        runs = 10
        ctype = 1
        rtype = 1
        vcycle = 1
        reconst = 0
        dbglvl = 0
        partitions = [zeros(Int, length(partition)) for i in 1:parallel_runs]
        cutsizes = zeros(Int, parallel_runs)
        @sync Threads.@threads for i in 1:parallel_runs
            local_hgr_name = hgr_file_name * "." * string(i)
            cmd = "cp " * hgr_file_name * " " * local_hgr_name
            run(`sh -c $cmd`, wait=true)
            hmetis(
                hgr_file_name = local_hgr_name,
                ub_factor = ub_factor,
                runs = runs,
                ctype = ctype,
                rtype = rtype,
                vcycle = vcycle,
                reconst = reconst,
                dbglvl = dbglvl
            )
        end
        for i in 1:parallel_runs
            local_hgr_name = hgr_file_name * "." * string(i)
            local_pfile_name = local_hgr_name * ".part." * string(num_parts)
            f = open(local_pfile_name, "r")
            itr = 0
            for ln in eachline(f)
                itr += 1
                p = parse(Int, ln)
                partitions[i][itr] = p
            end
            close(f)
            (cutsizes[i], ~) = golden_evaluator(hgraph, num_parts, partitions[i])
            rm_cmd = `rm $local_hgr_name $local_pfile_name`
            run(rm_cmd, wait=true)
        end
        ~, best_cut_idx = findmin(cutsizes)
        rm_cmd = `rm $hgr_file_name`
        run(rm_cmd, wait=true)
        partition = partitions[best_cut_idx]
    end
    return partition
end

function hmetis(; kwargs...)
    log_file = "$config_tmpDir/hmetis_log_$(time())"
    hmetis_string = "$EXTEND_PATH_COMMAND hmetis $(kwargs[:hgr_file_name]) $config_k $(kwargs[:ub_factor]) $(kwargs[:runs]) $(kwargs[:ctype]) $(kwargs[:rtype]) $(kwargs[:vcycle]) $(kwargs[:reconst]) $(kwargs[:dbglvl]) > $log_file"
    hmetis_command = `sh -c $hmetis_string`
    run(hmetis_command, wait=true)
    run(`rm $log_file`)
end

function ilp_part(; kwargs...)
    log_file = "$config_tmpDir/ilp_part_log_$(time())"
    ilp_string = EXTEND_PATH_COMMAND * " ilp_part" * " " * kwargs[:hgr_file_name] * " $config_k $(kwargs[:ub_factor]) > $log_file"
    ilp_command = `sh -c $ilp_string`
    inform("running ilp...")
    run(ilp_command, wait = true)
    run(`rm $log_file`)
end
