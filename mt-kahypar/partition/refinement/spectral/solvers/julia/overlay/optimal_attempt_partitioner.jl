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

    return check_balance(hgraph, blocks, ub_factor)
end

function check_balance(hgraph::__hypergraph__ ,part_weights::Vector{Int}, ub_factor::Int)
    total_weight = sum(hgraph.hwts)
    max_balance = max(floor(Int, total_weight * 0.5 * (1.0 + ub_factor / 100.)), ceil(total_weight / 2.))
    for w in part_weights
        if w > max_balance
            return false
        end
    end
    return true
end

function optimal_partitioner(hgraph::__hypergraph__, num_parts::Int, ub_factor::Int)
    partition = zeros(Int, hgraph.num_vertices)
    hgr_file_name = config_tmpDir * "/" * "coarse_$(time()).hgr"
    write_hypergraph_ol(hgraph, hgr_file_name)
    pfile = hgr_file_name * ".part.$num_parts"
    if (hgraph.num_hyperedges < 1500 && num_parts == 2 && ilp_part(hgr_file_name = hgr_file_name, ub_factor = ub_factor))
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
                dbglvl = 0,
                ub_factor = ub_factor
            )            
        end

        f = open(pfile, "r")
        itr = 0
        for ln in eachline(f)
            itr += 1
            p = parse(Int, ln)
            partition[itr] = p
        end
        close(f)
        
        rm_cmd = `rm -f $hgr_file_name`
        run(rm_cmd, wait=true)
        rm_cmd = `rm -f $pfile`
        run(rm_cmd, wait=true)
    elseif (hgraph.num_hyperedges < 300 && num_parts > 2 && ilp_part(hgr_file_name = hgr_file_name, ub_factor = ub_factor))
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
                dbglvl = 0,
                ub_factor = ub_factor
            )
        end

        f = open(pfile, "r")
        itr = 0
        for ln in eachline(f)
            itr += 1
            p = parse(Int, ln)
            partition[itr] = p
        end
        close(f)
        
        rm_cmd = `rm -f $hgr_file_name`
        run(rm_cmd, wait=true)
        rm_cmd = `rm -f $pfile`
        run(rm_cmd, wait=true)
    else
        inform("10 parallel runs of hMETIS")
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
            local_hgr_name = hgr_file_name * ".$i"
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
            rm_cmd = `rm -f $local_hgr_name $local_pfile_name`
            run(rm_cmd, wait=true)
        end
        ~, best_cut_idx = findmin(cutsizes)
        rm_cmd = `rm -f $hgr_file_name`
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
    if !config_verbose
        run(`rm -f $log_file`)
    end
end

function ilp_part(; kwargs...)
    try
        log_file = "$config_tmpDir/ilp_part_log_$(time())"
        ilp_string = "$EXTEND_PATH_COMMAND ilp_part $(kwargs[:hgr_file_name]) $config_k $(kwargs[:ub_factor]) >> $log_file 2>&1"
        log_cmd = "echo \"$ilp_string\" > $log_file"
        run(`sh -c $log_cmd`, wait=true)
        ilp_command = `sh -c $ilp_string`
        print("running ilp...")
        run(ilp_command, wait = true)
        if !config_verbose
            run(`rm -f $log_file`)
        end
        print(" - success\n")
        return true
    catch e
        print(" - fail\n")
        return false
    end
end
