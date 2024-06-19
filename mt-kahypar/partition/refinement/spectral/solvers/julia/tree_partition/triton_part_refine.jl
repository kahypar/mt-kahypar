function triton_part_refine(hypergraph_file::String, partition::AbstractArray{Int}, num_parts::Int, ub_factor::Int, seed::Int, identifier::String)
    partition_file = "$(config_tmpDir)/triton_part_input_$identifier.part.2" 
    write_partition(partition, partition_file)

    refined_partition = partition
    if triton_part_refine(hypergraph_file, partition_file, num_parts, ub_factor, seed, identifier)
        refined_partition = read_hint_file(partition_file)
    end
    
    if !config_verbose
        rm = "rm -f $partition_file"
        run(`sh -c $rm`, wait = true)
    end

    return refined_partition
end
    
function triton_part_refine(hypergraph_file::String, partition_file::String, num_parts::Int, ub_factor::Int, seed::Int, identifier::String)
    tcl_line = "triton_part_refine" *
            " -hypergraph_file $hypergraph_file" *
            " -partition_file $partition_file" *
            " -num_parts $num_parts" *
            " -balance_constraint $ub_factor" *
            " -seed $(abs(seed))"

    tcl_file = config_tmpDir * "/run_triton_part_refiner.$identifier.tcl"
    log_file = config_tmpDir * "/$identifier.log"
    sh_file = config_tmpDir * "/run_refiner_$identifier.sh"

    f = open(tcl_file, "w")
    println(f, tcl_line) 
    println(f, "exit")
    close(f)

    f = open(sh_file, "w")
    println(f, EXTEND_PATH_COMMAND)
    println(f, "openroad $tcl_file > $log_file")
    close(f)

    cmd = "chmod 777 " * sh_file 
    run(`sh -c $cmd`, wait=true)
    cmd = "chmod 777 " * tcl_file
    run(`sh -c $cmd`, wait = true)

    success = true
    try
        t = @elapsed run(`$sh_file`, wait=true)
        inform("triton_part took $(t)s")
    catch e
        inform("triton_part failed")
        success = false
    end

    if !config_verbose
        rm = "rm -f $tcl_file $sh_file $log_file"
        run(`sh -c $rm`, wait = true)
    end

    return success
end

triton_part_refine(h, p, n, u, s, id::Int) = triton_part_refine(h, p, n, u, s, split(h, "/")[end] * string(id))