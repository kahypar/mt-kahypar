function triton_part_refine(hypergraph_file::String, partition::AbstractArray{Int}, num_parts::Int, ub_factor::Int, seed::Int, id::Int)
    identifier = split(hypergraph_file, "/")[end] * string(id)

    partition_file = "$(config_tmpDir)/triton_part_input_$identifier.part.2" 
    write_partition(partition, partition_file)

    tcl_line = "triton_part_refine" *
            " -hypergraph_file $hypergraph_file" *
            " -partition_file $partition_file" *
            " -num_parts $num_parts" *
            " -balance_constraint $ub_factor" *
            " -seed $seed"

    tcl_file = config_tmpDir * "/run_triton_part_refiner.$identifier.tcl"
    log_file = config_tmpDir * "/$identifier.log"
    sh_file = config_tmpDir * "/run_refiner_$identifier.sh"

    f = open(tcl_file, "w")
    println(f, tcl_line) 
    println(f, "exit")
    close(f)

    f = open(sh_file, "w")
    line = "PATH=\$PATH:$(config_installationsPrefix)"
    println(f, line)
    line = "openroad $tcl_file > $log_file"
    println(f, line)
    close(f)

    cmd = "chmod 777 " * sh_file 
    run(`sh -c $cmd`, wait=true)
    cmd = "chmod 777 " * tcl_file
    run(`sh -c $cmd`, wait = true)

    triton_part_refiner_cmd = `$sh_file
                                $hypergraph_file
                                $partition
                                $num_parts
                                $ub_factor
                                $seed`
    run(triton_part_refiner_cmd, wait=true)

    refined_partition = read_hint_file(partition_file)

    rm = "rm $tcl_file $sh_file $partition_file $log_file"
    run(`sh -c $rm`, wait = true)

    return refined_partition
end