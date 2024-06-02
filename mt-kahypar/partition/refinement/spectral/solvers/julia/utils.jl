function inform(message::String)
    if (config_verbose)
        print("[julia]: " * message * "\n")
        flush(stdout)
    end
end

function inform(graph_size::Integer, big_graphs::Bool, message::String)
    if (graph_size < config_verbose_limits[1] && !big_graphs) || (graph_size > config_verbose_limits[2] && big_graphs)
        inform(message)
    end
end
