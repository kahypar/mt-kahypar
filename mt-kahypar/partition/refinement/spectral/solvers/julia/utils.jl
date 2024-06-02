function inform(message::String)
    if (config_verbose)
        print("[julia]: " * message * "\n")
        flush(stdout)
    end
end

function inform(graph_size::Integer, big_graphs::Bool, message::String)
    inform(graph_size, big_graphs, () -> message)
end

function inform(graph_size::Integer, big_graphs::Bool, message_getter::Function)
    if (graph_size < config_verbose_limits[1] && !big_graphs) || (graph_size > config_verbose_limits[2] && big_graphs)
        inform(message_getter())
    end
end

macro print_backtrace()
    quote 
        inform(sprint((io, v) -> show(io, "text/plain", v), stacktrace(catch_backtrace())))
    end
end
