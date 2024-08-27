using DelimitedFiles

# logging functions

function inform(message::String)
    print("[julia]: " * message[1:min(end,200)] * "\n")
    flush(stdout)
end

# inform depending on instance size being extreme
function inform(graph_size::Integer, big_graphs::Bool, message_getter::Function)
    if (graph_size < config_verbose_limits[1] && !big_graphs) || (graph_size > config_verbose_limits[2] && big_graphs)
        inform(message_getter())
    end
end

function inform(graph_size::Integer, big_graphs::Bool, message::String)
    inform(graph_size, big_graphs, () -> message)
end

inform_dbg = config_verbose ? inform : (x...) -> nothing

macro print_backtrace()
    quote 
        inform(sprint((io, v) -> show(io, "text/plain", v), stacktrace(catch_backtrace())))
    end
end

function pretty_print(A)
    str = IOBuffer()
    show(IOContext(str, :compact => false), "text/plain", A)
    return String(take!(str))
end

function read_hint_file(file_name::String)
    line_count = 0
    p = readdlm(file_name, Int)
    partition = p[:,1]
    return partition
end

# extend unix path for local installations
EXTEND_PATH_COMMAND = "PATH=\$PATH:$(config_installationsPrefix)/bin"
