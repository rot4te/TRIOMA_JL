module TriomaTypes

export TriomaClass, inspect, update_attribute!

abstract type TriomaClass end

function inspect(obj::TriomaClass; variable_name::Union{String,Nothing}=nothing,
                 indent::Int=0, io::IO=stdout)
    pad = "    " ^ indent
    for fname in fieldnames(typeof(obj))
        val = getfield(obj, fname)
        name_str = string(fname)
        if variable_name === nothing || lowercase(name_str) == lowercase(variable_name)
            if val isa TriomaClass
                println(io, "$(pad)$(name_str) is a $(typeof(val)), printing its variables:")
                inspect(val; variable_name=variable_name, indent=indent+1, io=io)
            else
                println(io, "$(pad)$(name_str): $(val)")
            end
        end
    end
end

function update_attribute!(obj::TriomaClass, attr_name, new_value)
    attr = attr_name isa Symbol ? attr_name : Symbol(attr_name)

    if attr in fieldnames(typeof(obj))
        current = getfield(obj, attr)
        if current === nothing
            # If the direct field is nothing, prefer updating a child that has it non-nothing
            for fname in fieldnames(typeof(obj))
                child = getfield(obj, fname)
                if child isa TriomaClass
                    try
                        update_attribute!(child, attr, new_value)
                        return
                    catch e
                        e isa ArgumentError || rethrow()
                    end
                end
            end
        end
        # Set on this object (either non-nothing, or no child had it)
        setfield!(obj, attr, new_value)
        if attr === :n_pipes
            for fname in fieldnames(typeof(obj))
                child = getfield(obj, fname)
                if child isa TriomaClass && :n_pipes in fieldnames(typeof(child))
                    setfield!(child, :n_pipes, new_value)
                end
            end
        end
        return
    end

    # attr not a direct field — recurse into TriomaClass children
    for fname in fieldnames(typeof(obj))
        child = getfield(obj, fname)
        if child isa TriomaClass
            try
                update_attribute!(child, attr, new_value)
                return
            catch e
                e isa ArgumentError || rethrow()
            end
        end
    end
    throw(ArgumentError("'$(attr_name)' is not an attribute of $(typeof(obj))"))
end

end # module TriomaTypes
