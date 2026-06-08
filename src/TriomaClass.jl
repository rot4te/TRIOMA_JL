"""
    TriomaClass

Abstract base type for all TRIOMA types.
Provides `inspect` and `update_attribute!` functionality analogous to
the Python TriomaClass.
"""
abstract type TriomaClass end

"""
    inspect(obj::TriomaClass; variable_name=nothing, indent=0)

Print all fields of `obj` (or only the field matching `variable_name`).
Recurses into nested `TriomaClass` instances.
"""
function inspect(obj::TriomaClass; variable_name::Union{String,Nothing}=nothing, indent::Int=0)
    pad = "    " ^ indent
    for fname in fieldnames(typeof(obj))
        val = getfield(obj, fname)
        name_str = string(fname)
        if variable_name === nothing || lowercase(name_str) == lowercase(variable_name)
            if val isa TriomaClass
                println("$(pad)$(name_str) is a $(typeof(val)), printing its variables:")
                inspect(val; variable_name=variable_name, indent=indent+1)
            else
                println("$(pad)$(name_str): $(val)")
            end
        end
    end
end

"""
    update_attribute!(obj::TriomaClass, attr_name::Symbol, new_value)

Set the field `attr_name` of `obj` (or a nested `TriomaClass` field) to
`new_value`. Mirrors the Python `update_attribute` behaviour, including
propagation of `n_pipes` to nested objects.

Throws `ArgumentError` if the attribute is not found anywhere.
"""
function update_attribute!(obj::TriomaClass, attr_name::Symbol, new_value)
    if attr_name in fieldnames(typeof(obj))
        setfield!(obj, attr_name, new_value)
        # Propagate n_pipes to nested TriomaClass children
        if attr_name === :n_pipes
            for fname in fieldnames(typeof(obj))
                child = getfield(obj, fname)
                if child isa TriomaClass && :n_pipes in fieldnames(typeof(child))
                    setfield!(child, :n_pipes, new_value)
                end
            end
        end
        return
    else
        # Recurse into nested TriomaClass fields
        for fname in fieldnames(typeof(obj))
            child = getfield(obj, fname)
            if child isa TriomaClass && attr_name in fieldnames(typeof(child))
                setfield!(child, attr_name, new_value)
                return
            end
        end
    end
    throw(ArgumentError("'$(attr_name)' is not an attribute of $(typeof(obj))"))
end
