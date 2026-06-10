module TriomaTypes

export inspect, update_attribute!

# True for any mutable struct — i.e. nested TRIOMA sub-objects.
# Primitive types (Float64, Bool, Nothing, String) and arrays return false.
_is_nested(val) = isstructtype(typeof(val)) && ismutabletype(typeof(val))

"""
    inspect(obj; variable_name=nothing, indent=0, io=stdout)

Recursively print all fields of `obj`. If `variable_name` is given, only
print the field with that name (case-insensitive). Nested mutable structs
are expanded in-place with additional indentation.
"""
function inspect(obj; variable_name::Union{String,Nothing}=nothing,
                 indent::Int=0, io::IO=stdout)
    pad = "    " ^ indent
    for fname in fieldnames(typeof(obj))
        val      = getfield(obj, fname)
        name_str = string(fname)
        if variable_name === nothing || lowercase(name_str) == lowercase(variable_name)
            if _is_nested(val)
                println(io, "$(pad)$(name_str) is a $(typeof(val)), printing its variables:")
                inspect(val; variable_name=variable_name, indent=indent+1, io=io)
            else
                println(io, "$(pad)$(name_str): $(val)")
            end
        end
    end
end

"""
    update_attribute!(obj, attr_name, new_value)

Set the field named `attr_name` on `obj` to `new_value`. If the field is not
directly on `obj`, the function recurses into nested mutable-struct children.
Setting `:n_pipes` propagates the value to all nested children that carry it.
Throws `ArgumentError` if no matching field is found anywhere.
"""
function update_attribute!(obj, attr_name, new_value)
    attr = attr_name isa Symbol ? attr_name : Symbol(attr_name)

    if attr in fieldnames(typeof(obj))
        current = getfield(obj, attr)
        if current === nothing
            # Prefer updating a non-nothing child rather than the nothing slot
            for fname in fieldnames(typeof(obj))
                child = getfield(obj, fname)
                if _is_nested(child)
                    try
                        update_attribute!(child, attr, new_value)
                        return
                    catch e
                        e isa ArgumentError || rethrow()
                    end
                end
            end
        end
        setfield!(obj, attr, new_value)
        if attr === :n_pipes
            for fname in fieldnames(typeof(obj))
                child = getfield(obj, fname)
                if _is_nested(child) && :n_pipes in fieldnames(typeof(child))
                    setfield!(child, :n_pipes, new_value)
                end
            end
        end
        return
    end

    # Attribute not on obj directly — recurse into mutable-struct children
    for fname in fieldnames(typeof(obj))
        child = getfield(obj, fname)
        if _is_nested(child)
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
