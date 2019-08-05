"""
An empty struct
"""
struct Zero
end

"""
Take any number of aruments and return an element of type Zero
"""
function get_zero(args...)
    return Zero()
end

# Overload sum
import Base.+

"""
 Sum any number to a Zero element, returning the number
"""
@inline function +(n::Number, z::Zero)
    return n
end

"""
Sum an array of numbers to a Zero element, returning the array
"""
@inline function +(A::Array{T} where T<:Number, z::Zero)
    return A
end

# Overload product
import Base.*

"""
 Multiply any number to a Zero element, returning the Zero element
"""
@inline function *(n::Number, z::Zero)
    return z
end

"""
 Multiply an array of numbers to a Zero element, returning the Zero element
"""
@inline function *(A::Array{T} where T<:Number, z::Zero)
    return z
end

# Make sum and prod with zero associative 
+(z::Zero, a) = +(a, z)
*(z::Zero, a) = *(a, z)

# Any broadcasted operation other than sum will return a Zero
@inline function Base.broadcasted(f, A::Array{T} where T<:Number, z::Zero)
    if f == +
        return A
    else
        return z
    end
end
