"""
An empty struct that behaves like a zero.
"""
struct Zero
end

"""
Take any number of aruments and return an element of type Zero
"""
get_zero(args...) = Zero()


# Overload sum
import Base.+

"""
 Sum any number to a Zero element, returning the number
"""
@inline +(n::Number, z::Zero) = n

"""
Sum an array of numbers to a Zero element, returning the array
"""
@inline +(A::Array{T} where T<:Number, z::Zero) = A

# Overload product
import Base.*

"""
 Multiply any number to a Zero element, returning the Zero element
"""
@inline *(n::Number, z::Zero) = z

"""
 Multiply an array of numbers to a Zero element, returning the Zero element
"""
@inline *(A::Array{T} where T<:Number, z::Zero) = z

# Make sum and prod with zero associative 
+(z::Zero, a) = +(a, z)
*(z::Zero, a) = *(a, z)

# A broadcasted sum of an array and a zero will return the array
@inline Base.Broadcast.broadcasted(+, A::Array{T} where T<:Number, z::Zero) = A

# Show "0::Zero" to make it different from 0.
Base.show(io::IO, z::Zero) = print(io, "0::Zero")
