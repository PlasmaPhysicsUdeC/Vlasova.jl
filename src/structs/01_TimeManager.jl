import Dates

mutable struct TimeManager
    final_time::Float64
    final_iteration::Int64
    dt::Float64
    start::Dates.DateTime
    
    TimeManager(final_time,
                dt ) = new(final_time,
                           floor(Int, final_time / dt + 1),
                           dt,
                           Dates.now() )
end
