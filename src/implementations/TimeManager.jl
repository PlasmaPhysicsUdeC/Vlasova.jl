"""
    Set time counter for a TimeManager object

    For example, if `time_manager` is of type TimeManager
    ```
    time_manager( start = Dates.now() )
    ```
"""
function (tm::TimeManager)(;start::Dates.DateTime)
    tm.start = start
end


function notify_progress(tm::TimeManager, t::Int64; filename = "/")
    accomplished = round( 100(t-1)/(tm.final_iteration - 1) ) # TODO: Int
    elapsed = round(Dates.now() - tm.start, Dates.Second )
    notify("\t$accomplished% accomplished in $elapsed", filename = filename) 
end
