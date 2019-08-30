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
