export vlasova_multithread

"""
    Allow Vlasova to use multithreading on the fourier transforms ( through FFTW )
    and the space and velocity advections ( written in Fortran, parallelized with OpenMP )
"""
function vlasova_multithread(; FFTW_NUM_THREADS::Integer = 0, OMP_NUM_THREADS::Integer = 0)
    @assert FFTW_NUM_THREADS >= 0 "Number of threads can't be negative"
    @assert OMP_NUM_THREADS >= 0 "Number of threads can't be negative"
    
    if (FFTW_NUM_THREADS == 0) || (FFTW_NUM_THREADS == 0)
        notify("0 threads provided. No change will be made in the number of threads.")
    end
    
    (FFTW_NUM_THREADS == 0) ? nothing : Vlasova.FFTW.set_num_threads( FFTW_NUM_THREADS )
    (OMP_NUM_THREADS == 0) ?  nothing : ( ENV["OMP_NUM_THREADS"] = OMP_NUM_THREADS )
    
    return 0;
end

"""
    Set the number of threads used in the FFT's and space/velocity advections.
"""
function vlasova_multithread(nthreads::Integer)
    
    vlasova_multithread(;FFTW_NUM_THREADS = nthreads,
                        OMP_NUM_THREADS = nthreads   )
        return 0;
end

