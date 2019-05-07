"""
Allows vlasova to use multithreading on the fourier transforms ( through FFTW )
and the space and velocity advections ( through OpenMP )
"""
function vlasova_multithread(; FFTW_NUM_THREADS = 0, OMP_NUM_THREADS = 0)
    
    (FFTW_NUM_THREADS == 0) ? nothing : Vlasova.FFTW.set_num_threads( FFTW_NUM_THREADS )
    (OMP_NUM_THREADS == 0) ?  nothing : ( ENV["OMP_NUM_THREADS"] = OMP_NUM_THREADS )

end
