export vlasova_multithread

"""
    Allow Vlasova to use multithreading on the fourier transforms ( through FFTW )
"""
function vlasova_multithread(nthreads::Integer)
    @assert (nthreads > 0) "Number of threads can't be negative nor zero"

    Vlasova.FFTW.set_num_threads( nthreads )
    
    return 0;
end

