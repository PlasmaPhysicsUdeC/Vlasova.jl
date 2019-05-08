#TODO: Use BinDeps

# Compile fortran libraries to deps/src/usr/lib
function compile_fortran(filename::String, prefix::String, postfix::String)
    run(`gfortran -Wall -shared -fPIC -O3 -fopenmp $prefix$filename.f90 -o $postfix$filename.so`)
end

# All files to be compiled
source_code = ["advections", "velocity_filter"]

# Paths
current_path = dirname(@__FILE__)
prefix = "$current_path/src/"
postfix = "$current_path/usr/lib/"

# Compile
mkpath(prefix)
mkpath(postfix)
for file in source_code
    compile_fortran(file, prefix, postfix)
end
