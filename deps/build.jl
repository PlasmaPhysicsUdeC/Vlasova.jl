#TODO: Use BinDeps!
# using BinDeps

# Compile fortran libraries to deps/src/usr/lib
function compile_fortran(filename::String)
    current_path = dirname(@__FILE__)
    prefix = "$current_path/src/"
    postfix = "$current_path/usr/lib/"
    run(`gfortran -Wall -shared -fPIC -O3 -fopenmp $prefix$filename.f90 -o $postfix$filename.so`)
end

# All files to be compiled
source_code = ["advections", "velocity_filter"]
for file in source_code
    compile_fortran(file)
end

# @BinDeps.setup

# for file in source_code
#     # Use gfortran to compile source code
#     compile_fortran(file)

#     lib = library_dependency("$file.so")
#     # Make library generated available to julia
#     provides(Binaries,
#              URI("file:$(dirname(@__FILE__))/usr/lib/$file.so"),
#              lib,
#              os = :Unix)
# end

# @BinDeps.install
