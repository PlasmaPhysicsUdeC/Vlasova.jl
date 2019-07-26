println("\nThis file will add the package Vlasova to your julia installation")
println("Since Vlasova.jl is a private repository, it will ask for your GitLab account and you are required to have been given access to it\n")

using Pkg
pkg"add https://gitlab.com/jgidi/Vlasova.jl"

using Vlasova
