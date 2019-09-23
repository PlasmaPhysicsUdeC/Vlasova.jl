## Installation

Currently, the [`official repository`](https://gitlab.com/jgidi/Vlasova.jl) is private.
If you want to install this package, you need a gitlab account (you can connect to gitlab using a github account, too) and permission to access the repository (which you can get by asking to the mail `jorgegidi@udec.cl`).


When you have access to read the repository, you can install Vlasova to julia by typing on the REPL:
```julia
]add https://gitlab.com/jgidi/Vlasova.jl
```

Alternatively, Vlasova may be installed by typing:
```julia
using Pkg; Pkg.add(PackageSpec(url = "https://gitlab.com/jgidi/Vlasova.jl"))
```

Using any of the methods shown above, and since the repository is private, julia will ask the `username` and `password` of the account with access to the repository to download the package.



Whenever you want to use Vlasova, you have to enter julia and type:

```julia
using Vlasova
```
