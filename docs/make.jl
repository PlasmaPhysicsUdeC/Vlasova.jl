using Documenter, Vlasova

(@isdefined format) ? nothing : (format = Documenter.HTML() )

makedocs(format = format,
         sitename="Vlasova.jl",
         authors = "Jorge Gidi",
         repo = "https://gitlab.com/jgidi/Vlasova.jl/{commit}{path}#L{line}",
         modules = [Vlasova],
         pages = [
             "Home" => "home.md",
             "Basic concepts" => "basic_concepts.md",
             " A simple example" => "simple_example.md",
             "API" => [
                 "Exported" => "API/exported.md",
                 "Not exported" => "API/not_exported.md"
             ]
             
         ]
         )
