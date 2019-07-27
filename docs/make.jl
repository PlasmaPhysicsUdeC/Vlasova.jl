using Documenter, Vlasova

# If format not defined, make HTML pages
(@isdefined format) ? nothing : (format = Documenter.HTML() )

makedocs(format = format,
         modules = [Vlasova],
         sitename = "Vlasova.jl",
         authors = "Jorge Gidi",
         repo = "https://gitlab.com/jgidi/Vlasova.jl/blob/{commit}{path}#{line}",
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
