using Documenter, Vlasova

# If format is not defined, make HTML pages
(@isdefined format) ? nothing : (format = Documenter.HTML() )

makedocs(format = format,
         modules = [Vlasova],
         sitename = "Vlasova.jl",
         authors = "Jorge Gidi",
         repo = "https://gitlab.com/jgidi/Vlasova.jl/blob/{commit}{path}#{line}",
         pages = [
             "Home" => "index.md", # The name index.md is required by Documenter
             "Basic concepts" => "basic_concepts.md",
             "Integrating a plasma" => "integrating_a_plasma.md",
             "Examples" => "examples.md",
             "API" => [
                 "Exported" => "API/exported.md",
                 "Not exported" => "API/not_exported.md"
             ]
         ]
         )
