using Documenter, Vlasova


(@isdefined format) ? nothing : (format = Documenter.HTML() )

makedocs(format = format,
         sitename="Vlasova.jl",
         authors = "Jorge Gidi",
         modules = [Vlasova],
         pages = [
             "Manual" => [
                 "Guide" => "guide.md"
             ],
             
             "Examples" => "examples.md",
             "API" => [
                 "Structs" => "structs.md"
                 "Functions" => "functions.md"
             ]
             
         ]
         )
