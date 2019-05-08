using Documenter, Vlasova

makedocs(sitename="Vlasova.jl",
         modules = [Vlasova],
         pages = [
             "index.md",
             "structs" => "structs.md",
             "functions" => [
                 "Exported functions" => "functions_public.md"
                 "Private functions" => "functions_private.md"
                 
             ]
         ])
