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
             "Installation" => "installation.md",
             "Basic concepts" => "basic_concepts.md",
             "Default Behavior" => "default_behavior.md",
             "Tools" => "tools.md",
             "Advanced Usage" => "advanced_usage.md",
             "Examples" => "examples.md",
             "Not Exported" => "not_exported.md"
         ]
         )
