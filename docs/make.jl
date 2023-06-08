using Documenter, Vlasova

# If format is not defined, make HTML pages
(@isdefined format) ? nothing : (format = Documenter.HTML() )

makedocs(format = format,
         modules = [Vlasova],
         sitename = "Vlasova.jl",
         authors = "Jorge Gidi",
         repo = "github.com/PlasmaPhysicsUdeC/Vlasova.jl.git",
         pages = [
             "Home" => "index.md", # The name index.md is required by Documenter
             "Installation" => "installation.md",
             "Basic concepts" => "basic_concepts.md",
             "Default Behavior" => "default_behavior.md",
             "Tools" => "tools.md",
             #"Examples" => "examples.md",
             "Advanced Usage" => "advanced_usage.md",
             "Not Exported" => "not_exported.md"
         ]
         )
