using Documenter, MathProgBase

makedocs(
    format = :html,
    sitename = "MathProgBase",
    pages = [
        "Introduction" => "index.md",
        "High-level Interfaces" => "highlevel.md",
        "Solver Interface" => [
             "Basics" => "basics.md",
             "Variables" => "variables.md",
             "Objectives" => "objectives.md",
             "Constraints" => "constraints.md",
             "Sets" => "sets.md",
             "Attributes" => "attributes.md",
             "Status Codes" => "statuscodes.md",
             "Duals" => "duals.md",
             "NLP" => "nlp.md"],
        "Choosing Solver" => "choosingsolver.md"
    ]
)

deploydocs(
    repo   = "github.com/JuliaOpt/MathProgBase.jl.git",
    target = "build",
    osname = "linux",
    julia  = "0.6",
    deps   = nothing,
    make   = nothing,
    latest = "break_everything"
)
