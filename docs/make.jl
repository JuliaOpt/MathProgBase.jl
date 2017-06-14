using Documenter, MathProgBase

makedocs(
    format = :html,
    sitename = "MathProgBase",
    pages = [
        "Introduction" => "index.md",
        "Solver Interface" => "solverinterface.md"
    ]
)

deploydocs(
    repo   = "github.com/JuliaOpt/MathProgBase.jl.git",
    target = "build",
    osname = "linux",
    julia  = "0.6",
    deps   = nothing,
    make   = nothing
)
