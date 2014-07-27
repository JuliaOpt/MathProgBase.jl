funcs = [:loadconicproblem!, :getconicdual]

for func in funcs
    @eval $(func)() = throw(MethodError($(func),()))
    eval(Expr(:export, func))
end
