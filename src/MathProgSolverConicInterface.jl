export loadconicproblem!,
       getconicdual

for func in [:loadconicproblem, :getconicdual]
    @eval $(func)() = throw(MethodError($(func),()))
end
