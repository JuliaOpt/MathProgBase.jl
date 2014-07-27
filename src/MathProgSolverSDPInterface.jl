funcs = [:addsdpvar!,
         :addsdpmatrix!,
         :addsdpconstr!,
         :setsdpobj!,
         :getsdpsolution]
for func in funcs
    @eval $(func)() = throw(MethodError($(func),()))
    eval(Expr(:export, func))
end
