export addsdpvar!,
       addsdpmatrix!,
       addsdpconstr!,
       setsdpobj!,
       getsdpsolution

for func in [:addsdpvar!,
             :addsdpmatrix!,
             :addsdpconstr!,
             :setsdpobj!,
             :getsdpsolution]
    @eval $(func)() = throw(MethodError($(func),()))
end
