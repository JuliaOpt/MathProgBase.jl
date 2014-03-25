using Base.Test
using MathProgBase
using MathProgBase.MathProgSolverInterface

function sdptest(solver=MathProgBase.defaultSDPsolver)
    println("Testing SDP interface with solver ", string(typeof(solver)))
    
    # min [2 1 0;
    #      1 2 1; ◦ X + y - z
    #      0 1 2]
    # s.t. I◦X + y == 1
    #      J◦X + z >= 1/2
    #      X(3,3) - psd
    #      0 <= y,z <= 1

    m = model(solver)
    id1 = addsdpvar!(m, 3)
    addvar!(m, 0.0, 1.0,  1.0)
    addvar!(m, 0.0, 1.0, -1.0)
    id2 = addsdpmatrix!(m, [2 1 0;1 2 1;0 1 2])
    setsdpobj!(m, [id1], [id2])
    id3 = addsdpmatrix!(m, eye(3,3))
    addsdpconstr!(m, [id1], [id3], [1], [1.0], 1.0, 1.0)
    id4 = addsdpmatrix!(m, ones(3,3))
    addsdpconstr!(m, [id1], [id4], [2], [1.0], 0.5, Inf)
    optimize!(m)
    stat = status(m)

    @test stat == :Optimal
    @test_approx_eq_eps getobjval(m) -0.4142 1e-3
    @test_approx_eq_eps norm(eig(getsdpsolution(m,id1))[1] .- [0.0,0.0,1.0]) 0.0 1e-3
    @test_approx_eq_eps norm(getsolution(m) - [0,1]) 0.0 1e-3
    println("Done")

end
