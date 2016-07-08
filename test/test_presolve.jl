using MathProgBase
using GLPKMathProgInterface

function make_lp(m::Int, n::Int, s::Float64)
    #m = rand((1:100))
    #n = rand((1:100))
    #m = 1000
    #n = 1000
    #enron data set contained 36,692 nodes and 183,831 nodes at sparsity of 0.00013
    #m = 36692
    #n = 36692
    c = ones(n)
    x = rand((1:1000),n)
    A = sprand(m,n,s)
    b = A*x
    lb = float(copy(x))
    ub = 1000000*ones(n)
    return m,n,c,A,b,lb,ub,x
end

# this is only needed because of the way I have generated sample LP instances. Not significant time cost anyway
function trim(x::Array{Float64,1}, y::Array{Float64,1})
    for i in 1:length(x)
        if(x[i] < y[i])
            x[i] = y[i]
        end
    end
    return x
end

#simplest test, fetch one nice instance and compare.
function correctness_test(in1::Int, in2::Int, in3::Float64)
    i=1
    j=1
    tol = 1e-3

    while(i<= 100)
        println("---------STARTING ITERATION $i---------")
        m,n,c,A,b,lb,ub,x = make_lp(in1,in2,in3)
        tolerance = tol * ones(n)
        @show x
        @show c
        #@show A
        @show b
        @show lb
        @show ub

        newc,newA,newb,newlb,newub,independentvar,pstack = presolver!(c,A,b,lb,ub)

        @show size(newA), typeof(newA)
        @show size(newc)
        @show size(newb)
        @show size(x), typeof(x)

        ans = Array{Float64,1}()
        if(length(find(independentvar))!=0)
            presol = linprog(newc, newA, '=', newb, GLPKSolverLP(presolve=false))
            #presol = linprog(newc, newA, '=', newb,newlb,newub)
            presol.status != :Optimal && error("Input feasible problem with an optimal solution but after presolving solver status is not optimal")
            ans = presol.sol
        end
        finalsol = return_postsolved(ans,independentvar,pstack)
        #@show finalsol, typeof(finalsol)
        finalsol = trim(finalsol,lb)

        if( ((x - finalsol) .< tolerance) == trues(n) )
            println("PASS!")
            j+=1
        else
            @show x - finalsol
            @show A*x
            @show A*finalsol
            error("DIDNT PASS!!!")
        end
        i+=1
    end
    if(i == j)
        println("------Presolve works subject to randomized testing-----")
    end
end

function time_test(in1::Int, in2::Int, in3::Float64)
    tol = 1e-3
    println("---------STARTING TIME PROFILE TEST--------")
    m,n,c,A,b,lb,ub,x = make_lp(in1,in2,in3)
    @show m,n

    tolerance = tol * ones(n)
    println("FIRST TIMING")
    @time begin
        answer = linprog(c,A,'=',b,lb,ub,GLPKSolverLP(presolve=true))
        answer.status != :Optimal && error("Input feasible problem with an optimal solution but solver status is not optimal")
    end
println("NEXT TIMING----------------")

@time begin
        newc,newA,newb,newlb,newub,independentvar,pstack = presolver!(c,A,b,lb,ub)
        @show size(newA)
        #(newA == A) && error("No Presolving Done")
        ans = Array{Float64,1}()
            if(length(find(independentvar))!=0)
            presol = linprog(newc, newA, '=', newb,newlb,newub, GLPKSolverLP(presolve=false))
            #presol = linprog(newc, newA, '=', newb,newlb,newub)
            presol.status != :Optimal && error("Input feasible problem with an optimal solution but after presolving solver status is not optimal")
            ans = presol.sol
            end
        finalsol = return_postsolved(ans,independentvar,pstack)
        finalsol = trim(finalsol,lb)
    end
    @show answer.sol
    @show finalsol
end

function do_tests(correctness::Bool, time::Bool)
    if(correctness)
        #correctness_tests
        correctness_test(1000,1000,0.001)
    end

    if(time)
        # Time-Profile tests
        time_test(1,1,0.3)
        println("AGAIN")
        time_test(10,10,0.3)
        println("AGAIN")
        time_test(100,100,0.01)
        println("AGAIN")
        time_test(1000,1000,0.001)
        println("AGAIN YO")
        time_test(10000,10000,0.0001)
    end
end

do_tests(true,false)

do_tests(false,true)
