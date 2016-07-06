using MathProgBase

function make_lp()
    m = rand((1:100))
    n = rand((1:100))
    c = ones(n)
    x = rand(n)
    A = sprand(m,n,0.3)
    b = A*x
    lb = zeros(n)
    ub = 1000*zeros(n)
    return m,n,c,A,b,lb,ub,x
end

function test1()
    i=1
    j=1
    tol = 1e-3

    while(i<= 100)
        #c,A,b,lb,ub = find_nice()
        m,n,c,A,b,lb,ub,x = make_lp()
        tolerance = tol * ones(n)

        #answer = linprog(c,A,'=',b)
        #answer.status != :Optimal && error("WHAT 1 and j is $j")
        @show x
        @show c
        @show A
        @show b
        @show i,j
        @show newc,newA,newb,newlb,newub,independentvar,pstack = presolver!(c,A,b,lb,ub)
        ans = Array{Float64,1}()
        if(length(find(independentvar))!=0)
            presol = linprog(newc, newA, '=', newb, GLPKSolverLP(method=:Exact, presolve=false))
            presol.status != :Optimal && error("WHAT 2 and j is $j")
            @show ans = presol.sol
        end
        finalsol = return_postsolved(ans,independentvar,pstack)

        if( ((x - finalsol) .< tolerance) == trues(n) )
            println("PASS!")
            j+=1
        end
        i+=1
    end
    @show i,j
end

test1()
