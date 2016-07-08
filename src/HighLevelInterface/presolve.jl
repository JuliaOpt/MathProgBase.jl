using MathProgBase
using GLPKMathProgInterface

#export presolver! (?)

# generates the unique key from row,col index for creating the dictionary
function rc(x::Int, y::Int, M::Int)
    return (x-1)*M + y
end

function roughly(x::Float64, y::Float64)
    if(abs(x-y) < 1e-3)
        return true
    else
        return false
    end
end

# A stack is needed to record information on redundancies removed
abstract PresolveStack

`A LinearDependency element represents -
x_index = val + \sum x[vec1[i]] * vec2[i]
if variable is fixed to a value then vec1,vec2 are null.
if variable is fixed to a linear dependency then
vec1 contains the indices of the elements in the constraint row that was removed.
vec2 contains the corresponding ratio of A matrix values.`
type LinearDependency <: PresolveStack
    index :: Int
    vec1 :: Vector{Int}
    vec2 :: Vector{Float64}
    value :: Float64

    function LinearDependency(ind::Int, val::Number)
        vec1 = Array{Int,1}()
        vec2 = Array{Float64,1}()
        new(ind,vec1,vec2,val)
    end
end

# function that will add the LinearDependency element to the stack.
function add_to_stack!(l::LinearDependency, independentvar::BitArray{1}, pstack::Array{PresolveStack,1})
    if (length(l.vec1) != length(l.vec2))
        error("vector1 size not equal to vector 2 size for LD element")
    end
    independentvar[l.index] = false # the variable at this index is not independent anymore
    push!(pstack,l)
end

# function that will post solve one LinearDependency element.
function post_solve!(post_solvedX::Array{Float64,1}, l::LinearDependency)
    post_solvedX[l.index] = l.value

    for i in 1:length(l.vec1)
        post_solvedX[l.index] += l.vec2[i]*post_solvedX[l.vec1[i]]
    end
end

# the x that is fed in has to be the solution obtained from the Solver.
function return_postsolved(x::Array{Float64,1}, independentvar::BitArray{1}, pstack :: Array{PresolveStack,1})
    postsolvedX = zeros(length(independentvar))
    newcols = find(independentvar)

    for i in 1:length(newcols)
        postsolvedX[newcols[i]] = x[i]
    end

    for i in reverse(collect(1:length(pstack)))
        post_solve!(postsolvedX,pstack[i])
    end
    return postsolvedX
end

# Type that will hold the details about the presolve problem we are constructing. The dictionary is made here.
type Presolve_Problem
    currentc :: Array{Float64,1}
    dictA :: Dict{Float64,Float64}
    currentb :: Array{Float64,1}
    currentlb :: Array{Float64,1}
    currentub :: Array{Float64,1}
    #rowindices :: Array{Array{Int64,1},1}
    #colindices :: Array{Array{Int64,1},1}
    independentvar :: BitArray{1}
    activeconstr :: BitArray{1}
    pstack :: Array{PresolveStack,1}
    originalm :: Int64
    originaln :: Int64
    bitmat :: BitArray{2}
    rowcounter :: Array{Float64,1}
    colcounter :: Array{Float64,1}

    #   Constructor that creates a Presolve_Problem
    function Presolve_Problem(c::Array{Float64,1},A::SparseMatrixCSC{Float64,Int64},b::Array{Float64,1},lb::Array{Float64,1},ub::Array{Float64,1})
    #   println("-----------INSIDE PRESOLVE CONSTRUCTOR------------")
        originalm,originaln = size(A)

    #   checks to ensure input problem is valid.
        originalm != length(b) && error("Wrong size of b wrt A")
        originaln != length(lb) && error("Wrong size of lb wrt A")
        originaln != length(ub) && error("Wrong size of ub wrt A")
        originaln != length(c) && error("Wrong size of c wrt A")

    #   SETUP
        independentvar = trues(originaln)
        activeconstr = trues(originalm)
        pstack = Array{PresolveStack,1}()
        dictA = Dict{Float64,Float64}()
        rowcounter = zeros(originalm)
        colcounter = zeros(originaln)
        #rowindices = Array{Int,1}[Int[] for i=1:originalm]
        #colindices = Array{Int,1}[Int[] for i=1:originaln]
        bitmat = falses(originalm,originaln)

    #   Iterating through the non-zeros of sparse matrix A to construct the dictionary
        Arows = rowvals(A)
        Avals = nonzeros(A)
        for j = 1:originaln
            for i in nzrange(A,j)
                r = Arows[i]
                v = Avals[i]
                rcval = rc(r,j,originaln)
                dictA[rcval] = v
                #push!(rowindices[r],j)
                #push!(colindices[j],r)
                rowcounter[r] += 1
                colcounter[j] += 1
                bitmat[r,j] = true
            end
        end

        #new(c,dictA,b,lb,ub,rowindices,colindices,independentvar,activeconstr,pstack,originalm,originaln,bitmat)
        new(c,dictA,b,lb,ub,independentvar,activeconstr,pstack,originalm,originaln,bitmat,rowcounter,colcounter)
    end
end

function presolver!(c::Array{Float64,1}, A::SparseMatrixCSC{Float64,Int64}, b::Array{Float64,1}, lb::Array{Float64,1}, ub::Array{Float64,1})
    #println("Making presolve")
    p = Presolve_Problem(c::Array{Float64,1}, A::SparseMatrixCSC{Float64,Int64}, b::Array{Float64,1}, lb::Array{Float64,1}, ub::Array{Float64,1})

    counter = 0
    while(counter >= 0)
        # println("Empty row call")
        # detects empty rows, throws infeasibilty error or marks it for removal
        er = empty_rows!(p::Presolve_Problem)
        #@show er
        # println("Empty col call")
        # detects emtpy col, throws unbounded error or marks it for removal
        ec = empty_cols!(p::Presolve_Problem)
        #@show ec
        # println("fixed variables call")
        #detects fixed variables, throws infeasiblity error or marks it for removal
        fv = fixed_variables!(p::Presolve_Problem)
        #@show fv
        #println("remove fv")
        # removes all fixed variable from the current copy of data.
        #remfv = remove_fixed!(p::Presolve_Problem)
        #@show remfv
        #println("singleton row call")
        # detects singleton rows
        sr = singleton_rows!(p::Presolve_Problem)
        #sr = false
        #@show sr
        #println("to next iteration")
        if(!(er||ec||fv||sr))
            break
        end
        counter+=1
    end

    #println("trying make new")
    #@time newc,newA,newb,newlb,newub = make_new(p::Presolve_Problem)
    c,A,b,lb,ub = make_new(p::Presolve_Problem)
    return c,A,b,lb,ub,p.independentvar,p.pstack
    #return newc,newA,newb,newlb,newub,p.independentvar,p.pstack
end

    # detecting and removing empty rows.
function empty_rows!(p::Presolve_Problem)
    empty_row = false
    for i in 1:p.originalm
        if(p.activeconstr[i] == false)
            continue
        end
        #if(length(p.rowindices[i])==0)
        if(p.rowcounter[i] == 0)
            #println("Detected Empty Row at $i")
            !(roughly(p.currentb[i],0.0)) && error("Empty Row Infeasibility at row $i and b[i] is - $(p.currentb[i])")
            p.activeconstr[i] = false
            empty_row = true
		end
    end
    #println("Exiting Empty Row and returning $empty_row")
    return empty_row
end

    # detecting and removing empty cols.
function empty_cols!(p::Presolve_Problem)
    empty_col = false
    for j in 1:p.originaln
        if(p.independentvar[j] == false)
            continue
        end
        #if(length(p.colindices[j])==0)
        if(p.colcounter[j] == 0)
                #println("Detected Empty col / Unrestricted variable at $j")
            (p.currentc[j] != 0)*(p.currentlb[j] == -Inf64) && error("Problem is unbounded.")
            p.independentvar[j] = false
            empty_col = true
        end
    end
    #println("Exiting Empty col and returning $empty_col")
    return empty_col
end

    # detecting an infeasible or fixed variable
function fixed_variables!(p::Presolve_Problem)
    detect_fixed = false
    for j in 1:p.originaln
        if(p.independentvar[j] == false)
            continue
        end
        (p.currentlb[j] > p.currentub[j]) && error("Infeasible bounds at $j")
        if(p.currentlb[j] == p.currentub[j])
            add_to_stack!(LinearDependency(j,p.currentlb[j]),p.independentvar,p.pstack)
            detect_fixed = true
        end
    end
    #println("Exiting fixed variable and returning $detect_fixed")
    return detect_fixed
end

    # Removing Fixed Variables
function remove_fixed!(p::Presolve_Problem)
    variable_remove = false
    for j in length(p.independentvar)
        if(p.independentvar[j])
            if((p.currentlb[j] != -Inf64) && p.currentlb[j] == p.currentub[j])
                #println("Found a fixed variable at $j")
                tmp = p.currentlb[j]
                add_to_stack!(LinearDependency(j,tmp),p.independentvar,p.pstack)
                #TODO : substitution into objective function

                #substitution into matrix
                for k in 1:N
                    for i in 1:length(p.rowindices[k])
                        if(p.rowindices[k][i]==j)
                            if length(p.rowindices[k])==1
                                #only x_j variable in this row
                                if((tmp - p.currentb[k]/p.dictA[rc(k,j)]) != 0)
                                    error("Infeasible Problem")
                                end
                                p.currentb[k] -= p.dictA[rc(k,j)]*tmp
                                splice!(p.rowindices[k],i)
                                delete!(p.dictA,rc(k,j));
                                break
                            end
                        end
                        if p.rowindices[k][i] > j
                                break
                        end
                end
            end
            variable_remove = true
            p.colindices[j]=Vector{Int64}[]
            end
        end
    end
    #println("Exiting variable remove and returning $variable_remove")
    return variable_remove
end

    # SINGLETON ROW
function singleton_rows!(p::Presolve_Problem)
    singleton_row = false
    #srows = falses(p.originalm)
    svariable = zeros(p.originaln)

    for i in 1:p.originalm
        if(p.activeconstr[i] == false)
            continue
        end
        #if(length(p.rowindices[i])==1)
        if(p.rowcounter[i] == 1)
            #println("found a row singleton at row $i")
            #j = p.rowindices[i][1]
            j = find(p.bitmat[i,:])[1]
            if(p.independentvar[j] == false)
                continue
            end
            aij = p.dictA[rc(i,j,p.originaln)]
            aij == 0 && error("Unexpected Zero")
            xj = p.currentb[i]/aij
            add_to_stack!(LinearDependency(j,xj),p.independentvar,p.pstack)
            p.activeconstr[i] = false
            singleton_row = true
                    #    srows[i] = true
            svariable[j] = xj
         end
    end

    k = find(p.bitmat)

    for ind in 1:length(k)
        i = k[ind] % p.originalm == 0 ? p.originalm : k[ind]%p.originalm
        j = Int((k[ind] - i)/p.originalm + 1)
        key = rc(i,j,p.originaln)
        if(svariable[j]!=0)
            p.currentb[i] -= svariable[j]*p.dictA[key]
            delete!(p.dictA,key)
            p.bitmat[i,j] = false
            p.rowcounter[i] -= 1
            #splice!(p.rowindices[i],j)
        end
    end
    return singleton_row
end

    # to remove singleton columns
function remove_singleton_cols(p::Presolve_Problem)
    singleton_col = false
    for j in 1:originaln
        if(length(colindices[j])==1)
            nnzrow = colindices[1]
            if(length(rowindices[nnzrow])==1)
                #variable is fixed, will be removed next iteration
                continue
            end
            aij = dictA[rc(nnzrow,j)]
            cj = c[j]
            lbj = lb[j]
            ubj = ub[j]
            is_LB_Unbounded = is_lb_Unbounded(lbj)
            is_UB_Unbounded = is_ub_Unbounded(ubj)

            if(is_LB_Unbounded || is_UB_Unbounded)
                if(is_LB_Unbounded)
                    if(is_UB_Unbounded)
                        zlb[j] = 0
                        zub[j] = 0
                        ylb[nnzrow] = cj/aij
                        yub[nnzrow] = cj/aij
                    elseif aij>0
                        zub[j] = 0
                        ylb[nnzrow] = cij/aij
                    else
                        zub[j] = 0
                        yub[nnzrow] = cij/aij
                    end
                else
                    if(is_UB_Unbounded)
                        if(aij>0)
                            zlb[j]=0
                            yub[nnzrow]=cj/aij
                        else
                            zlb[j]=0
                            ylb[nnzrow]=cj/aij
                        end
                    end
                end

                if(is_LB_Unbounded && is_UB_Unbounded)
                    #free column singletons
                    cnt = length(rowindices[nnzrow])-1
                    vec1 = vec2 = Vector{Int}()
                    #for nonzero rows.
                    for i in 1:currentm
                        #continue from here

                    end


                end
            end
        end
    end
    return singleton_col
end

function check_progress(currentc,dictA,curerntb,currentlb,currentub,rowindices,colindices)
    #checking rows
    #TODO .. remove rowindices,colindices usage
    for i in 1:N
        for nnz in 1:length(rowindices[i])
            dictA[rc(i,nnz)] == 0 && error("Zero erro in row")
        end
    end

    # checking columns
    for j in 1:M
        for nnz in 1:length(colindices[j])
            dictA[rc(nnz,j)] == 0 && error("Zero error in col")
        end
    end

    # checking Aij
    for key in keys(dictA)
        j = key%M == 0 ? 4 : key%M
        i = Int64((key -j)/M + 1)
        if(dictA[key] != 0)
            if(!in(rowindices[i]),j)
                error("$j not in row $i")
            end
            if(!in(colindices[j]),i)
                error("$i not in col $j")
            end
        end
    end

    # think of other helpful tests
end

    # to make c,A,sense,b,l,u
function make_new(p::Presolve_Problem)
    newc = Array{Float64,1}()
    for j in 1:p.originaln
        if(p.independentvar[j])
            push!(newc,p.currentc[j])
        end
    end

    newb = Array{Float64,1}()
    for j in 1:p.originalm
        if(p.activeconstr[j])
            push!(newb,p.currentb[j])
        end
    end

    newrows = find(p.activeconstr)
    newcols = find(p.independentvar)

    # constructing a sparse matrix
    #@show len = length(p.dictA)
    I = Array{Int64,1}()
    J = Array{Int64,1}()
    Val = []

    for rowiter in 1:length(newrows)
        for coliter in 1:length(newcols)
            if(haskey(p.dictA,rc(newrows[rowiter],newcols[coliter],p.originaln)))
                push!(I,rowiter)
                push!(J,coliter)
                push!(Val,p.dictA[rc(newrows[rowiter],newcols[coliter],p.originaln)])
            end
        end
    end
    newA = sparse(I,J,Val,length(newrows),length(newcols))

    newlb = zeros(length(newcols))
    for coliter in 1:length(newcols)
        newlb[coliter] = p.currentlb[newcols[coliter]]
    end

    newub = zeros(length(newcols))
    for coliter in 1:length(newcols)
        newub[coliter] = p.currentub[newcols[coliter]]
    end
    return newc,newA,newb,newlb,newub
end
