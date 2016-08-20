workspace()

using MathProgBase
using GLPKMathProgInterface

abstract Presolve_Element

type Presolve_Row <: Presolve_Element
    i :: Int
    b_val :: Float64
    aij :: Union{Presolve_Element,Nothing}
    prev :: Presolve_Row
    next :: Presolve_Row
    is_active :: Bool
    active_prev :: Presolve_Row
    active_next :: Presolve_Row
    function Presolve_Row()
        n = new()
        n.i = -1
        n.b_val = 0.0
        n.is_active = true
        n.aij = nothing
        n.prev = n
        n.next = n
        n.active_prev = n
        n.active_next = n
        n
    end
end

type Presolve_Col <: Presolve_Element
    j :: Int
    c_val :: Float64
    aij :: Union{Presolve_Element,Nothing}
    prev :: Presolve_Col
    next :: Presolve_Col
    is_independent :: Bool
    ind_prev :: Presolve_Col
    ind_next :: Presolve_Col
    function Presolve_Col()
        n = new()
        n.j = -1
        n.c_val = 0.0
        n.is_independent = true
        n.aij = nothing
        n.prev = n
        n.next = n
        n.ind_prev = n
        n.ind_next = n
        n
    end
end

type Presolve_Matrix <: Presolve_Element
    row :: Presolve_Row
    col :: Presolve_Col
    val :: Float64
    row_prev :: Presolve_Matrix
    row_next :: Presolve_Matrix
    col_prev :: Presolve_Matrix
    col_next :: Presolve_Matrix
    function Presolve_Matrix()
        n = new()
        row = Presolve_Row()
        col = Presolve_Col()
        val = 0.0
        row_prev = n
        row_next = n
        col_prev = n
        col_next = n
        n
    end
end

# A stack is needed to record information on redundancies removed
abstract PresolveStack

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

    function LinearDependency(ind::Int, vec1::Vector{Int}, vec2::Vector{Float64}, val::Number)
        new(ind,vec1,vec2,val)
    end

end


# Type that will hold the details about the presolve problem we are constructing. The dictionary is made here.
type Presolve_Problem
    currentc :: Array{Float64,1}
    currentb :: Array{Float64,1}
    currentlb :: Array{Float64,1}
    currentub :: Array{Float64,1}
    currentylb :: Array{Float64,1}
    currentyub :: Array{Float64,1}
    currentzlb :: Array{Float64,1}
    currentzub :: Array{Float64,1}
    independentvar :: BitArray{1}
    activeconstr :: BitArray{1}
    pstack :: Array{PresolveStack,1}
    originalm :: Int64
    originaln :: Int64
    rowcounter :: Array{Float64,1}
    colcounter :: Array{Float64,1}
    g :: Array{Float64,1}
    h :: Array{Float64,1}

    # new thingies
    dictrow :: Dict{Int64,Presolve_Row}
    dictcol :: Dict{Int64,Presolve_Col}
    dictaij :: Dict{Int64,Presolve_Matrix}

    rowptr :: Presolve_Row
    colptr :: Presolve_Col
    rowque :: Presolve_Row
    colque :: Presolve_Col

    finalm :: Int
    finaln :: Int
    finalrows :: Array{Int,1}
    finalcols :: Array{Int,1}
    #   Constructor that creates a Presolve_Problem
    function Presolve_Problem(verbose::Bool,m::Int,n::Int)
        verbose && println("-----------INSIDE PRESOLVE CONSTRUCTOR------------")
        originalm,originaln = m,n

        #   SETUP
        independentvar = trues(originaln)
        activeconstr = trues(originalm)
        pstack = Array{PresolveStack,1}()
        rowcounter = zeros(originalm)
        colcounter = zeros(originaln)
        g = fill(-Inf,originalm)
        h = fill(+Inf,originalm)
        ylb = fill(-Inf,originalm)
        yub = fill(+Inf,originalm)
        zlb = fill(-Inf,originaln)
        zub = fill(+Inf,originaln)
        c = zeros(originaln)
        b = zeros(originalm)
        lb = zeros(originalm)
        ub = zeros(originalm)

    #   SETUP for new variables
        dictrow = Dict{Int64,Presolve_Row}()
        dictcol = Dict{Int64,Presolve_Col}()
        dictaij = Dict{Int64,Presolve_Matrix}()

        rowptr = Presolve_Row()
        colptr = Presolve_Col()
        rowque = Presolve_Row()
        colque = Presolve_Col()

        finalm = originalm
        finaln = originaln
        finalrows = fill(-1,originalm)
        finalcols = fill(-1,originaln)

        new(c,b,lb,ub,ylb,yub,zlb,zub,independentvar,activeconstr,pstack,originalm,originaln,rowcounter,colcounter,g,h,dictrow,dictcol,dictaij,rowptr,colptr,rowque,colque,finalm,finaln,finalrows,finalcols)
    end
end

function add_row!(p::Presolve_Problem, i::Int, bval::Float64)
    println("Adding row : $(i)")

    row = Presolve_Row()
    row.i = i
    row.b_val = bval
    row.aij = nothing
    row.is_active = false
    if(i!=1)
        row.prev = p.dictrow[i-1]
        row.prev.next = row
    else
        row.prev = row
        p.rowptr = row
    end
    row.next = row
    enque_row!(p,row)
    p.dictrow[i] = row
end

function enque_row!(p::Presolve_Problem, row::Presolve_Row)
    println("Queueing row : $(row.i)")

    if(row.is_active == false)
        row.is_active = true
        row.active_prev = row.prev
        row.active_next = row
        if(p.rowque.i == -1)
            p.rowque = row
        end
        if(row.active_prev != row)
            row.active_prev.active_next = row
        end
    end
end

function deque_row!(p::Presolve_Problem, row::Presolve_Row)
    if(row.is_active == true)
        row.is_active = false
        if(row.active_prev == row)
            p.rowque = row.active_next
        else
            row.active_prev.active_next = row.active_next
        end
        if(row.active_next != row)
            row.active_next.active_prev = row.active_prev
        end
    end
end

function add_col!(p::Presolve_Problem, j::Int, cval::Float64)
    println("Adding col : $(j)")

    col = Presolve_Col()
    col.j = j
    col.c_val = cval
    col.aij = nothing
    col.is_independent = false
    if(j!=1)
        col.prev = p.dictcol[j-1]
        col.prev.next = col
    else
        col.prev = col
        p.colptr = col
    end
    col.next = col
    enque_col!(p,col)
    p.dictcol[j] = col
end

function enque_col!(p::Presolve_Problem, col::Presolve_Col)
    println("Queueing col : $(col.j)")

    if(col.is_independent == false)
        col.ind_prev = col.prev
        col.is_independent = true
        col.ind_next = col
        if(p.colque.j == -1)
            p.colque = col
        end
        if(col.ind_prev != col)
            col.ind_prev.ind_next = col
        end
    end
end

function deque_col!(p::Presolve_Problem, col::Presolve_Col)
    if(col.is_independent == true)
        col.is_independent = false
        if(col.ind_prev == col)
            p.colque = col.ind_next
        else
            col.ind_prev.ind_next = col.ind_next
        end
        if(col.ind_next != col)
            col.ind_next.ind_prev = col.ind_prev
        end
    end
end

function add_aij_normal!(p::Presolve_Problem, row_id::Int, col_id::Int, row_prev_id::Int, val::Float64)
    println("Adding mat element : $(row_id),$(col_id)")

    aij = Presolve_Matrix()
    aij.row = p.dictrow[row_id]
    aij.col = p.dictcol[col_id]
    aij.val = val
    aij.row_prev = aij
    aij.row_next = aij
    aij.col_next = aij
    aij.col_prev = aij

    if(row_prev_id != -1)
        prev_key = rc(row_prev_id,col_id,p.originaln)
        p.dictaij[prev_key].col_next = aij
        aij.col_prev = p.dictaij[prev_key]
    end
    #can also be done by checking if row_prev == -1
    if(p.dictcol[col_id].aij == nothing)
        p.dictcol[col_id].aij = aij
    end

    key = rc(row_id,col_id,p.originaln)
    p.dictaij[key] = aij
end

function add_aij_transpose!(p::Presolve_Problem, row_id::Int, col_id::Int, col_prev_id::Int, val::Float64)
    aij = p.dictaij[rc(row_id,col_id,p.originaln)]

    if(col_prev_id != -1)
        prev_key = rc(row_id,col_prev_id,p.originaln)
        p.dictaij[prev_key].row_next = aij
        aij.row_prev = p.dictaij[prev_key]
    end

    if(p.dictrow[row_id].aij == nothing)
        p.dictrow[row_id].aij = aij
    end
end

function remove_row!(p::Presolve_Problem, row::Presolve_Row)
    deque_row!(p,row)

    while(row.aij != nothing)
        println("INSIDE------ and aij is $(row.aij.row.i),$(row.aij.col.j)")
        tmp = row.aij
        key = rc(tmp.row.i,tmp.col.j,p.originaln)
        #enque_col!(p,tmp.col)
        row.aij = tmp.row_next

        if(tmp.col_prev == tmp)
            if(tmp.col_next == tmp)
                tmp.col.aij = nothing
            else
                tmp.col.aij = tmp.col_next
                tmp.col_next.col_prev = tmp.col_next
            end
        else
            if(tmp.col_next == tmp)
                tmp.col_prev.col_next = tmp.col_prev
            else
                tmp.col_prev.col_next = tmp.col_next
            end
        end

        delete!(p.dictaij,key)

        if(row.aij == tmp)
            row.aij = nothing
        end
        tmp = nothing
    end

    if(row.prev == row)
        if(row.next != row)
            p.rowptr = row.next
            row.next.prev = row.next
        else
            p.rowptr = Presolve_Row()
        end
    else
        if(row.next != row)
            row.prev.next = row.next
            row.next.prev = row.prev
        else
            row.prev.next = row.prev
        end
    end

    delete!(p.dictrow,row.i)
end

function remove_col!(p::Presolve_Problem, col::Presolve_Col)
    println("REMOVE COLUMN CALL")
    deque_col!(p,col)

    while(col.aij != nothing)
        tmp = col.aij
        key = rc(tmp.row.i,tmp.col.j,p.originaln)
        #enque_row!(p,tmp.row)
        col.aij = tmp.col_next

        if(tmp.row_prev == tmp)
            if(tmp.row_next == tmp)
                tmp.row.aij = nothing
            else
        #        println("here 1 ")
                tmp.row.aij = tmp.row_next
                tmp.row_next.row_prev = tmp.row_next
            end
        else
            if(tmp.row_next == tmp)
                tmp.row_prev.row_next = tmp.row_prev
            else
                tmp.row_prev.row_next = tmp.row_next
            end
        end

        delete!(p.dictaij,key)

        if(col.aij == tmp)
            col.aij = nothing
        end
        tmp = nothing
    end

    if(col.prev == col)
        if(col.next != col)
            p.colptr = col.next
            col.next.prev = col.next
        else
            p.colptr = Presolve_Col()
        end
    else
        if(col.next != col)
            col.prev.next = col.next
            col.next.prev = col.prev
        else
            col.prev.next = col.prev
        end
    end
    delete!(p.dictcol,col.j)
end

function make_presolve!(p::Presolve_Problem,c::Array{Float64,1}, A::SparseMatrixCSC{Float64,Int64}, b::Array{Float64,1}, lb::Array{Float64,1}, ub::Array{Float64,1})
#   checks to ensure input problem is valid.
    m,n = size(A)
    p.originalm != m && error("Wrong size of b wrt A")
    p.originaln != n && error("Wrong size of c wrt A")
    p.originaln != length(lb) && error("Wrong size of lb wrt A")
    p.originaln != length(ub) && error("Wrong size of ub wrt A")

    p.currentc = c
    p.currentb = b
    p.currentlb = lb
    p.currentub = ub
    println("Row SETUP ----- ")
    for i in 1:p.originalm
        add_row!(p,i,b[i])
    end

    println("COL SETUP -----")
    for j in 1:p.originaln
        add_col!(p,j,c[j])
    end

    #   Iterating through the non-zeros of sparse matrix A to construct the dictionary
    Arows = rowvals(A)
    println("MAT ELEMENT SETUP -----")
    Avals = nonzeros(A)
    for j = 1:p.originaln
        tmp = -1
        for i in nzrange(A,j)
            r = Arows[i]
            rcval = rc(r,j,p.originaln)
            #dictA[rcval] = Avals[i]
            p.rowcounter[r] += 1
            p.colcounter[j] += 1
            add_aij_normal!(p,r,j,tmp,Avals[i])
            tmp = r
        end
    end

    Arows = rowvals(A')
    Avals = nonzeros(A')
    for i = 1:p.originalm
        tmp = -1
        for c in nzrange(A',i)
            j = Arows[c]
            rcval = rc(i,j,p.originaln)
            add_aij_transpose!(p,i,j,tmp,Avals[c])
            tmp = j
        end
    end
end

function print_info(p::Presolve_Problem)
    println("Row Information--------------------------------------")
    for key in keys(p.dictrow)
        println("-----------")
        row = p.dictrow[key]
        @show row.i
        @show row.b_val
        if(row.aij != nothing)
            @show row.aij.row.i
            @show row.aij.col.j
            @show row.aij.val
        end
        @show row.prev.i
        @show row.next.i
        @show row.is_active
        @show row.active_prev.i
        @show row.active_next.i
        #println("Id : $(row.i)")
        #println("b_val : $row.b_val")
    end

    println("Col Information-------------------------------------------")
    for key in keys(p.dictcol)
        println("-----------")
        col = p.dictcol[key]
        @show col.j
        @show col.c_val
        if(col.aij != nothing)
            @show col.aij.row.i
            @show col.aij.col.j
            @show col.aij.val
        end
        @show col.prev.j
        @show col.next.j
        @show col.is_independent
        @show col.ind_prev.j
        @show col.ind_next.j
        #println("Id : $(row.i)")
        #println("b_val : $row.b_val")
    end

    println("MAT ELEMENT Information---------------------------")
    for key in keys(p.dictaij)
        println("-----------")
        aij = p.dictaij[key]
        @show aij.row.i
        @show aij.col.j
        @show aij.val

        @show aij.row_prev.row.i, aij.row_prev.col.j, aij.row_prev.val
        @show aij.row_next.row.i, aij.row_next.col.j, aij.row_next.val
        @show aij.col_prev.row.i, aij.col_prev.col.j, aij.col_prev.val
        @show aij.col_next.row.i, aij.col_next.col.j, aij.col_next.val
        #println("Id : $(row.i)")
        #println("b_val : $row.b_val")
    end
end

function presolver!(verbose::Bool,c::Array{Float64,1}, A::SparseMatrixCSC{Float64,Int64}, b::Array{Float64,1}, lb::Array{Float64,1}, ub::Array{Float64,1})
    v = verbose
    v && println("Making presolve")
    p = Presolve_Problem(v,length(b),length(c))
    make_presolve!(p,c,A,b,lb,ub)

    println("AFTER MAKE PRESOLVE -------------")
    #print_info(p)
    #remove_col!(p,p.dictcol[2])
    #print_info(p)
    #remove_row!(p,p.dictrow[1])
    #print_info(p)

    println("TIME FOR PRESOLVE ...............................")
    row = Presolve_Row()
    col = Presolve_Col()
    tmp = p.rowque

    while(tmp != nothing)
        row = tmp
        deque_row!(p,row)
        if(row.aij == nothing)
            println("EMPTY ROW FOUND AT $(row.i)")
            empty_row!(p,row,v)
        else
            if(row.aij.row_next == row.aij)
                println("SINGETONE ROW FOUND AT $(row.i)")
                singleton_row!(p,row,v)
            else
                println("happy for now")
                #forcing_constraints!(p,row,v)
            end
        end
        if(tmp.next == tmp)
            tmp = nothing
        else
            tmp = tmp.next
        end
    end

    v && println("trying make new")
    c,A,b,lb,ub = make_new(p::Presolve_Problem,v)
    @show c
    @show A
    @show b
    @show lb
    @show ub
    return c,A,b,lb,ub,p.independentvar,p.pstack
end

    # detecting and removing empty rows.
function empty_row!(p::Presolve_Problem, row::Presolve_Row, v::Bool)
    if(!roughly(row.b_val,0.0))
        error("Empty Row Infeasibility at row $row.i and b[i] is - $(row.b_val)")
    else
        remove_row!(p,row)
        p.activeconstr[row.i] = false
    end

    v && println("Exiting Empty Row")
end

    # SINGLETON ROW
function singleton_row!(p::Presolve_Problem, row::Presolve_Row, v::Bool)
    i = row.aij.row.i
    j = row.aij.col.j
    matval = row.aij.val
    bval = row.aij.row.b_val

    xj = bval/matval
    add_to_stack!(LinearDependency(j,xj),p.independentvar,p.pstack)

    remove_row!(p,row)
    p.activeconstr[row.i] = false
    aij = p.dictcol[j].aij
    while(aij != nothing)
        r = aij.row
        r.b_val -= xj*aij.val
        if(aij.col_next != aij)
            aij = aij.col_next
        else
            aij = nothing
        end
    end
    remove_col!(p,p.dictcol[j])
end

    # to make c,A,sense,b,l,u
function make_new(p::Presolve_Problem, v::Bool)
    n = 0
    newc = Array{Float64,1}()
    col = p.colptr
    if(col.j != -1)
        println("CONSTRUCTING newc")
        while(col != nothing)
            @show col.j
            push!(newc,col.c_val)
            n = n + 1
            p.finalcols[col.j] = n

            if(col.next == col)
                col = nothing
            else
                col = col.next
            end
        end
    end
    @show n

    m = 0
    newb = Array{Float64,1}()
    newlb = Array{Float64,1}()
    newub = Array{Float64,1}()
    row = p.rowptr
    if(row.i != -1)
        while(row != nothing)
            push!(newb,row.b_val)
            push!(newlb,p.currentlb[row.i])
            push!(newub,p.currentub[row.i])
            m = m + 1
            p.finalrows[row.i] = m
            if(row.next == row)
                row = nothing
            else
                row = row.next
            end
        end
    end
    @show m

    I = Array{Int64,1}()
    J = Array{Int64,1}()
    Val = Array{Float64,1}()
    col = p.colptr
    if(col.j != -1)
        while(col != nothing)
            tmp = col.aij
            while(tmp != nothing)
                push!(J,p.finalcols[tmp.col.j])
                push!(I,p.finalrows[tmp.row.i])
                push!(Val,tmp.val)
                if(tmp.row_next == tmp)
                    tmp = nothing
                else
                    tmp = tmp.row_next
                end
            end
            if(col.next == col)
                col = nothing
            else
                col = col.next
            end
        end
    end
    newA = sparse(I,J,Val,m,n)
    return newc,newA,newb,newlb,newub
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

    @show length(l.vec1)
    @show post_solvedX

    for i in 1:length(l.vec1)
        post_solvedX[l.index] += l.vec2[i]*post_solvedX[l.vec1[i]]
        println("made postsolved at $(l.index) to value $(post_solvedX[l.index])")
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
        @show pstack[i]
        post_solve!(postsolvedX,pstack[i])
    end
    return postsolvedX
end


function is_zero(i::Float64)
    if(abs(i-0.0) <= 1e-3)
        return true
    else
        return false
    end
end

function is_equal(a::Array{Float64,1}, b::Array{Float64,1})
    (length(a) != length(b)) && error("trying to determine equality of arrays of different sizes")
    for i in 1:length(a)
        if(is_zero(a[i]-b[i]) == false)
            return false
        end
    end
    return true
end

function is_lb_Unbounded(lb::Float64)
    if(lb == -Inf)
        return true
    else
        return false
    end
end


function is_ub_Unbounded(ub::Float64)
    if(ub == +Inf)
        return true
    else
        return false
    end
end

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
