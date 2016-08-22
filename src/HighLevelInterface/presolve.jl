#workspace()

`
This is an implementation of LP Presolving based on ideas from Andersen & Andersen paper - http://www.turma-aguia.com/davi/doc/Andersen.pdf
Presolving removes redundancies from the original problem given by the user and constructs a smaller equivalent problem
which is then fed to the solver.
The input to the presolver! function is the data fed into the linprog() function of MathProgBase API.
It returns the presolved problem along with a "Presolve" stack which contains information on each redundancy that was removed.
This was the objective of a 2016 Google Summer of Code project under Julia-Opt.
A blog describing the progress and issues can be found here - https://ramcha24.github.io/gsoc
Author : Ramchandran (https://github.com/ramcha24)
Mentor : Madeleiene Udell (https://github.com/madeleineudell/)
All queries and comments are welcome`

using MathProgBase

`
--- Introduction ---
The optimization problem being solved is of the form -
    min c*x
    s.t A*x = b
        lb <= x <= ub

    where x, c, lb, ub are n-dimensional vectors. b is a m-dimensional vector.
    A is a m x n constraint co-efficient matrix (usually spase)

The information we store about the problem is logically divided into - rows (m-dims), columns (n-dims) and constraint matrix (m x n dims)

Presolve_Element is the overall abstract type for storing the information about the LP problem.
The three subtypes are Presolve_Row , Presolve_Col , Presolve_Matrix
`

`
--- Active Rows and Independent Columns ---
The Presolve_Row and Presolve_Col data types defined below holds two doubly linked lists within it.
First - the normal previous and next references for rows(cols) of constraint matrix.
Second - a doubly linked list that signifies "active" rows ("independent cols"). This is explained here.

Active row means that row constraint is still to be processed.
Independent col means that corresponding variabls xj is still independent of other variables and we havent tested it yet.

Initially all rows are made "active" and all columns are made "Independent". The following is a description for rows. A similar workflow is adopted for columns.
In the function presolver! every row which is initially active is traversed and we first "deactivate" or "deque" them
After they are dequed from the active queue, they are then analyzed for possible redundancie.
If we are able to detect a redundancy we will appropriately remove the row.
If no redundancies can be detected, they have only been dequed and can still be accessed by the normal prev and next references.
At the end of the presolver! call, there are no active rows.
We make the new optimization problem by using the prev and next references among the rows and columns that havent been deleted yet.
`

abstract Presolve_Element

type Presolve_Row <: Presolve_Element
    # Each row represents a constraint row i
    i :: Int                                # Row number in the original problem.
    b_val :: Float64                        # b[i] in the original problem.
    aij :: Union{Presolve_Element,Nothing}  # Reference to the first non-zero matrix entry in row i. aij = nothing indicates empty row.
    prev :: Presolve_Row                    # Reference to previous row in the original problem. Row 5 holds ref to Row 4 etc
    next :: Presolve_Row                    # Reference to next row in the original problem. Row 5 holds ref to Row 6 etc
    is_active :: Bool                       # true if the row is in the active doubly linked list. false otherwise.
    active_prev :: Presolve_Row             # if in the active doubly linked list, reference to the previous row in the active doubly linked list.
    active_next :: Presolve_Row             # if in the active doubly linked list, reference to the next row in the active doubly linked list.
    function Presolve_Row()
        # Constructor that creates a default row which has an invalid row.i value.
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
    # Each column represents a variable xj
    j :: Int                                # Col number in the original problem.
    c_val :: Float64                        # c[j] in the original problem.
    lb :: Float64                           # lb[j] in the original problem
    ub :: Float64                           # ub[j] in the original problem
    aij :: Union{Presolve_Element,Nothing}  # Reference to the first non-zero matrix entry in col j. aij = nothing indicates empty col.
    prev :: Presolve_Col                    # Reference to previous Col in the original problem. Col 5 holds ref to Col 4 etc
    next :: Presolve_Col                    # Reference to next Col in the original problem. Col 5 holds ref to Col 6 etc
    is_independent :: Bool                  # true if the Col is in the independent doubly linked list. false otherwise.
    ind_prev :: Presolve_Col                # if in the independent doubly linked list, reference to the previous Col in the independent doubly linked list.
    ind_next :: Presolve_Col                # if in the independent doubly linked list, reference to the next Col in the independent doubly linked list.
    function Presolve_Col()
        # Constructor that creates a default row which has an invalid row.i value.
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
    row :: Presolve_Row                     # the row associated with element aij. for example, a(4,5) will hold row 4.
    col :: Presolve_Col                     # the col associated with element aij. for example, a(4,5) will hold col 5.
    val :: Float64                          # the matrix value in the constraint matrix
    row_prev :: Presolve_Matrix             # reference to the previous aij element along the same row.
    row_next :: Presolve_Matrix             # reference to the next aij element along the same row.
    col_prev :: Presolve_Matrix             # reference to the previous aij element along the same col.
    col_next :: Presolve_Matrix             # reference to the next aij element along the same col.
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

`
--- Presolve Stack ---
There are multiple redundancies possible , we combine those that are similar and can be resolved together
into a subtype of the abstract type Presolve_Stack.
As of now we have implemented Linear_Dependency which can help resolve - singleton rows, singleton columns and forcing constraints

The overall abstract type is called "Presolve_Stack" as we need to do the post-solving in the reverse order
and hence we refer to it as a stack for the LIFO logic.

Each subtype contains only as much information as is required for the postsolving.
`

abstract Presolve_Stack

`
--- Linear Dependency ---
Here we detect that a variable xj is linearly  dependent on some other variables by the equation -
xj = constant + ∑ xvalues * linear co-efficients

Consider a column singleton in column 5 and the corresponding element being a(4,5) (let dimensions be m,n = 6,6)
Suppose the fourth row of the constraint matrix looks like this -
                                b
4th row :   1 4 0 2 3 -1        10

This represents the equation -
1*x_1 + 4*x_2 + 0*x_3 + 2*x_4 + 3*x_5 + (-1)*x_6 = 10.

We can substitute x5 out of the problem as x5  = 10/3 - (1/3)*x_1 - (4/3)*x2 - (0/3)*x3 - (2/3)*x4 - (-1/3)*x6

We view this as
x[index] = value + ∑ x[vec1[i]]*vec2[i]
where,
index here is 5 for x_5
value here is 10/3
vec1 here is [1,2,4,6] (not 3 as coefficient of 3 is 0)
vec2 here is [-1/3, -4/3, -2/3, 1/3]

For a row singleton vec1 and vec2 are empty.
`

type Linear_Dependency <: Presolve_Stack
    index :: Int                            # index of the variable that is being presolved wrt original problem.
    vec1 :: Vector{Int}                     # indices of the variables that x[index] is dependent on
    vec2 :: Vector{Float64}                 # corresponding co-efficients of dependency. See above for an example
    value :: Float64                        # constant value of the Linear Dependence equation

    function Linear_Dependency(ind::Int, val::Number)
        vec1 = Array{Int,1}()
        vec2 = Array{Float64,1}()
        new(ind,vec1,vec2,val)
    end

    function Linear_Dependency(ind::Int, vec1::Vector{Int}, vec2::Vector{Float64}, val::Number)
        new(ind,vec1,vec2,val)
    end

end

`
--- Presolve Problem ---
We take the original problem and work with internally before reporting back the smaller resultant problem.
our internal workspace consist of problem type Presolve_Problem which holds information in the way we want for our internal functions.
This is never accessed by the user and its scope is the presolver! function call.

We have a constructor which initializes the variables to default values. A function make_presolve which
`

type Presolve_Problem
    # The dimensions of the original problem
    originalm :: Int64                      # number of rows in original problem
    originaln :: Int64                      # number of cols in original problem

    # Linked List storage
    dictrow :: Dict{Int64,Presolve_Row}     # Rows of the original problem.
    dictcol :: Dict{Int64,Presolve_Col}     # Cols of the original problem.
    dictaij :: Dict{Int64,Presolve_Matrix}  # Non-zero entries of the constraint matrix
    rowptr :: Presolve_Row                  # Reference to the first valid row. Will be updated as rows are deleted.
    colptr :: Presolve_Col                  # Reference to the first valid col. Will be updated as cols are deleted.
    rowque :: Presolve_Row                  # Reference to the first active row. Will be updated as rows are dequed.
    colque :: Presolve_Col                  # Reference to the first independent col. Will be updated as cols are dequed.

    # Boolean status fields
    independentvar :: BitArray{1}
    activeconstr :: BitArray{1}

    # counter variables for aij elements in row/col. Can be done away with probably.
    rowcounter :: Array{Float64,1}
    colcounter :: Array{Float64,1}

    # the stack that will be fed into the postsolver.
    pstack :: Array{Presolve_Stack,1}

    # map from the index of the inital rows to final rows. -1 if they are deleted.
    finalrows :: Array{Int,1}
    finalcols :: Array{Int,1}

    # Constructor that creates a Presolve_Problem
    function Presolve_Problem(verbose::Bool,m::Int,n::Int)
        verbose && println("-----------INSIDE PRESOLVE CONSTRUCTOR------------")

        originalm,originaln = m,n
        dictrow = Dict{Int64,Presolve_Row}()
        dictcol = Dict{Int64,Presolve_Col}()
        dictaij = Dict{Int64,Presolve_Matrix}()
        rowptr = Presolve_Row()
        colptr = Presolve_Col()
        rowque = Presolve_Row()
        colque = Presolve_Col()
        independentvar = trues(originaln)
        activeconstr = trues(originalm)
        rowcounter = zeros(originalm)
        colcounter = zeros(originaln)
        pstack = Array{Presolve_Stack,1}()
        finalrows = fill(-1,originalm)              # Initially everything is -1.
        finalcols = fill(-1,originaln)

        new(originalm,originaln,dictrow,dictcol,dictaij,rowptr,colptr,rowque,colque,independentvar,activeconstr,rowcounter,colcounter,pstack,finalrows,finalcols)
    end
end

`
--- Row and Column Operations ---
There are four functions for each - add, enque, deque and remove.
There are separate functions for enque and deque intentionally.
There are times we want to remove a row from the active list but not remove it from the problem.
Julia doesnt have NULL. Thus, self-referencing indicates end of the list in either direction.
Initially, the aij pointers are assigned "nothing". They are updated later by the add_aij functions.
`

# Creates row i with b[i] value b_val and adds it to Presolve_Problem p.
function add_row!(verbose::Bool, p::Presolve_Problem, i::Int, b_val::Float64)
    v = verbose
    v && println("Adding row : $(i)")

    row = Presolve_Row()
    row.i = i
    row.b_val = b_val
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
    enque_row!(v,p,row)
    p.dictrow[i] = row
end

# places the specified row in the active list.
function enque_row!(verbose::Bool, p::Presolve_Problem, row::Presolve_Row)
    v = verbose
    v && println("Queueing row : $(row.i)")

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

# removes the specified row from the active list.
function deque_row!(verbose::Bool, p::Presolve_Problem, row::Presolve_Row)
    v = verbose
    v && println("Dequeueing row : $(row.i)")

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

# Creates col j with c[i] value c_val and adds it to Presolve_Problem p.
function add_col!(verbose::Bool, p::Presolve_Problem, j::Int, c_val::Float64, lb::Float64, ub::Float64)
    v = verbose
    v && println("Adding col : $(j)")

    col = Presolve_Col()
    col.j = j
    col.c_val = c_val
    col.lb = lb
    col.ub = ub
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
    enque_col!(v,p,col)
    p.dictcol[j] = col
end

# places the specified col in the independent list.
function enque_col!(verbose::Bool, p::Presolve_Problem, col::Presolve_Col)
    v = verbose
    v && println("Queueing col : $(col.j)")

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

# removes the specified col from the independent list.
function deque_col!(verbose::Bool, p::Presolve_Problem, col::Presolve_Col)
    v = verbose
    v && println("Dequeueing col : $(col.j)")

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

function remove_row!(verbose::Bool, p::Presolve_Problem, row::Presolve_Row)
    v = verbose
    deque_row!(v,p,row)

    while(row.aij != nothing)
        v && println("INSIDE------ and aij is $(row.aij.row.i),$(row.aij.col.j)")
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

function remove_col!(verbose::Bool, p::Presolve_Problem, col::Presolve_Col)
    v = verbose
    v && println("REMOVE COLUMN CALL")
    deque_col!(v,p,col)

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

`
--- Matrix Element Operations ---
Matrix elements have only two functions - add_aij_normal and add_aij_transpose
They are never separately removed outside of removing rows or columns.
The input is a sparse matrix stored in the CSC format.
It can efficiently be accessed only in a column major order.
We need our matrix elements to act as two doubly linked lists.
One along the row and one along the column.
add_aij_normal creates the links along the column.
add_aij_transpose creates the links along the row.

We traverse the CSC matrix twice. First in the regular column major order and call add_aij_normal.
The second time we traverse the transpose(A) in column major order and call add_aij_transpose.
Note that the matrix element is already created by the time we call the second function.
`

function add_aij_normal!(verbose::Bool, p::Presolve_Problem, row_id::Int, col_id::Int, row_prev_id::Int, val::Float64)
    v = verbose
    v && println("Adding mat element (normal): $(row_id),$(col_id)")

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

function add_aij_transpose!(verbose::Bool, p::Presolve_Problem, row_id::Int, col_id::Int, col_prev_id::Int, val::Float64)
    v = verbose
    v && println("Adding mat element (transpose): $(row_id),$(col_id)")

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

`
--- Presolve Setup and Cleanup ---
make_presolve       : Sets up the linked list connections.
print_info          : Prints all the linked list information. Useful for debugging.
make_new            : Converts the final presolve problem data into the format of the original problem
`

# The function make_presolve! "makes" the links between the necessary row or col or matrix elements
function make_presolve!(verbose::Bool, p::Presolve_Problem,c::Array{Float64,1}, A::SparseMatrixCSC{Float64,Int64}, b::Array{Float64,1}, lb::Array{Float64,1}, ub::Array{Float64,1})
    v = verbose
    m,n = size(A)
    # checks to ensure input problem is valid.
    p.originalm != m && error("Wrong size of b wrt A")
    p.originaln != n && error("Wrong size of c wrt A")
    p.originaln != length(lb) && error("Wrong size of lb wrt A")
    p.originaln != length(ub) && error("Wrong size of ub wrt A")

    v && println("Row SETUP ----- ")
    for i in 1:p.originalm
        add_row!(v,p,i,b[i])
    end

    v && println("COL SETUP -----")
    for j in 1:p.originaln
        add_col!(v,p,j,c[j],lb[j],ub[j])
    end

    # Iterating through the non-zeros of sparse matrix A to construct the dictionary
    Arows = rowvals(A)
    v && println("MAT ELEMENT SETUP -----")
    Avals = nonzeros(A)
    for j = 1:p.originaln
        tmp = -1
        for i in nzrange(A,j)
            r = Arows[i]
            rcval = rc(r,j,p.originaln)
            #dictA[rcval] = Avals[i]
            p.rowcounter[r] += 1
            p.colcounter[j] += 1
            add_aij_normal!(v,p,r,j,tmp,Avals[i])
            tmp = r
        end
    end

    B = A'
    Arows = rowvals(B)
    Avals = nonzeros(B)
    for i = 1:p.originalm
        tmp = -1
        for c in nzrange(B,i)
            j = Arows[c]
            rcval = rc(i,j,p.originaln)
            add_aij_transpose!(v,p,i,j,tmp,Avals[c])
            tmp = j
        end
    end
end

# For Debugging.
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
    end
end

# Constructs the reduced problem
function make_new(verbose::Bool, p::Presolve_Problem)
    v = verbose

    currentn = 0
    newc = Array{Float64,1}()
    newlb = Array{Float64,1}()
    newub = Array{Float64,1}()

    col = p.colptr
    if(col.j != -1)
        v && println("Constructing newc,newlb,newub")
        while(col != nothing)
            v && @show col.j
            push!(newc,col.c_val)
            push!(newlb,col.lb)
            push!(newub,col.ub)

            currentn = currentn + 1
            p.finalcols[col.j] = currentn

            if(col.next == col)
                col = nothing
            else
                col = col.next
            end
        end
    end
    v && @show currentn

    currentm = 0
    newb = Array{Float64,1}()
    row = p.rowptr
    if(row.i != -1)
        v && println("Constructing newb")
        while(row != nothing)
            push!(newb,row.b_val)
            currentm = currentm + 1
            p.finalrows[row.i] = currentm
            if(row.next == row)
                row = nothing
            else
                row = row.next
            end
        end
    end
    v && @show currentm

    v && println(p.finalcols)
    v && println(p.finalrows)

    I = Array{Int64,1}()
    J = Array{Int64,1}()
    Val = Array{Float64,1}()
    col = p.colptr
    if(col.j != -1)
        v && println("Constructing the new A matrix")
        while(col != nothing)
            tmp = col.aij
            while(tmp != nothing)
                v && @show tmp.row.i , tmp.col.j
                v && @show tmp.val
                push!(J,p.finalcols[tmp.col.j])
                push!(I,p.finalrows[tmp.row.i])
                push!(Val,tmp.val)
                if(tmp.col_next == tmp)
                    tmp = nothing
                else
                    tmp = tmp.col_next
                end
            end
            if(col.next == col)
                col = nothing
            else
                col = col.next
            end
        end
    end
    newA = sparse(I,J,Val,currentm,currentn)
    return newc,newA,newb,newlb,newub
end

`
--- Presolver Core ---
empty_row!          : Processes the empty row. Removes it or reports an Infeasibility
presolver!          : Traverses the active list of rows and detects if redundancies are found. Call approporiate functions to handle redundancies.
singleton_row!      : Processes the singleton row. Deletes the row and makes changes to the constraint matrix appropriately
other functions will be added here in the future.
`

function presolver!(verbose::Bool,c::Array{Float64,1}, A::SparseMatrixCSC{Float64,Int64}, b::Array{Float64,1}, lb::Array{Float64,1}, ub::Array{Float64,1})
    v = verbose
    v && println("Making the Presolve Problem")
    p = Presolve_Problem(v,length(b),length(c))
    make_presolve!(v,p,c,A,b,lb,ub)

    v && println("PRESOLVE ROUTINES...............................")
    row = Presolve_Row()
    col = Presolve_Col()
    tmp = p.rowque

    while(tmp != nothing)
        row = tmp
        deque_row!(v,p,row)
        if(row.aij == nothing)
            empty_row!(v,p,row)
        else
            if(row.aij.row_next == row.aij)
                singleton_row!(v,p,row)
            else
                v && println("happy for now")
                #forcing_constraints!(p,row,v)
            end
        end
        if(tmp.next == tmp)
            tmp = nothing
        else
            tmp = tmp.next
        end
    end

    v && println("Making the reduced Problem")
    c,A,b,lb,ub = make_new(v,p)
    return c,A,b,lb,ub,p.independentvar,p.pstack
end

function empty_row!(verbose::Bool, p::Presolve_Problem, row::Presolve_Row)
    v = verbose
    v && println("EMPTY ROW FOUND AT $(row.i)")

    if(!roughly(row.b_val,0.0))
        error("Empty Row Infeasibility at row $row.i and b[i] is - $(row.b_val)")
    else
        remove_row!(v,p,row)
        p.activeconstr[row.i] = false
    end

    v && println("Exiting Empty Row")
end

function singleton_row!(verbose::Bool, p::Presolve_Problem, row::Presolve_Row)
    v = verbose
    v && println("SINGETONE ROW FOUND AT $(row.i)")

    i = row.aij.row.i
    j = row.aij.col.j
    matval = row.aij.val
    b_val = row.aij.row.b_val

    xj = b_val/matval
    add_to_stack!(Linear_Dependency(j,xj),p.independentvar,p.pstack)
    remove_row!(v,p,row)
    p.activeconstr[row.i] = false
    if(!haskey(p.dictcol,j))
        error("dictcol key error")
    end
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
    remove_col!(v,p,p.dictcol[j])
end

`
--- PostSolving Utilities ---
add_to_stack!       : function that will add the Linear_Dependency element to the stack.
post_solve!         : function that will post solve one Linear_Dependency element.
return_postsolved   : function that will take in the solution from solver for reduced problem and returns solution for original problem
`

function add_to_stack!(l::Linear_Dependency, independentvar::BitArray{1}, pstack::Array{Presolve_Stack,1})
    if (length(l.vec1) != length(l.vec2))
        error("vector1 size not equal to vector 2 size for LD element")
    end
    independentvar[l.index] = false # the variable at this index is not independent anymore
    push!(pstack,l)
end

function post_solve!(post_solvedX::Array{Float64,1}, l::Linear_Dependency)
    post_solvedX[l.index] = l.value

    for i in 1:length(l.vec1)
        post_solvedX[l.index] += l.vec2[i]*post_solvedX[l.vec1[i]]
        #println("made postsolved at $(l.index) to value $(post_solvedX[l.index])")
    end
end

function return_postsolved(x::Array{Float64,1}, independentvar::BitArray{1}, pstack :: Array{Presolve_Stack,1})
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

`
--- Miscellaneous Utilities ---
`

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
