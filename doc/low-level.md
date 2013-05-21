Working documentation for standardized low-level LP/MIP/QP/MIQP+ interface:

Conventions:
- All indices are one-based.
- It should be easy to access the solver-specific model object in order to call solver-specific functions.

---

```
model(;kwargs...)
```

Create a new empty model, optionally passing in parameters for the solver.

Notes:
- Solver-specific parameters are allowed. Eventually we will also have a standard set of parameters.
- In order to load a problem from a file in Gurobi, the file name needs to be passed to the GRBModel constructor. I suggest that this function create a GRBEnv object and wait until later to create a GRBModel.

---

```
loadproblem(m, filename::String)
```

Load a problem from a file into the model ``m``. File extension should be used to determine the format.

Notes:
- Different solvers support different input formats. An error should be thrown if an unsupported format is used.

```
loadproblem(m, A, collb, colub, obj, rowlb, rowub)
```

Load a problem from the given input data. ``A`` is the constraint matrix, ``collb`` are lower bounds on variables, ``colub`` are upper bounds on variables, ``rowlb`` are lower bounds on constraints, ``rowub`` are upper bounds on constraints.

Notes:
- The solver may choose what type of matrix to accept for ``A``. As a minimum, this should include ``Array{Float64,2}`` and ``SparseMatrixCSC{Float64,Int}``. Exact solvers may want to allow matrices with rational entries.
- ``Inf`` and ``-Inf`` indicate no corresponding upper or lower bound, respectively. Equal lower and upper bounds are used to indicate equalities.

---

```
writeproblem(m, filename::String)
```

Write problem to a file. File extension should be used to determine the format.

Notes:
- Different solvers support different output formats. An error should be thrown if an unsupported format is used.

---

```
getvarLB(m)
setvarLB(m, collb)
getvarUB(m)
setvarUB(m, colub)
getconstrLB(m)
setconstrLB(m, rowlb)
getconstrUB(m)
setconstrUB(m, rowub)
getobj(m)
setobj(m, obj)
```

Setters and getters for vector properties.

---

```
addvar(m, rowidx, rowcoef, collb, colub, objcoef)
```

Add a column to the problem.  ``rowcoef`` specifies non-zero coefficients in the new column in corresponding indices specified in ``rowcoef``. ``collb``, ``colub``, and ``objcoef`` are scalars specifying the lower bound, upper bound, and objective coefficient for the new variable.

---

```
addconstr(m, colidx, colcoef, rowlb, rowub)
```

Add a row to the problem.

---

```
updatemodel(m)
```

Update the model after making changes. Only some solvers need this.

---

```
setsense(m,sense)
getsense(m)
```
Set/get problem sense. Either ``:Min`` or ``:Max``.

---

```
numvar(m)
numconstr(m)
```

Get number of variables and constraints in the model.

---

```
optimize(m)
```

Solve (or resolve) the current model.

---

```
status(m)
```

Get termination status symbol, one of ``:Optimal``, ``:PrimalInfeasible``, ``:DualInfeasible``, ``:UserLimit`` (iteration limit or timeout), ``:Error``, maybe others.

---

```
getobjval(m)
```

Objective value of optimal (or best feasible) solution. 

```
getobjbound(m)
```

Best dual bound (lower bound, when minimizing) on optimal objective value.

```
getsolution(m)
```

Primal solution.

```
getconstrsolution(m)
```

Constraint matrix times solution vector, also known as "row activity".

```
getreducedcosts(m)
```

Dual solution corresponding to (active) variable bounds.

```
getconstrduals(m)
```

Dual solution corresponding to (active) row constraints.

Notes:
- By convention the dual values should have the sign following the interpretation of marginal change in the objective when the corresponding active right-hand side is increased. This corresponds to the standard that reduced costs should be nonnegative when a variable is at a lower bound and nonpositive when a variable is at an upper bound. Different solvers might have different conventions for the row duals, so transformations might be needed.

---

```
getrawsolver(m)
```

Return a solver-specific object that can be used directly with the solver's low-level interface.

---

```
setvartype(m,v)
getvartype(m)
```

Get/set vector of variable types. ``v`` should be a ``Char`` vector where `'C'` indicates continuous and `'I'` indicates integer. (Binaries should be specified as `'I'` with lower bound 0 and upper bound 1). 

---


What's missing:
- Setting/getting solution method
- Add multiple constraints/variables. Remove constraints/variables.
- Setting/getting simplex basis (if applicable)
- Standard subset of parameters (Timeout, output level, iteration limit, MIP gap, etc.)
- SOS constraints for integer problems
- Quadratic objectives

