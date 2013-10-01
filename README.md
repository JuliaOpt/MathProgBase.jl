This package provides one-shot functions for linear and mixed-integer programming. To use it, you must install an external solver.

Supported LP solvers are ``Clp``, ``Gurobi``, and ``GLPK``. Supported MIP solvers are ``Cbc``, ``Gurobi``, and ``GLPK``. The corresponding solver packages (of the same name, except for ``GLPKMathProgInterface``) must be installed first. Solvers are chosen by passing solver objects as optional final parameters to any call to ``linprog`` or ``mixintprog``.

---
For linear programming:

```julia
solution = linprog(c, A, sense, b, lb, ub, solver=LPSolver())
```
where
- ``c`` is the objective vector, always in the sense of minimization
- ``A`` is the constraint matrix 
- ``sense`` is the vector of constraint sense characters: ``'<'``, ``'='``, and ``'>'``. 
- ``b`` is the right-hand side vector.
- ``lb`` is the vector of lower bounds on the variables
- ``ub`` is the vector of upper bounds on the variables

A scalar is accepted for the ``b``, ``sense``, ``lb``, and ``ub`` arguments, in which case its value is replicated. ``-Inf`` and ``Inf`` are interpreted to mean that there is no corresponding lower or upper bound.


A shortened version is available as 
```julia
linprog(c, A, sense, b, solver = LPSolver()) = linprog(c, A, sense, b, 0, Inf, solver)
```

Second version based on range constraints:

```julia
solution = linprog(c, A, rowlb, rowub, lb, ub, solver=LPSolver())
```
where
- ``c`` is the objective vector, always in the sense of minimization
- ``A`` is the constraint matrix
- ``rowlb`` is the vector of row lower bounds.
- ``rowub`` is the vector of row upper bounds.
- ``lb`` is the vector of lower bounds on the variables
- ``ub`` is the vector of upper bounds on the variables


``options`` contains solution parameters as named arguments. This is only partially implemented; options may be ignored.

``solution`` is an instance of
```julia
type LinprogSolution
    status
    objval
    sol
    attrs
end
```
where
``status`` is a termination status symbol, one of ``:Optimal``, ``:Infeasible``, ``:Unbounded``, ``:UserLimit`` (iteration limit or timeout), ``:Error``, maybe others.

If ``status`` is ``:Optimal``, the other members have the following values

``objval`` - optimal objective value

``sol`` - primal solution vector

``attrs`` - a dictionary to contain other relevant attributes such as:

---
``redcost`` - dual multipliers for active variable bounds (zero if inactive)

``lambda`` - dual multipliers for active linear constraints (equalities are always active)

``colbasis`` - optimal simplex basis statuses for the variables (columns) if available. Possible values are ``:NonbasicAtLower``, ``:NonbasicAtUpper``, ``:Basic``,  ``:Superbasic``. (not yet implemented)

``rowbasis`` - optimal simplex basis statuses for the constraints (rows) if available. Same statuses. (not yet implemented)

---

By convention the dual multipliers have the sign following the interpretation of marginal change in the objective when the corresponding active right-hand side is increased. This corresponds to the standard that reduced costs should be nonnegative when a variable is at a lower bound and nonpositive when a variable is at an upper bound when minimizing. 

---

For mixed-integer programming:

```julia
solution = mixintprog(c, A, sense, b, vartypes, lb, ub, solver=MIPSolver())
```
where
- ``c`` is the objective vector, always in the sense of minimization
- ``A`` is the constraint matrix 
- ``sense`` is the vector of constraint sense characters: ``'<'``, ``'='``, and ``'>'``. 
- ``b`` is the right-hand side vector.
- ``vartypes`` is a ``Vector{Char}`` with ``'I'`` indicating integer variables and ``'C'`` indicating continous variables. Binary variables should be coded as ``'I'`` with lower bounds of 0 and upper bounds of 1.
- ``lb`` is the vector of lower bounds on the variables
- ``ub`` is the vector of upper bounds on the variables

``solution`` is an instance of
```julia
type MixintprogSolution
    status
    objval
    sol
    attrs
end
```
where
``status`` is as above. 

If ``status`` is ``:Optimal``, the other members have the following values

``objval`` - optimal objective value

``sol`` - primal solution vector

``attrs`` - a dictionary to contain other relevant attributes such as:

``objbound`` - Best known lower bound on objective value.


Shortened and range-constraint versions are available as well.

---

The ``LPSolver`` and ``MIPSolver`` constructors accept the name of the solver as a symbol and an optional list of solver parameters. For example, ``MIPSolver(:Cbc,LogLevel=1)`` will create a solver object that when passed to ``mixintprog`` will solve the problem with Cbc with verbose output. Solver parameters are entirely solver dependent and have not yet been well documented.
