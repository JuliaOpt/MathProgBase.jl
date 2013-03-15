This package provides one-shot functions for linear programming. To use it, you must install an LP solver first with 
``Pkg.add("Clp")``.

```julia
solution = linprog(c, A, sense, b, lb, ub[, kwargs])
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
linprog(c, A, b, sense[, kwargs]) = linprog(c, A, b, sense, 0, Inf[, kwargs])
```

Second version based on range constraints:

```julia
solution = linprog(c, A, rowlb, rowub, lb, ub[, kwargs])
```
where
- ``c`` is the objective vector, always in the sense of minimization
- ``A`` is the constraint matrix
- ``rowlb`` is the vector of row lower bounds.
- ``rowub`` is the vector of row upper bounds.
- ``lb`` is the vector of lower bounds on the variables
- ``ub`` is the vector of upper bounds on the variables


``kwargs`` contains solution parameters as keyword arguments. This isn't implemented yet, we're waiting for keyword arguments to appear in Julia.

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

If ``status`` is "Optimal", the other members have the following values

``objval`` - optimal objective value

``sol`` - primal solution vector

``attrs`` - a dictionary to contain other relevant attributes such as:

---
``redcost`` - dual multipliers for active variable bounds (zero if inactive)

``lambda`` - dual multipliers for active linear constraints (equalities are always active)

``colbasis`` - optimal simplex basis statuses for the variables (columns) if available. Possible values are ``:NonbasicAtLower``, ``:NonbasicAtUpper``, ``:Basic``,  ``:Superbasic``

``rowbasis`` - optimal simplex basis statuses for the constraints (rows) if available. Same statuses.

---

By convention the dual multipliers should have the sign following the interpretation of marginal change in the objective when the corresponding active right-hand side is increased. This corresponds to the standard that reduced costs should be nonnegative when a variable is at a lower bound and nonpositive when a variable is at an upper bound. Different solvers might have different conventions for the ``lambda`` vector, so transformations might be needed.

