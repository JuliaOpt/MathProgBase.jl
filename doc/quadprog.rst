---------------------
Quadratic Programming
---------------------

.. function:: quadprog(c, Q, A, sense, b, l, u, solver)

Solves the quadratic programming problem:

.. math::
    \min_{x}\, &\frac{1}{2}x^TQx + c^Tx\\
    s.t.     &a_i^Tx \text{ sense}_i \, b_i \forall\,\, i\\
             &l \leq x \leq u\\

where:

*    ``c`` is the objective vector, always in the sense of minimization
*    ``Q`` is the Hessian matrix of the objective
*    ``A`` is the constraint matrix, with rows :math:`a_i` (viewed as column-oriented vectors)
*    ``sense`` is a vector of constraint sense characters ``'<'``, ``'='``, and ``'>'``
*    ``b`` is the right-hand side vector
*    ``l`` is the vector of lower bounds on the variables
*    ``u`` is the vector of upper bounds on the variables, and
*    ``solver`` specifies the desired solver, see :ref:`choosing solvers <choosing-solvers>`.
 
A scalar is accepted for the ``b``, ``sense``, ``l``, and ``u`` arguments, in which case its value is replicated. The values ``-Inf`` and ``Inf`` are interpreted to mean that there is no corresponding lower or upper bound.

.. note::
    Quadratic programming solvers extensively exploit the sparsity of the Hessian matrix ``Q`` and the constraint matrix ``A``. While both dense and sparse matrices are accepted, for large-scale problems sparse matrices should be provided if permitted by the problem structure.  

The ``quadprog`` function returns an instance of the type::
    
    type QuadprogSolution
        status
        objval
        sol
        attrs
    end

where ``status`` is a termination status symbol, one of ``:Optimal``, ``:Infeasible``, ``:Unbounded``, ``:UserLimit`` (iteration limit or timeout), ``:Error`` (and maybe others).

If ``status`` is ``:Optimal``, the other members have the following values:

* ``objval`` -- optimal objective value
* ``sol`` -- primal solution vector
* ``attrs`` -- a dictionary that may contain other relevant attributes (not currently used).

Analogous shortened and range-constraint versions are available as well.

We can solve the three-dimensional QP (see ``test/quadprog.jl``):

.. math::
    \min_{x,y,z}\, &x^2+y^2+z^2+xy+yz\\
    s.t.         &x + 2y + 3z \geq 4\\
                 &x + y \geq 1

by::

    using MathProgBase, Ipopt
    
    sol = quadprog([0., 0., 0.],[2. 1. 0.; 1. 2. 1.; 0. 1. 2.],[1. 2. 3.; 1. 1. 0.],'>',[4., 1.],-Inf,Inf, IpoptSolver())
    if sol.status == :Optimal
        println("Optimal objective value is $(sol.objval)")
        println("Optimal solution vector is: [$(sol.sol[1]), $(sol.sol[2]), $(sol.sol[3])]")
    else
        println("Error: solution status $(sol.status)")
    end
