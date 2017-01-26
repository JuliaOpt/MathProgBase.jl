-------------------------
Mixed-integer Programming
-------------------------

.. function:: mixintprog(c, A, sense, b, vartypes, lb, ub, solver)

Solves the same optimization problem as ``linprog`` above, except variables
are additionally constrained to take only integer values if the corresponding
entry in the ``varypes`` vector is the symbol ``:Int``. Continuous
variables are indicated by the value ``:Cont``, binary variables should be specified by ``:Bin``, 
semicontinuous by ``:SemiCont``, and semi-integer by ``:SemiInt``.

A scalar is accepted for the ``sense``, ``b``, ``vartypes``, ``lb``, and ``ub`` arguments, in which case its value is replicated. The values ``-Inf`` and ``Inf`` are interpreted to mean that there is no corresponding lower or upper bound. 

The ``mixintprog`` function returns an instance of the type::
    
    type MixintprogSolution
        status
        objval
        sol
        attrs
    end

where ``status`` takes the same values as with ``linprog``.

If ``status`` does not indicate error or infeasiblity, the other members have the following values:

* ``objval`` -- optimal objective value
* ``sol`` -- primal solution vector
* ``attrs`` -- a dictionary that may contain other relevant attributes such as:

  - ``objbound`` -- Best known lower bound on the objective value

Analogous shortened and range-constraint versions are available as well.

We can solve a `binary knapsack problem <http://en.wikipedia.org/wiki/Knapsack_problem>`_

.. math::
    max\, &5x_1 + 3x_2 + 2x_3 + 7x_4 + 4x_5\\
    s.t.  &2x_1 + 8x_2 + 4x_3 + 2x_4 + 5x_5 \leq 10\\
          & (x_1, x_2, x_3, x_4, x_5) \in \{0,1\}^5

with the code::
    
    mixintprog(-[5.,3.,2.,7.,4.],[2. 8. 4. 2. 5.],'<',10,:Int,0,1,CbcSolver())
