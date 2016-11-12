----------------
Solver callbacks
----------------

This section describes the standardized methods for solver callbacks. This callback interface is targeted primarily at interacting with branch-and-bound solvers which are based on solving continous relaxations to find bounds, combined with heuristics to find solutions which satisfy integrality constraints. Currently there is support for adding:

*    Lazy constraints (only added to model if violated by integer-feasible solution)
*    Cut callbacks (only cuts off non-integer feasible solutions)
*    Heuristic callbacks (proposes heuristically constructed integer-feasible solutions)

A more detailed description of the three types of supported callbacks can be found in the JuMP documentation `here <http://jump.readthedocs.org/en/latest/callbacks.html>`_.

The ``SolverInterface`` module exports an abstract type ``MathProgCallbackData`` which represents the solver-specific data needed to implement the callback.

If a callback function returns ``:Exit``, the solver is expected to terminate with ``UserLimit`` status.

.. function:: setlazycallback!(m::AbstractMathProgModel,f)

   Adds lazy constraint callback ``f`` to the model. Function ``f`` takes as argument only a ``MathProgCallbackData`` object.

.. function:: setcutcallback!(m::AbstractMathProgModel,f)

   Adds cut callback ``f`` to the model. Function ``f`` takes as argument only a ``MathProgCallbackData`` object.

.. function:: setheuristiccallback!(m::AbstractMathProgModel,f)

   Adds heuristic callback ``f`` to the model. Function ``f`` takes as argument only a ``MathProgCallbackData`` object.

.. function:: setinfocallback!(m::AbstractMathProgModel,f)

   Adds informational callback ``f`` to the model. Function ``f`` takes as argument only a ``MathProgCallbackData`` object.

.. function:: cbgetmipsolution(d::MathProgCallbackData[, output])

   Grabs current best integer-feasible solution to the model. The optional second argument specifies an output vector.

.. function:: cbgetlpsolution(d::MathProgCallbackData[, output])

   Grabs current relaxation solution to the model. The optional second argument specifies an output vector.

.. function:: cbgetobj(d::MathProgCallbackData)

   Grabs objective value for current best integer-feasible solution.

.. function:: cbgetbestbound(d::MathProgCallbackData)

   Grabs best bound for objective function found so far (lower bound when minimizing, upper bound when maximizing).

.. function:: cbgetexplorednodes(d::MathProgCallbackData)

   Returns number of nodes that have been explored so far in the solve process.

.. function:: cbgetstate(d::MathProgCallbackData)

   Returns current location in solve process: ``:MIPNode`` if at node in branch-and-cut tree, ``:MIPSol`` at an integer-feasible solution, and ``:Intermediate`` otherwise. 
   
   *    ``MIPNode``: when we are at a node in the branch-and-cut tree. This is generally used for access to the solution of a relaxation (via ``cbgetlpsolution()``), or for adding cuts (via ``cbaddcut!()`` or ``cbaddcutlocal!()``).
   
   *    ``MIPSol``: when we have found a new MIP incumbent (an integer-feasible solution). This is generally used for keeping track of the intermediate solutions generated (via ``cbgetmipsolution()``) inside the branch-and-cut tree.
   
   *    ``Intermediate``: when we are still in the process of MIP or during iterations of a continuous solver. For MIPs, this is generally be used for keeping track of pessimistic and optimistic bounds (via ``cbgetobj()`` and ``cbgetbestbound()`` respectively), or the number of explored nodes (via ``cbgetexplorednodes()``) in the branch-and-cut tree.

.. function:: cbaddcut!(d::MathProgCallbackData,varidx,varcoef,sense,rhs)

   Adds cut to model. The coefficient values are represented sparsely, with (one-indexed) indices in ``varidx`` and values in ``varcoef``. The constraint sense ``sense`` is a character taking value ``<``, ``>``, or ``=``, and the right-hand side value is ``rhs``.

.. function:: cbaddcutlocal!(d::MathProgCallbackData,varidx,varcoef,sense,rhs)

   Adds local cut to model. It works as ``cbaddcut!`` but the cut is local in the sense that it only applies to the current node and the subtree rooted at this node.

.. function:: cbaddlazy!(d::MathProgCallbackData,varidx,varcoef,sense,rhs)

   Adds lazy constraint to model. The coefficient values are represented sparsely, with (one-indexed) indices in ``varidx`` and values in ``varcoef``. The constraint sense ``sense`` is a character taking value ``<``, ``>``, or ``=``, and the right-hand side value is ``rhs``.
 
.. function:: cbaddlazylocal!(d::MathProgCallbackData,varidx,varcoef,sense,rhs)

   Adds local lazy constraint to model. It works as ``cbaddlazy!`` but the lazy constraint is local in the sense that it only applies to the current node and the subtree rooted at this node.

.. function:: cbaddsolution!(d::MathProgCallbackData)

   Submit a (possibly partially defined) heuristic solution for the model. Should reset the solution stored in ``d`` to the original state at the start of callback.

.. function:: cbsetsolutionvalue!(d::MathProgCallbackData,varidx,value)

   Sets the value of a variable with (one-based) index ``varidx`` to ``value`` in the current partial solution being constructed by a user heuristic.
