-------------
MIP Callbacks
-------------
MathProgBase supports a standardized and abstracted way to implement common MIP callbacks on the model. Currently there is support for adding:

*    Lazy constraints (only added to model if violated by integer-feasible solution)
*    Cut callbacks (only cuts off non-integer feasible solutions)
*    Heuristic callbacks (proposes heuristically constructed integer-feasible solutions at MIP nodes)

A more detailed description of the three types of supported callbacks can be found in the JuMP documentation `here <https://jump.readthedocs.org/en/latest/jump.html#solver-callbacks>`_.

The ``MathProgSolverInterface`` exports an abstract type ``MathProgCallbackData`` which represents the solver-specific data needed to implement the callback.

.. function:: setlazycallback!(m::AbstractMathProgModel,f)

   Adds lazy constraint callback ``f`` to the model. Function ``f`` takes as argument only a ``MathProgCallbackData`` object.
   
.. function:: setcutcallback!(m::AbstractMathProgModel,f)

   Adds cut callback ``f`` to the model. Function ``f`` takes as argument only a ``MathProgCallbackData`` object.
   
.. function:: setheuristiccallback!(m::AbstractMathProgModel,f)

   Adds heuristic callback ``f`` to the model. Function ``f`` takes as argument only a ``MathProgCallbackData`` object.
   
.. function:: cbgetmipsolution(d::MathProgCallbackData[, output])

   Grabs current best integer-feasible solution to the model. The optional second argument specifies an output vector.
   
.. function:: cbgetlpsolution(d::MathProgCallbackData[, output])

   Grabs current best linear relaxation solution to the model. The optional second argument specifies an output vector.
   
.. function:: cbgetobj(d::MathProgCallbackData)

   Grabs objective value for current best integer-feasible solution.
   
.. function:: cbgetbestbound(d::MathProgCallbackData) 

   Grabs best bound for objective function found so far (lower bound when minimizing, upper bound when maximizing).
   
.. function:: cbgetexplorednodes(d::MathProgCallbackData)

   Returns number of nodes that have been explored so far in the solve process.

.. function:: cbgetstate(d::MathProgCallbackData)

   Returns current location in solve process: ``:MIPNode`` if at node in branch-and-cut tree, ``:MIPSol`` at an integer-feasible solution, and ``:Other`` otherwise.

.. function:: cbaddcut!(d::MathProgCallbackData,varidx,varcoef,sense,rhs) 

   Adds cut to model. The coefficient values are represented sparsely, with (one-indexed) indices in ``varidx`` and values in ``varcoef``. The constraint sense ``sense`` is a character taking value ``<``, ``>``, or ``=``, and the right-hand side value is ``rhs``.
   
.. function:: cbaddlazy!(d::MathProgCallbackData,varidx,varcoef,sense,rhs)

   Adds lazy constraint to model. The coefficient values are represented sparsely, with (one-indexed) indices in ``varidx`` and values in ``varcoef``. The constraint sense ``sense`` is a character taking value ``<``, ``>``, or ``=``, and the right-hand side value is ``rhs``.

.. function:: cbaddsolution!(d::MathProgCallbackData)

   Submit a (possibly partially defined) heuristic solution for the model. Should reset the solution stored in ``d`` to the original state at the start of callback.

.. function:: cbsetsolutionvalue!(d::MathProgCallbackData,varidx,value)

   Sets the value of a variable with (one-based) index ``varidx`` to ``value`` in the current partial solution being constructed by a user heuristic.
