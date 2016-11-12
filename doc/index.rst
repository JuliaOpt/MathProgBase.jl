===============
MathProgBase.jl
===============

.. module:: MathProgBase
   :synopsis: Solver-independent functions and low-level interfaces for Mathematical Programming

`MathProgBase.jl <https://github.com/JuliaOpt/MathProgBase.jl>`_ provides high-level
one-shot functions for linear and mixed-integer programming,
as well as a solver-independent low-level interface for implementing advanced techniques
that require efficiently solving a sequence of linear programming problems.

To use MathProgBase, an external solver must be installed. See :ref:`choosing solvers <choosing-solvers>`.


Contents
--------

.. toctree::
    :maxdepth: 2
    
    linprog.rst
    mixintprog.rst
    quadprog.rst
    solverinterface.rst
    solvers.rst
    lpqcqp.rst
    conic.rst
    nlp.rst
    callbacks.rst


