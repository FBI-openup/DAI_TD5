import sys
from itertools import chain
from fractions import Fraction

from pysmt.shortcuts import *
from pysmt.logics import QF_LRA
from pysmt.typing import REAL
from pyvmt.environment import reset_env, get_env

# from pysmt.solvers.msat import MSatEnv, MSatConverter  # not needed
from pysmt.solvers.z3 import Z3Converter, Z3Solver

from simplex import Context, Simplex

def main():
    context = Context()
    context.env.enable_infix_notation = True

    x,y,z,v,w = [Symbol(a,REAL) for a in ['x','y','z','v','w']]

    systems = [
        [ 2 * y - x - 2 <= 0,
          -2 * y - x + 4 <= 0],
        [ x - y >= -1,
          y <= 4,
          x + y >= 6,
          3 * x - y <= 7],
        [ x - y >= 1,
          y >= 1,
          y + x <= 0],
        [ x >= 1,
          y >= 1,
          y + x <= 0],
        [1*x + 1*y <= 4,
         1*x <= 3,
         1*y <= 3,
         -1*x <= 0,
         -1*y <= 0],
        [1*x + 1*y <= 2,
         1*x + 1*y >= 5,
         -1*x <= 0,
         -1*y <= 0],
        [2*x + 2*y <= 8,
         1*x + 1*y <= 4,
         1*x <= 3,
         1*y <= 3,
         -1*x <= 0,
         -1*y <= 0],
        [1*x + 1*y <= 4,
         1*x <= 2,
         1*y <= 2,
         -1*x <= 0,
         -1*y <= 0],
        [1*x - 1*y <= 2,
         -1*y <= 0,
         -1*x <= 0],
        [1*x + 1*y + 1*z <= 6,
         1*x <= 4,
         1*y <= 4,
         1*z <= 4,
         -1*x <= 0,
         -1*y <= 0,
         -1*z <= 0],
        [1*x + 1*y <= 3,
         1*x >= 4,
         -1*y <= 0],
        [1*x + 1*y <= 5,
         -1*x - 1*y <= -5,
         -1*x <= 0,
         -1*y <= 0],
        [1*x - 11*y - 5*z + 18*v <= 0,
         1*x - 3*y - 1*z + 2*v <= 0,
         1*x <= 1,
         -1*x <= 0,
         -1*y <= 0,
         -1*z <= 0,
         -1*v <= 0],
    ]

    for system in systems:
        s = Simplex()
        s.preprocess(context, system)
        s.print()

        res = s.solve()
        if res:
            print("feasible")
        else:
            print("unfeasible")

        ground = is_sat(And(system))
        assert res == ground


if __name__ == "__main__":
    main()
