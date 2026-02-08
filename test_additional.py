#!/usr/bin/env python3
# ZHANG Boyuan DAI TD5
# Additional test cases for simplex solver

import sys
from fractions import Fraction
from pysmt.shortcuts import *
from pysmt.typing import REAL
from simplex import Context, Simplex

def test_case(context, system, expected):
    """Run a single test case with the same style as driver.py"""
    s = Simplex()
    s.preprocess(context, system)
    s.print()
    result = s.solve()
    print('feasible' if result else 'unfeasible')
    assert result == expected, f"Expected {'feasible' if expected else 'unfeasible'}, got {'feasible' if result else 'unfeasible'}"

def main():
    context = Context()
    context.env.enable_infix_notation = True
    
    x, y, z = [Symbol(a, REAL) for a in ['x', 'y', 'z']]
    
    # Test 1: Single constraint
    test_case(context, [x <= 10], True)
    
    # Test 2: Equality constraint (x + y = 5)
    test_case(context, [
        x + y <= 5,
        x + y >= 5,
        x >= 0,
        y >= 0
    ], True)
    
    # Test 3: Unbounded feasible region
    test_case(context, [
        x >= 0,
        y >= 0,
        x - y <= 3
    ], True)
    
    # Test 4: Tight constraints (unique point)
    test_case(context, [
        x + y >= 10,
        x + y <= 10,
        x >= 5,
        y >= 5
    ], True)
    
    # Test 5: Contradictory constraints
    test_case(context, [
        x <= 5,
        x >= 10
    ], False)
    
    # Test 6: Three variables
    test_case(context, [
        x + y + z <= 10,
        x >= 1,
        y >= 2,
        z >= 3,
        2*x + y <= 8
    ], True)

if __name__ == "__main__":
    main()
