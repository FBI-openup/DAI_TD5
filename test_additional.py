#!/usr/bin/env python3
# ZHANG Boyuan DAI TD5
# Additional test cases for simplex solver

import sys
from fractions import Fraction
from pysmt.shortcuts import *
from pysmt.typing import REAL
from simplex import Context, Simplex

def test_single_constraint():
    """Test with a single simple constraint"""
    print("\n=== Test 1: Single Constraint ===")
    context = Context()
    context.env.enable_infix_notation = True
    
    x = Symbol('x', REAL)
    system = [x <= 10]
    
    s = Simplex()
    s.preprocess(context, system)
    result = s.solve()
    
    print(f"System: x <= 10")
    print(f"Result: {'feasible' if result else 'unfeasible'}")
    print(f"Expected: feasible")
    assert result == True
    print("✓ PASSED\n")

def test_equality_constraints():
    """Test with equality constraints (converted to two inequalities)"""
    print("=== Test 2: Equality Constraints ===")
    context = Context()
    context.env.enable_infix_notation = True
    
    x, y = [Symbol(a, REAL) for a in ['x', 'y']]
    # x + y = 5 is represented as x + y <= 5 AND x + y >= 5
    system = [
        x + y <= 5,
        x + y >= 5,
        x >= 0,
        y >= 0
    ]
    
    s = Simplex()
    s.preprocess(context, system)
    result = s.solve()
    
    print(f"System: x + y = 5, x >= 0, y >= 0")
    print(f"Result: {'feasible' if result else 'unfeasible'}")
    print(f"Expected: feasible")
    assert result == True
    print("✓ PASSED\n")

def test_unbounded_region():
    """Test with unbounded feasible region"""
    print("=== Test 3: Unbounded Feasible Region ===")
    context = Context()
    context.env.enable_infix_notation = True
    
    x, y = [Symbol(a, REAL) for a in ['x', 'y']]
    system = [
        x >= 0,
        y >= 0,
        x - y <= 3
    ]
    
    s = Simplex()
    s.preprocess(context, system)
    result = s.solve()
    
    print(f"System: x >= 0, y >= 0, x - y <= 3")
    print(f"Result: {'feasible' if result else 'unfeasible'}")
    print(f"Expected: feasible (unbounded)")
    assert result == True
    print("✓ PASSED\n")

def test_tight_constraints():
    """Test with very tight constraints"""
    print("=== Test 4: Tight Constraints ===")
    context = Context()
    context.env.enable_infix_notation = True
    
    x, y = [Symbol(a, REAL) for a in ['x', 'y']]
    system = [
        x + y >= 10,
        x + y <= 10,  # Forces x + y = 10
        x >= 5,
        y >= 5
    ]
    
    s = Simplex()
    s.preprocess(context, system)
    result = s.solve()
    
    print(f"System: x + y = 10, x >= 5, y >= 5")
    print(f"Result: {'feasible' if result else 'unfeasible'}")
    print(f"Expected: feasible (unique point: x=5, y=5)")
    assert result == True
    print("✓ PASSED\n")

def test_contradictory():
    """Test obviously contradictory constraints"""
    print("=== Test 5: Contradictory Constraints ===")
    context = Context()
    context.env.enable_infix_notation = True
    
    x = Symbol('x', REAL)
    system = [
        x <= 5,
        x >= 10
    ]
    
    s = Simplex()
    s.preprocess(context, system)
    result = s.solve()
    
    print(f"System: x <= 5 AND x >= 10")
    print(f"Result: {'feasible' if result else 'unfeasible'}")
    print(f"Expected: unfeasible (contradiction)")
    assert result == False
    print("✓ PASSED\n")

def test_three_variables():
    """Test with 3 variables"""
    print("=== Test 6: Three Variables ===")
    context = Context()
    context.env.enable_infix_notation = True
    
    x, y, z = [Symbol(a, REAL) for a in ['x', 'y', 'z']]
    system = [
        x + y + z <= 10,
        x >= 1,
        y >= 2,
        z >= 3,
        2*x + y <= 8
    ]
    
    s = Simplex()
    s.preprocess(context, system)
    result = s.solve()
    
    print(f"System: x + y + z <= 10, x >= 1, y >= 2, z >= 3, 2x + y <= 8")
    print(f"Result: {'feasible' if result else 'unfeasible'}")
    print(f"Expected: feasible")
    assert result == True
    print("✓ PASSED\n")

def main():
    print("="*60)
    print("     Additional Simplex Solver Test Cases")
    print("="*60)
    
    try:
        test_single_constraint()
        test_equality_constraints()
        test_unbounded_region()
        test_tight_constraints()
        test_contradictory()
        test_three_variables()
        
        print("="*60)
        print("✓✓✓ All additional tests PASSED! ✓✓✓")
        print("="*60)
    except AssertionError as e:
        print(f"\n✗✗✗ Test FAILED! ✗✗✗")
        sys.exit(1)
    except Exception as e:
        print(f"\n✗✗✗ Error: {e} ✗✗✗")
        sys.exit(1)

if __name__ == "__main__":
    main()
