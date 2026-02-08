import sys
from itertools import chain
from fractions import Fraction

from pysmt.shortcuts import *
from pysmt.logics import QF_LRA
from pysmt.typing import REAL
from pyvmt.environment import reset_env, get_env

# from pysmt.solvers.msat import MSatEnv, MSatConverter  # not needed, msat solver not installed
from pysmt.solvers.z3 import Z3Converter, Z3Solver


class Context:
    def __init__(self):
        self.env = get_env()
        self.z3 = Z3Solver(self.env, logic=QF_LRA)
        self.convert = Z3Converter(self.env, self.z3.z3.ctx)

    def norm(self, f):
        normalized = self.convert.back(self.convert.convert(f))
        return normalized

class Simplex:

    """
    Represents the tableau
    """
    def __init__(self):
        self.system = []
        self.basic = []
        self.non_basic = []
        self.tableaux = []

        # Maps from variable to lower bounds, upper bounds, and assignments
        self.lb = {}
        self.ub = {}
        self.assignments = {}

    def build_tableaux_row(self, nb_index, eq):
        # fill the tableau
        def proc_lc(lc, change_sign, line):
            if lc.is_plus():
                proc_lc(lc.args()[0],change_sign, line)
                proc_lc(lc.args()[1],change_sign, line)
                return
            if lc.is_minus():
                proc_lc(lc.args()[0],change_sign, line)
                proc_lc(lc.args()[1],not change_sign, line)
                return
            if lc.is_real_constant():
                coeff = lc.constant_value()
                symb = None
            elif lc.is_symbol():
                coeff = Fraction(1)
                symb = lc
            elif lc.is_times():
                i_c = 0 if lc.args()[0].is_real_constant() else 1
                i_s = 1 if lc.args()[1].is_symbol() else 0
                assert i_c != i_s
                coeff = lc.args()[i_c].constant_value()
                symb = lc.args()[i_s]
            else:
                raise Exception("Uknown type for " + str(lc))

            assert not coeff is None

            if not symb is None:
                if change_sign:
                    coeff = Fraction(-1) * coeff
                line[nb_index[symb]] = coeff
            else:
                if not change_sign:
                    coeff = Fraction(-1) * coeff
                self.ub[bv] = coeff

        # add basic variable
        bv = Symbol("_s_%d" % len(self.basic), REAL)
        self.basic.append(bv)

        self.lb[bv] = None
        self.ub[bv] = None

        line = [Fraction(0) for v in self.non_basic]
        proc_lc(eq.args()[0], False, line)
        proc_lc(eq.args()[1], True, line)
        self.tableaux.append(line)


    def preprocess(self, context, system):
        """
        Parse a linear system and builds the tableau
        """
        self.system = [context.norm(s) for s in system]

        fv = set()
        for f in self.system:
            fv |= f.get_free_variables()
        self.non_basic = list(fv)
        self.non_basic.sort()

        nb_index = {}
        for i in range(len(self.non_basic)):
            nb_index[self.non_basic[i]] = i
            self.lb[self.non_basic[i]] = None
            self.ub[self.non_basic[i]] = None

        for s in self.system:
            self.build_tableaux_row(nb_index, s)

        for v in chain(self.basic, self.non_basic):
            self.assignments[v] = Fraction(0)

    def lt(self, var, is_upper=True):
        # var < upper or lower bound
        # Note: use "None" as value for infinite
        val = self.assignments[var]
        bound = self.ub[var] if is_upper else self.lb[var]

        if is_upper and bound is None:
            return True # val < +oo
        elif not is_upper and bound is None:
            return False # -oo < val
        else:
            return val < bound

    def gt(self, var, is_upper=True):
        # var > upper or lower bound
        # Note: use "None" as value for infinite
        val = self.assignments[var]
        bound = self.ub[var] if is_upper else self.lb[var]

        if is_upper and bound is None:
            return False # val < +oo
        elif not is_upper and bound is None:
            return True # -oo < val
        else:
            return val > bound

    def in_bound(self, var):
        return not (self.gt(var, True) or
                    self.lt(var, False))

    def check_invariant(self):
        def check_eq():
            # assignments must satisfy equations at all times
            for i in range(len(self.basic)):
                basic_val = self.assignments[self.basic[i]]

                res = Fraction(0)
                for j in range(len(self.non_basic)):
                    nb_val = self.assignments[self.non_basic[j]]
                    nb_coeff = self.tableaux[i][j]
                    res = res + nb_coeff * nb_val

                if res != basic_val:
                    print("Check eq failed", self.basic[i],"=",basic_val,"but sum is", res)
                    return False
            return True

        def check_bounds():
            # all bounds on non-basic variables must hold
            for j in range(len(self.non_basic)):
                nb_var = self.non_basic[j]
                nb_val = self.assignments[nb_var]

                if not self.in_bound(nb_var):
                    print("Bound does not hold for",nb_var,
                          self.lb[nb_var],nb_val,self.ub[nb_var])
                    return False
            return True

        return check_eq() and check_bounds()

    def print_tableau(self):
        print("Basic", ",".join([str(v) for v in self.basic]))
        print("Non basic", ",".join([str(v) for v in self.non_basic]))

        print("Tableaux:")
        for t in self.tableaux:
            print(t)

        print("Bounds:")
        for v in chain(self.basic, self.non_basic):
            def get_bound_str(b, is_ub):
                if b is None:
                    return "+oo" if is_ub else "-oo"
                else:
                    return str(b)
            print(get_bound_str(self.lb[v], False),"<=",v,"<=",get_bound_str(self.ub[v],True))

    def print(self):
        print("System")
        for s in self.system:
            print(s.serialize())
        self.print_tableau()

    def solve(self):
        # main loop: keep going until all variables are within bounds or we find it's impossible
        max_iterations = 10000  # safety limit to prevent infinite loops during debugging
        iteration = 0
        
        # debug flag - set to True to see detailed output
        debug = False  # change to True to enable debugging
        
        while True:
            assert self.check_invariant()
            iteration += 1
            if iteration > max_iterations:
                # safety check - shouldn't happen if algorithm is correct
                print(f"Warning: reached max iterations {max_iterations}")
                return False
            
            # Step 1: find a basic variable that violates its bounds
            # (basic vars can temporarily be out of bounds, we need to fix them)
            violating_basic_idx = None
            for i in range(len(self.basic)):
                bvar = self.basic[i]
                if not self.in_bound(bvar):
                    violating_basic_idx = i
                    break
            
            # if no violation found, we're done! system is feasible
            if violating_basic_idx is None:
                return True
            
            # ok, we found a basic variable that's out of bounds, let's fix it
            bvar = self.basic[violating_basic_idx]
            bval = self.assignments[bvar]
            
            # figure out if it's too high or too low
            is_upper_violation = self.gt(bvar, is_upper=True)  # bvar > upper bound
            
            # Step 2: pick a non-basic variable to adjust
            # we need to find one that can help bring the basic var back in bounds
            selected_nb_idx = None
            
            for j in range(len(self.non_basic)):
                nbvar = self.non_basic[j]
                coeff = self.tableaux[violating_basic_idx][j]
                
                # skip if coefficient is zero (this nbvar doesn't affect bvar)
                if coeff == 0:
                    continue
                
                # the coefficient tells us how changing nbvar affects bvar
                # if coeff > 0: increasing nbvar increases bvar
                # if coeff < 0: increasing nbvar decreases bvar
                
                if is_upper_violation:
                    # bvar is too high, we need to decrease it
                    if coeff > 0:
                        # increase nbvar -> bvar goes up (bad!)
                        # so we need to decrease nbvar
                        # can we decrease? check if nbvar > lower bound
                        if self.gt(nbvar, is_upper=False):
                            selected_nb_idx = j
                            break
                    else:  # coeff < 0
                        # increase nbvar -> bvar goes down (good!)
                        # can we increase? check if nbvar < upper bound
                        if self.lt(nbvar, is_upper=True):
                            selected_nb_idx = j
                            break
                else:
                    # bvar is too low, we need to increase it
                    if coeff > 0:
                        # increase nbvar -> bvar goes up (good!)
                        # can we increase? check if nbvar < upper bound
                        if self.lt(nbvar, is_upper=True):
                            selected_nb_idx = j
                            break
                    else:  # coeff < 0
                        # increase nbvar -> bvar goes down (bad!)
                        # so we need to decrease nbvar
                        # can we decrease? check if nbvar > lower bound
                        if self.gt(nbvar, is_upper=False):
                            selected_nb_idx = j
                            break
            
            # if we can't find any non-basic variable to help, system is infeasible
            if selected_nb_idx is None:
                if debug:
                    print(f"  Iteration {iteration}: No suitable non-basic variable found -> INFEASIBLE")
                return False
            
            if debug:
                nbvar_name = self.non_basic[selected_nb_idx]
                nbvar_val = self.assignments[nbvar_name]
                coeff_val = self.tableaux[violating_basic_idx][selected_nb_idx]
                print(f"  Iteration {iteration}: Selected {nbvar_name} (current={nbvar_val}, coeff={coeff_val})")
            
            # Step 3: calculate how much to change the non-basic variable
            nbvar = self.non_basic[selected_nb_idx]
            coeff = self.tableaux[violating_basic_idx][selected_nb_idx]
            
            # we want to move bvar to its bound
            if is_upper_violation:
                target_bval = self.ub[bvar]
            else:
                target_bval = self.lb[bvar]
            
            # calculate the change needed for nbvar to bring bvar to its bound
            # bvar + coeff * delta_nb = target_bval
            # delta_nb = (target_bval - bvar) / coeff
            delta_nb = (target_bval - bval) / coeff
            new_nbval = self.assignments[nbvar] + delta_nb
            
            # Step 4: check if nbvar would hit its own bound
            # if yes, we need to pivot; if no, just update
            needs_pivot = False
            
            if delta_nb > 0:
                # we're increasing nbvar
                if self.ub[nbvar] is not None and new_nbval > self.ub[nbvar]:
                    # we hit the upper bound before reaching target
                    new_nbval = self.ub[nbvar]
                    needs_pivot = True
            else:
                # we're decreasing nbvar
                if self.lb[nbvar] is not None and new_nbval < self.lb[nbvar]:
                    # we hit the lower bound before reaching target
                    new_nbval = self.lb[nbvar]
                    needs_pivot = True
            
            # update nbvar assignment
            old_nbval = self.assignments[nbvar]
            self.assignments[nbvar] = new_nbval
            actual_delta = new_nbval - old_nbval
            
            # update all basic variables affected by this change
            # (they all depend on this non-basic variable)
            for i in range(len(self.basic)):
                bv = self.basic[i]
                c = self.tableaux[i][selected_nb_idx]
                self.assignments[bv] = self.assignments[bv] + c * actual_delta
            
            # Step 5: if we need to pivot, swap the variables
            if needs_pivot:
                # swap basic[violating_basic_idx] with non_basic[selected_nb_idx]
                # this is the trickiest part - we need to update the tableau
                
                # save the variables we're swapping
                old_basic = self.basic[violating_basic_idx]
                old_nb = self.non_basic[selected_nb_idx]
                
                # the pivot coefficient
                pivot_coeff = self.tableaux[violating_basic_idx][selected_nb_idx]
                
                # rewrite the pivot row: solve for the new basic variable (old_nb)
                # old row: old_basic = sum(a[j] * nb[j])
                # pivot element is a[selected_nb_idx] = pivot_coeff
                # new row: old_nb = (1/pivot_coeff) * old_basic - sum((a[j]/pivot_coeff) * nb[j] for j != selected)
                
                pivot_row = self.tableaux[violating_basic_idx]
                new_pivot_row = []
                for j in range(len(self.non_basic)):
                    if j == selected_nb_idx:
                        # this position will hold the old basic variable
                        # coefficient is 1/pivot_coeff
                        new_pivot_row.append(Fraction(1) / pivot_coeff)
                    else:
                        # coefficient is -a[j]/pivot_coeff
                        new_pivot_row.append(-pivot_row[j] / pivot_coeff)
                
                # update all other rows: substitute the new expression for old_nb
                # old row: basic[i] = sum(a[i][j] * nb[j])
                # we're replacing nb[selected_nb_idx] with the new expression
                
                for i in range(len(self.basic)):
                    if i == violating_basic_idx:
                        # this is the pivot row, already handled
                        continue
                    
                    old_coeff = self.tableaux[i][selected_nb_idx]
                    if old_coeff == 0:
                        # this row doesn't depend on the variable we're pivoting, skip
                        continue
                    
                    # substitute: nb[selected] = new_pivot_row expression
                    # basic[i] = ... + old_coeff * (new_pivot_row) + ...
                    for j in range(len(self.non_basic)):
                        if j == selected_nb_idx:
                            # this will become the coefficient for old_basic
                            self.tableaux[i][j] = old_coeff * new_pivot_row[j]
                        else:
                            # add the contribution from the substitution
                            self.tableaux[i][j] = self.tableaux[i][j] + old_coeff * new_pivot_row[j]
                
                # finally, update the pivot row and swap variables
                self.tableaux[violating_basic_idx] = new_pivot_row
                self.basic[violating_basic_idx] = old_nb
                self.non_basic[selected_nb_idx] = old_basic

        self.check_invariant()

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
