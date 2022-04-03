import numpy as np
from scipy.sparse import csc_matrix
from scipy.sparse.linalg import spsolve


class PowerFlow:

    def __init__(self,
                 case_name,
                 tol,
                 max_iters,
                 enable_limiting):
        """Initialize the PowerFlow instance.

        Args:
            case_name (str): A string with the path to the test case.
            tol (float): The chosen NR tolerance.
            max_iters (int): The maximum number of NR iterations.
            enable_limiting (bool): A flag that indicates if we use voltage limiting or not in our solver.
        """
        # Clean up the case name string
        case_name = case_name.replace('.RAW', '')
        case_name = case_name.replace('testcases/', '')

        self.case_name = case_name
        self.tol = tol
        self.max_iters = max_iters
        self.enable_limiting = enable_limiting

    def solve(self, Y, J):
        return spsolve(Y, J)

    def apply_limiting(self):
        pass

    def check_error(self, v, v_sol):
        return np.amax(np.abs(v - v_sol))

    def stamp_linear(self, branch, transformer, shunt, slack, v_init):
        size_Y = v_init.shape[0]
        nnz = 100*size_Y
        Ylin_row = np.zeros(nnz, dtype=int)
        Ylin_col = np.zeros(nnz, dtype=int)
        Ylin_val = np.zeros(nnz, dtype=np.double)
        Jlin_row = np.zeros(4*size_Y, dtype=int)
        Jlin_val = np.zeros(4*size_Y, dtype=np.double)
        
        idx_Y = 0
        idx_J = 0

        for ele in branch:
            (idx_Y, idx_J) = ele.stamp(v_init, Ylin_val, Ylin_row, Ylin_col, Jlin_val, Jlin_row, idx_Y, idx_J)
        for ele in transformer:
            (idx_Y, idx_J) = ele.stamp(v_init, Ylin_val, Ylin_row, Ylin_col, Jlin_val, Jlin_row, idx_Y, idx_J)
        for ele in shunt:
            (idx_Y, idx_J) = ele.stamp(v_init, Ylin_val, Ylin_row, Ylin_col, Jlin_val, Jlin_row, idx_Y, idx_J)
        for ele in slack:
            (idx_Y, idx_J) = ele.stamp(v_init, Ylin_val, Ylin_row, Ylin_col, Jlin_val, Jlin_row, idx_Y, idx_J)
        
        nnz_indices = np.nonzero(Ylin_val)[0]
        Ylin = csc_matrix((Ylin_val[nnz_indices], (Ylin_row[nnz_indices], Ylin_col[nnz_indices])), shape=(size_Y, size_Y), dtype=np.float64)
        nnz_indices = np.nonzero(Jlin_val)[0]
        Jlin_col = np.zeros(Jlin_row.shape, dtype=np.int)
        Jlin = csc_matrix((Jlin_val, (Jlin_row, Jlin_col)), shape=(size_Y, 1), dtype=np.float64)
        return (Ylin, Jlin)


    def stamp_nonlinear(self, generator, load, v_init):
        size_Y = v_init.shape[0]
        nnz = 100*size_Y
        Ynlin_row = np.zeros(nnz, dtype=int)
        Ynlin_col = np.zeros(nnz, dtype=int)
        Ynlin_val = np.zeros(nnz, dtype=np.double)
        Jnlin_row = np.zeros(4*size_Y, dtype=int)
        Jnlin_val = np.zeros(4*size_Y, dtype=np.double)
        
        idx_Y = 0
        idx_J = 0

        for ele in generator:
            (idx_Y, idx_J) = ele.stamp(v_init, Ynlin_val, Ynlin_row, Ynlin_col, Jnlin_val, Jnlin_row, idx_Y, idx_J)
        for ele in load:
            (idx_Y, idx_J) = ele.stamp(v_init, Ynlin_val, Ynlin_row, Ynlin_col, Jnlin_val, Jnlin_row, idx_Y, idx_J)

        nnz_indices = np.nonzero(Ynlin_val)[0]
        Ynlin = csc_matrix((Ynlin_val[nnz_indices], (Ynlin_row[nnz_indices], Ynlin_col[nnz_indices])), shape=(size_Y, size_Y), dtype=np.float64)
        nnz_indices = np.nonzero(Jnlin_val)[0]
        Jlin_col = np.zeros(Jnlin_row.shape, dtype=np.int)
        Jnlin = csc_matrix((Jnlin_val, (Jnlin_row, Jlin_col)), shape=(size_Y, 1), dtype=np.float64)
        return (Ynlin, Jnlin)
    
    def calc_resid(self, v, generator, load, slack, branch, transformer, shunt):
        resid = np.zeros(v.shape)
        for ele in slack:
            ele.calc_residuals(resid, v)
        for ele in generator:
            ele.calc_residuals(resid, v)
        for ele in load:
            ele.calc_residuals(resid, v)
        for ele in branch:
            ele.calc_residuals(resid, v)
        for ele in transformer:
            ele.calc_residuals(resid, v)
        for ele in shunt:
            ele.calc_residuals(resid, v)
        return resid

    def run_powerflow(self,
                      v_init,
                      bus,
                      slack,
                      generator,
                      transformer,
                      branch,
                      shunt,
                      load):
        """Runs a positive sequence power flow using the Equivalent Circuit Formulation.

        Args:
            v_init (np.array): The initial solution vector which has the same number of rows as the Y matrix.
            bus (list): Contains all the buses in the network as instances of the Buses class.
            slack (list): Contains all the slack generators in the network as instances of the Slack class.
            generator (list): Contains all the generators in the network as instances of the Generators class.
            transformer (list): Contains all the transformers in the network as instance of the Transformers class.
            branch (list): Contains all the branches in the network as instances of the Branches class.
            shunt (list): Contains all the shunts in the network as instances of the Shunts class.
            load (list): Contains all the loads in the network as instances of the Load class.

        Returns:
            v(np.array): The final solution vector.

        """

        # # # Copy v_init into the Solution Vectors used during NR, v, and the final solution vector v_sol # # #
        v = np.copy(v_init)
        v_sol = np.copy(v)

        # # # Stamp Linear Power Grid Elements into Y matrix # # #
        # TODO: PART 1, STEP 2.1 - Complete the stamp_linear function which stamps all linear power grid elements.
        #  This function should call the stamp_linear function of each linear element and return an updated Y matrix.
        #  You need to decide the input arguments and return values.
        Ylin, Jlin = self.stamp_linear(branch, transformer, shunt, slack, v_init)

        # # # Initialize While Loop (NR) Variables # # #
        # TODO: PART 1, STEP 2.2 - Initialize the NR variables
        err_max = 10  # maximum error at the current NR iteration
        tol = self.tol  # chosen NR tolerance
        NR_count = 0  # current NR iteration

        # # # Begin Solving Via NR # # #
        # TODO: PART 1, STEP 2.3 - Complete the NR While Loop
        while err_max > tol:

            # # # Stamp Nonlinear Power Grid Elements into Y matrix # # #
            # TODO: PART 1, STEP 2.4 - Complete the stamp_nonlinear function which stamps all nonlinear power grid
            #  elements. This function should call the stamp_nonlinear function of each nonlinear element and return
            #  an updated Y matrix. You need to decide the input arguments and return values.
            (Ynlin, Jnlin) = self.stamp_nonlinear(generator, load, v)

            # # # Solve The System # # #
            # TODO: PART 1, STEP 2.5 - Complete the solve function which solves system of equations Yv = J. The
            #  function should return a new v_sol.
            #  You need to decide the input arguments and return values.
            Y = Ynlin + Ylin
            J = Jnlin + Jlin
            zero_rows = []
            zero_cols = []
            for i in range(Y.shape[0]):
                if len(Y[i,:].data) == 0:
                    zero_rows.append(i)
                if len(Y[:,i].data) == 0:
                    zero_cols.append(i)
            v_sol = self.solve(Y, J)
            NR_count += 1
            # # # Compute The Error at the current NR iteration # # #
            # TODO: PART 1, STEP 2.6 - Finish the check_error function which calculates the maximum error, err_max
            #  You need to decide the input arguments and return values.
            err_max = self.check_error(v, v_sol)
            print("Iter: %d, max error: %.3e" % (NR_count, err_max))
            print(np.argmax(np.abs(v-v_sol)))
            # # # Compute The Error at the current NR iteration # # #
            # TODO: PART 2, STEP 1 - Develop the apply_limiting function which implements voltage and reactive power
            #  limiting. Also, complete the else condition. Do not complete this step until you've finished Part 1.
            #  You need to decide the input arguments and return values.
            if self.enable_limiting and err_max > tol:
                v = self.apply_limiting(v, v_sol)
            else:
                v = np.copy(v_sol)

        resid = self.calc_resid(v_sol, generator, load, slack, branch, transformer, shunt)
        max_resid = np.amax(np.abs(resid))
        max_resid_ind = np.argmax(np.abs(resid))
        print("Powerflow converged in %d iterations" % (NR_count))
        print("Maximum residual in system is %.3e" % (max_resid))

        return v
