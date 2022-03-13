from scripts.Solve import solve
from models.global_vars import global_vars

# path to the grid network RAW file
casename = 'testcases/IEEE-14_prior_solution.RAW'

# the settings for the solver
settings = {
    "Tolerance": 1E-05,
    "Max Iters": 1000,
    "Limiting":  False
}

# run the solver
solve(casename, settings)