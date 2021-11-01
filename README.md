# Code for "Credit Cycles with Market-Based Household Leverage"
## William Diamond and Tim Landvoigt<sup>1</sup>


Fast Steps to Reproduce Benchmark Model
=======================================

-   In Matlab, execute the script `main_create_env`.

     -  `main_create_env` will create a file `env_bench_ini0` that
        contains the experiment definition for the benchmark economy on
        the initial coarse grid.

        -   `expername = 'bench_ini0';`

        -   `guess_mode = 'no_guess';`

-   Execute script `main_run_exper` to run the benchmark on the coarse
    grid for up to 100 iterations.

    -   You can set the number of parallel workers to be started.

    -   Set to zero if you want to run it with a single process.

    -   On a computer with sixteen cores (and 16 parallel workers) this
        should take about 30 minutes.

    -   `main_run_exper` creates a results file named
        `res_[current_date_time]` that contains the converged policy
        functions.

    -   Rename this file to `res_20200904_bench_ini0.mat`.

-   Execute `main_create_env` again using results file from previous
    step as input.

    -   Configure `main_create_env` as follows (leave other variables as
        is):

        -   `expername = 'bench';`

        -   `guess_mode = 'guess';`

        -   `guess_path = 'res_20200904_bench_ini0';`

    -   `main_create_env` will create a file `env_bench` that contains
        the experiment definition for the benchmark economy on the fine
        grid, using the resulting policy functions from the coarse grid
        iterations in `res_20200904_bench_ini0.mat` as initial guess.

-   Execute script `main_run_exper` to run the benchmark on the fine
    grid for up to 100 iterations.

    -   Configure `main_run_exper` as follows (leave other variables as
        is):

        -   `exper_path = 'env_bench.mat';`

    -   On a computer with sixteen cores (and 16 parallel workers) this
        should take about 45 minutes.

    -   `main_run_exper` creates a results file named
        "`res_[current_date_time]`" that contains the converged policy
        functions.

    -   Rename this file to `res_20200904_bench.mat`.

-   Simulate the model using `sim_stationary` and `sim_trans_irf`.

    -   `sim_stationary` simulates the model policies contained in
        `res_20200904_bench.mat` for 5,000 periods and writes out
        the resulting time-series and several statistics. The main
        output is a file named `sim_res_20200904_bench.mat`.

    -   `sim_trans_irf` reads both `res_20200904_bench.mat` and
        `sim_res_20200904_bench.mat`, and simulates generalized
        IRFs.

    -   To plot IRFs, run `plot_trans_irf`.
	
**For More Details See readme_replication.txt**
	

<sup>1</sup>: Diamond: University of Pennsylvania Wharton
    School; email:
    <diamondw@wharton.upenn.edu>. Landvoigt: University of Pennsylvania Wharton
    School; email: <timland@wharton.upenn.edu>.
