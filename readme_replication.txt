This file provides step by step instructions to replicate all results in the paper, as submitted for publication to the Journal of Financial Economics.




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        Part 1: Solve and Simulate the Model
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

This part can be done locally or on a HPCC cluster. Instructions below are for local execution. If recomputing all economies (=paramater combinations) required for the paper, it is much faster to do it on the cluster.


------------------------------
Part 1a: Computing locally
------------------------------

Repeat the steps below for each economy in experdef_20200904.txt. The definitions of these economies are in experdef_20200904.m. Let the placeholder "econ" represent the name of the economy from this text file.

First, run the economy for up to 100 iterations on the coarse grid.
1.  Configure main_create_env.m as follows (leave other variables as is):
	- expername = 'econ_ini0';
	- guess_mode = 'no_guess';
2.  Run main_create_env.m and it will produce a file called env_econ_ini0.mat
3.  Configure main_run_exper.m as follows (leave other variables as is):
	- no_par_processes = XX; % where XX is the number of processors on the machine you're running it on
	- exper_path = 'env_econ_ini0.mat';
	- maxit = 100;
4.  Run main_run_exper.m. On a machine with 16 cores, this should take about 30 min. It will create a file named res_TIMESTAMP.mat, where TIMESTAMP is the time at which the calcuation is finished expressed as YYYY_MM_DD_hh_mm
5.  Rename the res_TIMESTAMP.mat file to res_20200904_econ_ini0.mat.

Next, run the economy for up to 100 iterations on the fine grid
6.  Configure main_create_env.m as follows (leave other variables as is):
	- expername = 'econ';
	- guess_mode = 'guess';
	- guess_path = 'res_20200910_econ_i100';
7.  Run main_create_env.m and it will produce a file called env_econ.mat.
8.  Configure main_run_exper.m as follows (leave other variables as is):
	- no_par_processes = XX; % where XX is the number of processors on the machine you're running it on
	- exper_path = 'env_econ.mat';
	- maxit = 100;
9.  Run main_run_exper.m. On a machine with 16 cores, this should take about 20 min. It will create a file named res_TIMESTAMP.mat, where TIMESTAMP is the time at which the calcuation is finished expressed as YYYY_MM_DD_hh_mm
10. Rename the res_TIMESTAMP.mat file to res_20200904_econ.mat.

Next, simulate. The model must be solved and res* file must exist.
11. Configure sim_stationary.m as follows (leave other variables as is):
	- resfile = 'res_20200904_econ';
12. Run sim_stationary.m. If you are running sim_stationary having run it before during the current MATLAB session, run "clear" beforehand. This operation will create the following files:
	- sim_res_20200904_econ.mat: all simulation results incl. full series, statistics, errors, and parameters
	- Results/statsexog_res_20200904_econ.xls: statistics using exogenous subsampling to define crises (one sample per worksheet)
	- Results/errstats_res_20200904_econ.xls: statistics of EE, VF, and TF errors

Next, compute impulse response functions. Model must be solved and simulated. Both res* and sim_res* files must exist.
13. Configure sim_trans_irf.m as follows (leave other variables as is):
	- resfile = 'res_20200904_econ';	
14. Run sim_trans_irf.m. On a machine with 16 cores, this should take about 5 min. It will create the following file:
	- GR_res_20200904_econ.mat: mean, median, and sd of IRF paths for each of 3 shocks (no shock, non-housing rec, housing rec)
	
Next, compute IRFs for the MIT shock to bank capital. Model must be solved and simulated. Both res* and sim_res* files must exist.
13. Configure sim_trans_mit.m as follows (leave other variables as is):
	- resfile = 'res_20200904_econ';	
14. Run sim_trans_mit.m. On a machine with 16 cores, this should take about 5 min. It will create the following file:
	- ST_20200904_MITloss.mat: mean, median, and sd of IRF paths 
	
Next, compute and simulate the transition paths for the boom-bust simulations. Model must be solved and simulated. Both res* and sim_res* files must exist. These instructions are for the baseline path caused by higher demand for assets. 
13. Configure main_run_trend.m as follows (leave other variables as is):
	- res_path='res_20200904_lowratefinal';   % this is the file containing the policy functions of the final economy
	- experseq = { 'trans01_eta';  % these are the experiment definitions in "experdef_20200904" containing the intermediate steps
              'trans02_eta';
              'trans03_eta';
              'trans04_eta';
              'trans05_eta';
              'trans06_eta';
              'trans07_eta';
              'trans08_eta'};
	- npp=XX;  % number of cores
14. Run main_run_trend.m. This will produce as many res_* files as there are intermediate steps in experseq. Should take about 60 minutes on 16 cores. 
15. Configure sim_trans_trend.m as follows:
	- experdef='20200904';
	- experinit = ['res_',experdef,'_bench'];
	- experfinal_list = {['res_',experdef,'_lowratefinal']};
	- translabel='_eta';
	- translabelfin='';
	- res_names = {['res_',experdef,'_trans01',translabel],...
             ['res_',experdef,'_trans02',translabel],...
             ['res_',experdef,'_trans03',translabel],...
             ['res_',experdef,'_trans04',translabel],...
             ['res_',experdef,'_trans05',translabel],...
             ['res_',experdef,'_trans06',translabel],...
             ['res_',experdef,'_trans07',translabel],...
             ['res_',experdef,'_trans08',translabel],...
             ['res_',experdef,'_lowratefinal',translabelfin]};
	Run sim_trans_trend.m. On a machine with 16 cores, this should take about 1 min. It will create the following file:
	- ST_20200904_eta_simtrend: mean, median, and sd of transition paths 	
	 
------------------------------------
Part 1b: Computing on a cluster
------------------------------------

The configuration needed to run the numerical experiments as batch job on a cluster will vary substantially based on the specific cluster configuration. However, the basic structure should be transferable.  




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
			Part 2: Generating the graphs and tables
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

This part pre-supposes that the res*, sim_res*, and other results files described in part 1 have been computed.

The steps below re-create the key results used in the paper.

1.  Model numbers in Tables 2 and 3: these numbers are taken from the simulation output of the baseline economy (steps 11 and 12 above). They are contained in Results/statsexog_res_20200904_econ.xls.
	
2.  Create benchmark economy IRFs (Figure 4). Configure plot_trans_irf.m as follows:
	- resfile='res_20200904_bench';
	Run plot_trans_irf.m. This will create figures in Results/GR_res_20200904_bench_IRF#.(pdf|eps).

3. Create a figure showing the IRFs after an unanticipated shock to intermediary wealth (Figure 5).  Configure plot_trans_mit.m as follows:
	- resfiles = {'20200904'};
	- res_suffix='_MITloss';
	- plot_exper=1;
	Run plot_trans_mit.m. This will create figures in Results/transitions_20200904_MITloss_#.(pdf|eps).


4.  Create figure with transitions from benchmark to boom-bust (Figure 6). This requires that the transitions for the baseline case with only increased asset demand (ST_20200904_eta_simtrend.mat), and the two extensions with "securitization" (ST_20200904_secur_simtrend.mat) and optimistic beliefs about credit risk (ST_20200904_belief_simtrend.mat). Given these solutions, run plot_trans_trend.m as configured in the repository. This will produce figures in Results/transitions_20200904_belief_simtrend_#.(pdf|eps).



