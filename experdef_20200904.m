params=struct;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% mean-reverting AR(1) growth rate %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

params.sig_y=.027;
params.rho_y=.45;
params.mu_y=-.5*params.sig_y^2/(1-params.rho_y^2);
exp_g =.02;
params.mu_g=exp_g;
params.mu_G=exp(params.mu_g);

%%%%%%%%%%%%%%%
% Preferences %
%%%%%%%%%%%%%%%

params.gamma=1;
params.beta=.925;
params.theta=0.135;
params.psiY=0.015;
% Liquidity preference M
params.psiM=0.015;

%%%%%%%%%%%%%%
% Life-cycle %
%%%%%%%%%%%%%%

params.nu=0.45;
params.depY=true; % What is this?
params.piM=0.05;
params.piY=0.05;
params.phi=2;
params.delta_eta=0.94;

%%%%%%%%%%%%%%%%%%%%%%%%%
% Housing and mortgages %
%%%%%%%%%%%%%%%%%%%%%%%%%

params.deltaH=.025; % Maintenance cost.
params.betaO=1;


% hard borrowing constraints?
params.LTV=false;
% LTV Parameters
params.thetaM_LTV=.45;
params.thetaY_LTV=.68;

params.PTI=false;
% PTI Parameters
params.thetaM_PTI=0.83;
params.thetaY_PTI=1.14;

% rental market
params.rental=false;
params.kappaY = .1;

% Mortgages and default
params.Hbar=1;

% Default cost for each agent
params.lambdaY=0.01;
params.lambdaM=0.01;

% Foreclosure loss to bank
params.xi=0.35;
params.xiDWL=0; 

% Values for the housing risk shock
params.sig_epsYlow = .19;
params.sig_epsYhigh = .32;
params.sig_epsMlow = .19;
params.sig_epsMhigh = .32;
params.sig_epsY=[params.sig_epsYlow,params.sig_epsYhigh];
params.sig_epsM=[params.sig_epsMlow,params.sig_epsMhigh];

% Values for the transition matrix of the housing risk shock
Epsprob_recess=[0.95, 0.05;
                0.2, 0.8];
            
params.Epsprob_recess=Epsprob_recess;  

Epsprob_boom=[0.95, 0.05;
                0.2, 0.8];
            
params.Epsprob_boom=Epsprob_boom;   

%%%%%%%%%%%%%%%%
% Intermediary %
%%%%%%%%%%%%%%%%
params.chi=850;
params.tau=0.068;
params.ebarY=0.08; % for steady state
params.ebarM=0.01; % for steady state
params.zeta=0.014; 
params.tax=0.012;
params.alpha=0;
params.MITshockloss=0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Grids for State Variables %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%
% Initial grid %
%%%%%%%%%%%%%%%%

startgrid.label = 'ini0';
startgrid.wMpts=linspace(1.8,2.5,6);
startgrid.eIpts=linspace(0.03, 0.18,6);

finegrid.label = '';
finegrid.wMpts=[1.95,2.0,2.05,2.1,2.12, 2.14,2.16,2.18,2.2,2.25,2.3,2.4];
finegrid.eIpts=linspace(0.035, 0.12,9);

benchgrid=struct('startgrid', startgrid, 'finegrid', finegrid);

%%%%%%%%%%%%%%
% Final grid %
%%%%%%%%%%%%%%

lowrategrid=startgrid;

lowrategrid.wMpts=linspace(2.1,2.9,6);
lowrategrid.eIpts=linspace(0.01, 0.16,6);
lowratefinegrid=finegrid;
lowratefinegrid.wMpts=linspace(2.1,2.9,12);
lowratefinegrid.eIpts=linspace(0.01, 0.16,10);

beliefgrid=startgrid;
beliefgrid.wMpts=linspace(2.4,3.15,6);
beliefgrid.eIpts=linspace(0.01, 0.13,6);
belieffinegrid=finegrid;
belieffinegrid.wMpts=linspace(2.4,3.15,12);
belieffinegrid.eIpts=linspace(0.01, 0.13,9);

lowrategridrental=startgrid;
lowrategridrental.wMpts=linspace(2.4,3.2,6);
lowrategridrental.eIpts=linspace(0.01, 0.13,6);
rentalfinegrid=finegrid;
rentalfinegrid.wMpts=linspace(2.4,3.2,12);
rentalfinegrid.eIpts=linspace(0.01, 0.13,10);


lowrategrid=struct('startgrid', lowrategrid, 'finegrid', lowratefinegrid);
beliefgrid=struct('startgrid', beliefgrid, 'finegrid', belieffinegrid);
lowrategridrental=struct('startgrid', lowrategridrental, 'finegrid', rentalfinegrid);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Transition grid: Fixed grid method %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

medrategrid=startgrid;
medrategrid.wMpts=linspace(1.8,2.95,15);
medrategrid.eIpts=linspace(0.02,0.18,12);
medrategrid=struct('startgrid', medrategrid, 'finegrid', finegrid);

% grid for experiment with high cost of default
hilamgrid=startgrid;
hilamgrid.wMpts=linspace(1.5,2.4,6);
hilamgrid.eIpts=linspace(0.03, 0.25,6);
hilamfinegrid=finegrid;
hilamfinegrid.wMpts=linspace(1.5,2.4,12);
hilamfinegrid.eIpts=linspace(0.03, 0.25,9);
hilamgrid=struct('startgrid', hilamgrid, 'finegrid', hilamfinegrid);


% Put all grids together.
grids = struct('benchgrid',benchgrid, ...
    'lowrategrid',lowrategrid,...
    'beliefgrid',beliefgrid,...
    'lowrategridrental',lowrategridrental,...
    'medrategrid',medrategrid,...
    'hilamgrid',hilamgrid);

%%%%%%%%%%%%%%%%%%%%%
% EXPERIMENT VALUES %
%%%%%%%%%%%%%%%%%%%%%

% Grids for the experiments
% Base
transpath_piM=linspace(0.04,0.04,10)';
transpath_psiM=linspace(0.015,0.04,10)';
transpath_eta=linspace(0.94,0.99,10)';

% Deregulation of Banks
transpath_chi=linspace(850,400,10)';
transpath_ebarY=linspace(.08,.04,10)';
transpath_ebarM=linspace(.01,0,10)';
transpath_kappaY=linspace(.17,0.2,10)';

% Housing beliefs
transpath_epsYhigh=linspace(.32,.25,9)';
transpath_epsMhigh=linspace(.32,.25,9)';

% Define experiments as modifications to base
Epsprob_exper=[0.98, 0.02;
                0.2, 0.8];
            
Epsprob_exper01=[0.99, 0.01;
                0.2, 0.8];
  

          
expers_calib = {'bench',{};
                'highgamma',{'gamma',3; 'beta',0.925; 'grid', 'hilamgrid'};
                'bench_LTV',{'LTV',true};
                'bench_DTI',{'PTI',true};
                'bench_rental',{'rental',true; 'kappaY',.17; 'theta',0.115;};
                'lowratefinal',{'delta_eta',.99; 'psiM', .04; 'grid','lowrategrid'};
                
                'lowratefinal_LTV',{'delta_eta',.99; 'psiM', .04; 'LTV',true; 'grid','lowrategrid'};
                'lowratefinal_DTI',{'delta_eta',.99; 'psiM', .04; 'PTI',true; 'grid','lowrategrid'};
                'lowratefinal_rental',{'delta_eta',.99; 'psiM', .04; 'rental',true; 'kappaY',0.2; 'theta',0.115; 'grid','lowrategrid'};
                'lowratefinal_bank',{'delta_eta',.99; 'psiM', .025; 'ebarY', .04; 'ebarM', 0; 'grid','lowrategrid'};
                'lowratefinal_secur',{'delta_eta',.99; 'psiM', .04; 'ebarY', .04; 'ebarM', 0; 'chi',400;  'grid','lowrategrid'};
                'lowratefinal_belief',{'delta_eta',.99; 'psiM', .04; 'ebarY', .04; 'ebarM', 0; 'chi',400; 'sig_epsYhigh',0.25; 'sig_epsMhigh',0.25; 'grid','beliefgrid'};
                'lowratefinal_secur_LTV',{'delta_eta',.99; 'psiM', .04; 'ebarY', .04; 'ebarM', 0; 'chi',400; 'LTV', true; 'grid','lowrategrid'};
                'lowratefinal_secur_DTI',{'delta_eta',.99; 'psiM', .04; 'ebarY', .04; 'ebarM', 0; 'chi',400; 'PTI', true; 'grid','lowrategrid'};
                'lowratefinal_secur_rental',{'delta_eta',.99; 'psiM', .04; 'ebarY', .04; 'ebarM', 0; 'chi',400; 'rental', true; 'kappaY',0.2; 'theta',0.115; 'grid','lowrategrid'};
                                           
               
                'trans01_eta',{'delta_eta',transpath_eta(2);'psiM',transpath_psiM(2); 'grid','lowrategrid'};
                'trans02_eta',{'delta_eta',transpath_eta(3);'psiM',transpath_psiM(3); 'grid','lowrategrid'};
                'trans03_eta',{'delta_eta',transpath_eta(4);'psiM',transpath_psiM(4); 'grid','lowrategrid'};
                'trans04_eta',{'delta_eta',transpath_eta(5);'psiM',transpath_psiM(5); 'grid','lowrategrid'};
                'trans05_eta',{'delta_eta',transpath_eta(6);'psiM',transpath_psiM(6); 'grid','lowrategrid'};
                'trans06_eta',{'delta_eta',transpath_eta(7);'psiM',transpath_psiM(7); 'grid','lowrategrid'};
                'trans07_eta',{'delta_eta',transpath_eta(8);'psiM',transpath_psiM(8); 'grid','lowrategrid'};
                'trans08_eta',{'delta_eta',transpath_eta(9);'psiM',transpath_psiM(9); 'grid','lowrategrid'};
 
               'trans01_eta_LTV',{'delta_eta',transpath_eta(2);'psiM',transpath_psiM(2); 'LTV', true; 'grid','medrategrid'};
               'trans02_eta_LTV',{'delta_eta',transpath_eta(3);'psiM',transpath_psiM(3); 'LTV', true; 'grid','medrategrid'};
               'trans03_eta_LTV',{'delta_eta',transpath_eta(4);'psiM',transpath_psiM(4); 'LTV', true; 'grid','medrategrid'};
               'trans04_eta_LTV',{'delta_eta',transpath_eta(5);'psiM',transpath_psiM(5); 'LTV', true; 'grid','medrategrid'};
               'trans05_eta_LTV',{'delta_eta',transpath_eta(6);'psiM',transpath_psiM(6); 'LTV', true; 'grid','medrategrid'};
               'trans06_eta_LTV',{'delta_eta',transpath_eta(7);'psiM',transpath_psiM(7); 'LTV', true; 'grid','medrategrid'};
               'trans07_eta_LTV',{'delta_eta',transpath_eta(8);'psiM',transpath_psiM(8); 'LTV', true; 'grid','medrategrid'};
               'trans08_eta_LTV',{'delta_eta',transpath_eta(9);'psiM',transpath_psiM(9); 'LTV', true; 'grid','medrategrid'};
               
                            
               'trans01_secur',{'delta_eta',transpath_eta(2);'psiM',transpath_psiM(2); 'ebarY', transpath_ebarY(2,1); 'ebarM', transpath_ebarM(2,1); 'chi', transpath_chi(2,1); 'grid','lowrategrid'};
               'trans02_secur',{'delta_eta',transpath_eta(3);'psiM',transpath_psiM(3); 'ebarY', transpath_ebarY(3,1); 'ebarM', transpath_ebarM(3,1); 'chi', transpath_chi(3,1); 'grid','lowrategrid'};
               'trans03_secur',{'delta_eta',transpath_eta(4);'psiM',transpath_psiM(4); 'ebarY', transpath_ebarY(4,1); 'ebarM', transpath_ebarM(4,1); 'chi', transpath_chi(4,1); 'grid','lowrategrid'};
               'trans04_secur',{'delta_eta',transpath_eta(5);'psiM',transpath_psiM(5); 'ebarY', transpath_ebarY(5,1); 'ebarM', transpath_ebarM(5,1); 'chi', transpath_chi(5,1); 'grid','lowrategrid'};
               'trans05_secur',{'delta_eta',transpath_eta(6);'psiM',transpath_psiM(6); 'ebarY', transpath_ebarY(6,1); 'ebarM', transpath_ebarM(6,1); 'chi', transpath_chi(6,1); 'grid','lowrategrid'};
               'trans06_secur',{'delta_eta',transpath_eta(7);'psiM',transpath_psiM(7); 'ebarY', transpath_ebarY(7,1); 'ebarM', transpath_ebarM(7,1); 'chi', transpath_chi(7,1); 'grid','lowrategrid'};
               'trans07_secur',{'delta_eta',transpath_eta(8);'psiM',transpath_psiM(8); 'ebarY', transpath_ebarY(8,1); 'ebarM', transpath_ebarM(8,1); 'chi', transpath_chi(8,1); 'grid','lowrategrid'};
               'trans08_secur',{'delta_eta',transpath_eta(9);'psiM',transpath_psiM(9); 'ebarY', transpath_ebarY(9,1); 'ebarM', transpath_ebarM(9,1); 'chi', transpath_chi(9,1); 'grid','lowrategrid'};
 
               'trans01_belief',{'delta_eta',transpath_eta(2);'psiM',transpath_psiM(2); ...
                                       'ebarY', transpath_ebarY(2,1); 'ebarM', transpath_ebarM(2,1); 'chi', transpath_chi(2,1);...
                                       'sig_epsYhigh',transpath_epsYhigh(2,1); 'sig_epsMhigh',transpath_epsMhigh(2,1); 'grid','lowrategrid'};
               'trans02_belief',{'delta_eta',transpath_eta(3);'psiM',transpath_psiM(3); ...
                                       'ebarY', transpath_ebarY(3,1); 'ebarM', transpath_ebarM(3,1); 'chi', transpath_chi(3,1);...
                                       'sig_epsYhigh',transpath_epsYhigh(3,1); 'sig_epsMhigh',transpath_epsMhigh(3,1); 'grid','lowrategrid'};
               'trans03_belief',{'delta_eta',transpath_eta(4);'psiM',transpath_psiM(4);...
                                       'ebarY', transpath_ebarY(4,1); 'ebarM', transpath_ebarM(4,1); 'chi', transpath_chi(4,1);...
                                       'sig_epsYhigh',transpath_epsYhigh(4,1); 'sig_epsMhigh',transpath_epsMhigh(4,1); 'grid','lowrategrid'};
               'trans04_belief',{'delta_eta',transpath_eta(5);'psiM',transpath_psiM(5);...
                                       'ebarY', transpath_ebarY(5,1); 'ebarM', transpath_ebarM(5,1); 'chi', transpath_chi(5,1);...
                                       'sig_epsYhigh',transpath_epsYhigh(5,1); 'sig_epsMhigh',transpath_epsMhigh(5,1); 'grid','lowrategrid'};
               'trans05_belief',{'delta_eta',transpath_eta(6);'psiM',transpath_psiM(6); ...
                                       'ebarY', transpath_ebarY(6,1); 'ebarM', transpath_ebarM(6,1); 'chi', transpath_chi(6,1);...
                                       'sig_epsYhigh',transpath_epsYhigh(6,1); 'sig_epsMhigh',transpath_epsMhigh(6,1); 'grid','lowrategrid'};
               'trans06_belief',{'delta_eta',transpath_eta(7);'psiM',transpath_psiM(7);...
                                       'ebarY', transpath_ebarY(7,1); 'ebarM', transpath_ebarM(7,1); 'chi', transpath_chi(7,1);...
                                       'sig_epsYhigh',transpath_epsYhigh(7,1); 'sig_epsMhigh',transpath_epsMhigh(7,1); 'grid','beliefgrid'};
               'trans07_belief',{'delta_eta',transpath_eta(8);'psiM',transpath_psiM(8);...
                                       'ebarY', transpath_ebarY(8,1); 'ebarM', transpath_ebarM(8,1); 'chi', transpath_chi(8,1);...
                                       'sig_epsYhigh',transpath_epsYhigh(8,1); 'sig_epsMhigh',transpath_epsMhigh(8,1); 'grid','beliefgrid'};
               'trans08_belief',{'delta_eta',transpath_eta(9);'psiM',transpath_psiM(9);...
                                       'ebarY', transpath_ebarY(9,1); 'ebarM', transpath_ebarM(9,1); 'chi', transpath_chi(9,1);...
                                       'sig_epsYhigh',transpath_epsYhigh(9,1); 'sig_epsMhigh',transpath_epsMhigh(9,1); 'grid','beliefgrid'};

%                'trans01_secur_LTV',{'delta_eta',transpath_eta(2);'psiM',transpath_psiM(2); ...
%                                        'ebarY', transpath_ebarY(2,1); 'ebarM', transpath_ebarM(2,1); 'chi', transpath_chi(2,1); 'LTV',true; 'grid','lowrategrid'};
%                'trans02_secur_LTV',{'delta_eta',transpath_eta(3);'psiM',transpath_psiM(3); ...
%                                        'ebarY', transpath_ebarY(3,1); 'ebarM', transpath_ebarM(3,1); 'chi', transpath_chi(3,1);'LTV',true; 'grid','lowrategrid'};
%                'trans03_secur_LTV',{'delta_eta',transpath_eta(4);'psiM',transpath_psiM(4);...
%                                        'ebarY', transpath_ebarY(4,1); 'ebarM', transpath_ebarM(4,1); 'chi', transpath_chi(4,1); 'LTV',true;'grid','lowrategrid'};
%                'trans04_secur_LTV',{'delta_eta',transpath_eta(5);'psiM',transpath_psiM(5);...
%                                        'ebarY', transpath_ebarY(5,1); 'ebarM', transpath_ebarM(5,1); 'chi', transpath_chi(5,1);'LTV',true; 'grid','lowrategrid'};
%                'trans05_secur_LTV',{'delta_eta',transpath_eta(6);'psiM',transpath_psiM(6); ...
%                                        'ebarY', transpath_ebarY(6,1); 'ebarM', transpath_ebarM(6,1); 'chi', transpath_chi(6,1);'LTV',true; 'grid','lowrategrid'};
%                'trans06_secur_LTV',{'delta_eta',transpath_eta(7);'psiM',transpath_psiM(7);...
%                                        'ebarY', transpath_ebarY(7,1); 'ebarM', transpath_ebarM(7,1); 'chi', transpath_chi(7,1); 'LTV',true;'grid','lowrategrid'};
%                'trans07_secur_LTV',{'delta_eta',transpath_eta(8);'psiM',transpath_psiM(8);...
%                                        'ebarY', transpath_ebarY(8,1); 'ebarM', transpath_ebarM(8,1); 'chi', transpath_chi(8,1);'LTV',true; 'grid','lowrategrid'};
%                'trans08_secur_LTV',{'delta_eta',transpath_eta(9);'psiM',transpath_psiM(9);...
%                                        'ebarY', transpath_ebarY(9,1); 'ebarM', transpath_ebarM(9,1); 'chi', transpath_chi(9,1);'LTV',true;'grid','lowrategrid'};
% 
%                'trans01_secur_DTI',{'delta_eta',transpath_eta(2);'psiM',transpath_psiM(2); ...
%                                        'ebarY', transpath_ebarY(2,1); 'ebarM', transpath_ebarM(2,1); 'chi', transpath_chi(2,1);'PTI',true; 'grid','lowrategrid'};
%                'trans02_secur_DTI',{'delta_eta',transpath_eta(3);'psiM',transpath_psiM(3); ...
%                                        'ebarY', transpath_ebarY(3,1); 'ebarM', transpath_ebarM(3,1); 'chi', transpath_chi(3,1); 'PTI',true; 'grid','lowrategrid'};
%                'trans03_secur_DTI',{'delta_eta',transpath_eta(4);'psiM',transpath_psiM(4);...
%                                        'ebarY', transpath_ebarY(4,1); 'ebarM', transpath_ebarM(4,1); 'chi', transpath_chi(4,1); 'PTI',true;'grid','lowrategrid'};
%                'trans04_secur_DTI',{'delta_eta',transpath_eta(5);'psiM',transpath_psiM(5);...
%                                        'ebarY', transpath_ebarY(5,1); 'ebarM', transpath_ebarM(5,1); 'chi', transpath_chi(5,1);'PTI',true; 'grid','lowrategrid'};
%                'trans05_secur_DTI',{'delta_eta',transpath_eta(6);'psiM',transpath_psiM(6); ...
%                                        'ebarY', transpath_ebarY(6,1); 'ebarM', transpath_ebarM(6,1); 'chi', transpath_chi(6,1);'PTI',true; 'grid','lowrategrid'};
%                'trans06_secur_DTI',{'delta_eta',transpath_eta(7);'psiM',transpath_psiM(7);...
%                                        'ebarY', transpath_ebarY(7,1); 'ebarM', transpath_ebarM(7,1); 'chi', transpath_chi(7,1); 'PTI',true;'grid','lowrategrid'};
%                'trans07_secur_DTI',{'delta_eta',transpath_eta(8);'psiM',transpath_psiM(8);...
%                                        'ebarY', transpath_ebarY(8,1); 'ebarM', transpath_ebarM(8,1); 'chi', transpath_chi(8,1);'PTI',true; 'grid','lowrategrid'};
%                'trans08_secur_DTI',{'delta_eta',transpath_eta(9);'psiM',transpath_psiM(9);...
%                                        'ebarY', transpath_ebarY(9,1); 'ebarM', transpath_ebarM(9,1); 'chi', transpath_chi(9,1);'PTI',true;'grid','lowrategrid'};
% 
%                'trans01_secur_rental',{'delta_eta',transpath_eta(2);'psiM',transpath_psiM(2); ...
%                                        'ebarY', transpath_ebarY(2,1); 'ebarM', transpath_ebarM(2,1); 'chi', transpath_chi(2,1); ...
%                                        'rental',true; 'kappaY',transpath_kappaY(2); 'theta',0.115; 'grid','lowrategrid'};
%                'trans02_secur_rental',{'delta_eta',transpath_eta(3);'psiM',transpath_psiM(3); ...
%                                        'ebarY', transpath_ebarY(3,1); 'ebarM', transpath_ebarM(3,1); 'chi', transpath_chi(3,1);...
%                                        'rental',true; 'kappaY',transpath_kappaY(3); 'theta',0.115;  'grid','lowrategrid'};
%                'trans03_secur_rental',{'delta_eta',transpath_eta(4);'psiM',transpath_psiM(4);...
%                                        'ebarY', transpath_ebarY(4,1); 'ebarM', transpath_ebarM(4,1); 'chi', transpath_chi(4,1);...
%                                        'rental',true; 'kappaY',transpath_kappaY(4); 'theta',0.115; 'grid','lowrategrid'};
%                'trans04_secur_rental',{'delta_eta',transpath_eta(5);'psiM',transpath_psiM(5);...
%                                        'ebarY', transpath_ebarY(5,1); 'ebarM', transpath_ebarM(5,1); 'chi', transpath_chi(5,1);...
%                                        'rental',true; 'kappaY',transpath_kappaY(5); 'theta',0.115; 'grid','lowrategrid'};
%                'trans05_secur_rental',{'delta_eta',transpath_eta(6);'psiM',transpath_psiM(6); ...
%                                        'ebarY', transpath_ebarY(6,1); 'ebarM', transpath_ebarM(6,1); 'chi', transpath_chi(6,1);...
%                                        'rental',true; 'kappaY',transpath_kappaY(6); 'theta',0.115; 'grid','lowrategrid'};
%                'trans06_secur_rental',{'delta_eta',transpath_eta(7);'psiM',transpath_psiM(7);...
%                                        'ebarY', transpath_ebarY(7,1); 'ebarM', transpath_ebarM(7,1); 'chi', transpath_chi(7,1);...
%                                        'rental',true; 'kappaY',transpath_kappaY(7); 'theta',0.115; 'grid','lowrategrid'};
%                'trans07_secur_rental',{'delta_eta',transpath_eta(8);'psiM',transpath_psiM(8);...
%                                        'ebarY', transpath_ebarY(8,1); 'ebarM', transpath_ebarM(8,1); 'chi', transpath_chi(8,1);...
%                                        'rental',true; 'kappaY',transpath_kappaY(8); 'theta',0.115; 'grid','lowrategrid'};
%                'trans08_secur_rental',{'delta_eta',transpath_eta(9);'psiM',transpath_psiM(9);...
%                                        'ebarY', transpath_ebarY(9,1); 'ebarM', transpath_ebarM(9,1); 'chi', transpath_chi(9,1);...
%                                        'rental',true; 'kappaY',transpath_kappaY(9); 'theta',0.115;'grid','lowrategrid'};
          };

      
% Parameter Sensitivity (Semi) Elasticities
vary_expers = {{'sig_epsYlow','sig_epsMlow'},{'sig_epsYhigh','sig_epsMhigh'},...
               'beta','theta',{'psiM','psiY'},'delta_eta','phi','tax', ...
	           {'lambdaY','lambdaM'}, 'zeta','chi'};
           
eps_default = 1e-3;
% Modify interval size for specific parameters
eps_struct = struct; 

% Compute full elasticities (i.e. percent parameter changes)
fullElasticities = {};


expers_gs=cell(0,6);
for ii=1:numel(vary_expers)
	currExper = vary_expers{ii};
	try
		currEps = eps_struct.(currExper);
	catch
		currEps = eps_default;
    end
    if ~iscell(currExper)
        currExper={currExper};
    end
   
    currExperdef = createGSExper(currExper,currEps,params,fullElasticities);
    
    expers_gs = [ expers_gs; currExperdef ];
end
      
      
      
% Second, vary other parameters
expers_other = { };
expers = [expers_calib; expers_other; expers_gs(:,1:2)];


%% Create each experiment definition (fully dynamic below this line)
N_EXP=size(expers,1);
N_GRIDS = length(fieldnames(benchgrid));

% Write list of defined experiments to file
fid = fopen([mfilename,'.txt'],'w');
for i=1:N_EXP
   fprintf(fid,expers{i,1});
   fprintf(fid,'\n');
end
fclose(fid);

% Create experiment table
% Basically is putting together the grids and the
% parameters of all the experiments and the benchmark model together.

baserow = {'test','ini0',{startgrid.wMpts},{startgrid.eIpts}};
basetable = cell2table(baserow);          
basetable.Properties.VariableNames = {'exper','grid','wMpts','eIpts'};
basetable = [basetable, struct2table(params,'AsArray',true)];

state_range=3:4;

expertable = repmat(basetable,N_EXP*N_GRIDS,1);

benchgrid_cell=struct2cell(benchgrid);

for i=1:N_GRIDS
    fnames = fieldnames(benchgrid_cell{i});
    gridvectors=struct2cell(benchgrid_cell{i})';
    [~,~,order]=intersect(expertable.Properties.VariableNames(state_range),fnames, ...
            'stable');
    expertable{N_EXP*(i-1)+1 : N_EXP*i,state_range} = ...
        repmat( gridvectors(order), N_EXP, 1);
    expertable(N_EXP*(i-1)+1 : N_EXP*i,2) = ...
        {benchgrid_cell{i}.label};   
end

gridOrder = @(grid_as_cell,order)grid_as_cell(order);
for i=1:N_EXP
   expertable.exper([i,N_EXP+i]) = expers(i,1);
   changes = expers{i,2};
   for j=1:size(changes,1)
      varname=changes{j,1};
      if strcmp(varname,'grid')
          newgrid = struct2cell(grids.(changes{j,2}));
          fnames = fieldnames(newgrid{1});
          [~,~,order]=intersect(expertable.Properties.VariableNames(state_range),fnames, ...
            'stable');
        for k=1:N_GRIDS
          expertable{i+N_EXP*(k-1),state_range} = ...
              gridOrder(struct2cell(newgrid{k}),order)';
        end
      else
          if ischar(changes{j,2}) || ~isscalar(changes{j,2})
            expertable.(varname)([i,N_EXP+i],:)=repmat(changes(j,2),N_GRIDS,1);
          else
            expertable.(varname)([i,N_EXP+i],:)=repmat(changes{j,2},N_GRIDS,1);    
          end
      end
   end
end

% Name each row
row_id = cell(2*N_EXP,1);
idx_suffix = find(~strcmp(expertable.grid,''));
idx_nosuffix = find(strcmp(expertable.grid,''));
row_id(idx_suffix) = arrayfun(@(idx){[expertable.exper{idx},'_',expertable.grid{idx}]},idx_suffix,'UniformOutput',false);
row_id(idx_nosuffix) = arrayfun(@(idx)expertable.exper(idx),idx_nosuffix,'UniformOutput',false);
expertable.Properties.RowNames = [row_id{:}]';
%expertable.Properties.RowNames(1) = {'xpbase'};

% Create a struct of structs
allexpers=cell2struct(expertable.Properties.RowNames,expertable.Properties.RowNames);
for i=1:size(expertable,1)
    s = table2struct(expertable(i,state_range));
    s.params = table2struct(expertable(i,state_range(end)+1:end));
    name = expertable.Properties.RowNames(i);
    name = name{:};
    allexpers.(name) = s;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function GSexper = createGSExper(currExper,currEps,params,fullElasticities)

    expername_U = sprintf('GS%sU',[currExper{1},'_']);
    expername_D = sprintf('GS%sD',[currExper{1},'_']);
	N_paramDim = length( currExper );	
	
    experdevs_U = {};
    experdevs_D = {};
	for jj=1:N_paramDim
       	currValue = params.(currExper{jj});

		newval_U = currValue;
		newval_D = currValue;
		if ismember(currExper,fullElasticities)
			newval_U=newval_U*exp(currEps);
			newval_D=newval_D*exp(-currEps);
		else
			newval_U=newval_U+currEps;
			newval_D=newval_D-currEps;
		end
		
		experdevs_U = [experdevs_U; {currExper{jj},newval_U}];
		experdevs_D = [experdevs_D; {currExper{jj},newval_D}];

    end

	GSexper = {expername_U,experdevs_U,currExper{1},'U'; 
               expername_D,experdevs_D,currExper{1},'D'};

end

