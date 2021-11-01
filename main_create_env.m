function [outfname,mobj]=main_create_env(varargin)

if nargin>0
    experdef_file=varargin{1};
    expername=varargin{2};
    guess_mode=varargin{3};
    if nargin>3
        guess_path=varargin{4};
        if nargin>4
            ststonly=varargin{5};
        end
    end    
end


% ===========================================
% program control
% ===========================================

% name of file with experiment definitions
if ~exist('experdef_file','var')
    experdef_file='experdef_20200904';
end
% name of experiment
if ~exist('expername','var')
%    expername='bench';
    expername='lowratefinal';
end 

% possible values 'no_guess', 'guess'
if ~exist('guess_mode','var')
    guess_mode='guess';
end

% path to file with initial guess; not used if guess_mode='no_guess'
if ~exist('guess_path','var')
    guess_path='res_20200904_lowratefinal_ini0.mat';
end

% path and file name for output file
outfname=['env_',expername,'.mat'];

% only compute steady state (without writing initial model object)
if ~exist('ststonly','var')
    ststonly = 0;
end

% ===================================================
% non-stochastic steady state
% ===================================================
warning('off','MATLAB:table:RowsAddedExistingVars');
% load experiment definitions
run(experdef_file);
disp(experdef_file);
% Look inside allexpers to see all the experiments (different parameter
% combbinations)
expdef=allexpers.(expername);

% transform for log-normality
sig_epsY = [expdef.params.sig_epsYlow, expdef.params.sig_epsYhigh];
sig_epsM = [expdef.params.sig_epsMlow, expdef.params.sig_epsMhigh];

expdef.params.sig2_epsY=log(1+sig_epsY.^2);
expdef.params.mu_epsY=-.5*expdef.params.sig2_epsY;
expdef.params.sig_epsY=sqrt(expdef.params.sig2_epsY);
expdef.params.sig2_epsM=log(1+sig_epsM.^2);
expdef.params.mu_epsM=-.5*expdef.params.sig2_epsM;
expdef.params.sig_epsM=sqrt(expdef.params.sig2_epsM);

modelclassname='OLGIntermediaryModel';
instfun=str2func(modelclassname);


guessfun=str2func([modelclassname,'.assignGuess']);

% I will name that function ststfun in this file.
ststfun=str2func([modelclassname,'.compStSt']);


% Variable names and order in the files.
% 1: log(P)
% 2: r 
% 3: log(q)
% 4: log(pY)
% 5: log(pM)
% 6: log(cY)
% 7: log(hY)
% 8: log(mY)
% 9: log(cM)
% 10: log(dM)
% 11: log(dY)
% 12: log(mM)
% 13: xI
% 14: log(BO)
% 15: log(vM)
% 16: log(vY)
% 17: log(Rebate)
% 18: muI
% 19. muI_rent

% Variable values to calculate the SS (initial guess).

% Initial economy
rental = expdef.params.rental;

if rental
gvec= [       0.8914
    0.0348
   -1.8343
    1.6893
    0.6558
   -0.7990
    0.3560
   -0.4673
   -0.8918
    0.2882
   -2.4635
   -0.1565
    0.0032
   -1.8438
    2.2529
    2.4849
   -5.0276
    0.3540
   -0.7087
];
else
gvec= [  
    0.7321
    0.0320
   -1.7384
    1.7326
    0.7577
   -0.7551
   -0.8428
   -0.4230
   -0.9141
    0.1272
   -2.4041
   -0.4569
    0.0032
   -1.8718
    2.2529
    2.4849
   -4.6603
    0.3603
    ];
end

options=optimset('Display','iter','TolX',1e-10,'TolFun',1e-10,'MaxIter',500);
options_step=optimset('Display','off','TolX',1e-10,'TolFun',1e-10,'MaxIter',500);

[nd,wgt]=model.HelperCollection.lgwt(30,-1,1);

%%%%%%%%%%%%%%
% Compute SS %
%%%%%%%%%%%%%%

% The function will read expdef which contains the particular
% parametrization of any experiment.
fh_compStSt=@(x)ststfun(x,expdef.params,0,nd,wgt);
[solvec,~,exfl]=fsolve(fh_compStSt,gvec,options);

if exfl<1
    disp('!! Problem computing steady state');
end
[~,stvtmp]=ststfun(solvec,expdef.params,1,nd,wgt);

stv=stvtmp;


if ststonly
    solvec
    return;
end


% ===================================================
% Parameters of stochastic model
% ===================================================

% for growth rate of income
N_y=3;
mu_y=expdef.params.mu_y;
sig_y=expdef.params.sig_y;
rho_y=expdef.params.rho_y;
Skew_y=0;
[Yprob,Y] = model.DSGEModel.rouwen(rho_y,mu_y,sig_y,Skew_y,N_y);
Yprob=Yprob';
Y=exp(Y);

% for shock to std.dev. of borrowers house values
N_eps=2;
EpsY=expdef.params.sig_epsY';
EpsM=expdef.params.sig_epsM';
Epsprob_recess=expdef.params.Epsprob_recess;
Epsprob_boom=expdef.params.Epsprob_boom;

% first find recession states
threshold=mu_y;
recess=find(Y<1+threshold,1,'last');
boom=N_y-recess;
epscomb=grid.StateSpaceGrid.makeCombinations([N_y,N_eps]);
mpts_perm=[ Y(epscomb(:,1)), EpsY(epscomb(:,2)), EpsM(epscomb(:,2))];

% transition matrix
rectrans=kron(Yprob(1:recess,:),Epsprob_recess);
boomtrans=kron(Yprob(recess+1:end,:),Epsprob_boom);
mtrans=[rectrans; boomtrans];
% Number of different combinations of the shock: in this case 6.
exnpt=size(mpts_perm,1);
% Putting the shocks and the last vector together.
mpts_all=[mpts_perm, (1:exnpt)'];
  
%%%%%%%%%%%%%%
% Simulation %
%%%%%%%%%%%%%%

% simulate exog process and produce histogram
Nsim=100000;
simst=model.DSGEModel.hitm_s(mtrans',rand(Nsim,1));
% Plotting an histogram of the values
% histogram(simst)
% Getting the fraction of the simulations in each state
simfrac=histc(simst,1:exnpt)/Nsim;
disp('-----------------------------------------');
disp('Unconditional prob of each state: ');
disp(num2str([(1:exnpt)',mpts_perm,simfrac]));
disp(['"Recession" states: ', num2str(sum(simfrac(1:N_eps*recess)))]); 
omhi_ind=N_eps*(1:N_y);
disp(['High uncertainty states: ', num2str(sum(simfrac(omhi_ind)))]);

% Average length of spell
forecl=ismember(simst,omhi_ind);
forestend=forecl(2:end)-forecl(1:end-1);
firstend=find(forestend==-1,1,'first');
% NO ES: laststart=find(forestend==-1,1,'last')?
laststart=find(forestend==1,1,'last');
forestend=forestend(firstend+1:laststart-1);
forecep=sum(forestend==1);
disp(['Uncertainty episodes (%): ',num2str(forecep/Nsim)]);
forest=find(forestend==1);
foreend=find(forestend==-1);
forelen=foreend-forest;
disp(['Average length: ',num2str(mean(forelen))]);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Exog. Var as a structure %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% variable lists
% exogenous state variable names
exogenv=struct;
exogenv.exnames={'Y','sig_epsY','sig_epsM'};
exogenv.exnpt=exnpt;
exogenv.pts_perm=mpts_perm;
exogenv.pts_all=mpts_all;
exogenv.mtrans=mtrans;

          
% ===========================================
% create model object
% ===========================================

% endogenous state variables
endogenv=struct;
endogenv.ennames={'whatM','eI'};

% stv is where we save the solution of the SS, here we call it to
% generate a matrix of initial guesses.

% assign values for initial guess
[solguessvec,Vguessvec,Vnames]=guessfun(stv,expdef.params.LTV,expdef.params.PTI,expdef.params.rental);
            
% save name lists
endogenv.solnames=fieldnames(stv.Sol);
endogenv.solbase=solguessvec;
endogenv.addnames=fieldnames(stv.Add);
endogenv.Vnames=Vnames;
endogenv.condnames = {'expR_MM', 'expR_MY', 'expR_H', 'expR_Bext', 'expR_Bin',...
                      'min_eInext', 'cyieldM', 'cyieldY', 'EDefY', 'EDefM'};                        


% initial definition of grids based on steady state values of state vars
unigrids={1:exogenv.exnpt,expdef.wMpts,expdef.eIpts};
basegrid=grid.TensorGrid(unigrids);
            
if isequal(guess_mode,'no_guess')
    inigrid=basegrid;
    solmatguess=repmat(solguessvec',[basegrid.Npt,1]);
    Vmatguess=repmat(Vguessvec',[basegrid.Npt,1]);
    Forecguess=repmat(stv.Rebate,[basegrid.Npt,exogenv.exnpt]);
    transguess=kron(basegrid.Pointmat(:,2:end), ones(1,exogenv.exnpt));
else
    %inigrid=grid.TensorGrid(basegrid.StateBounds,[exnpt,20,20,15]);
    inigrid=basegrid;
    guessobj=load(guess_path,'mobj');
    solmatguess=guessobj.mobj.evaluatePol(inigrid.Pointmat)';
    Vmatguess=guessobj.mobj.evaluateVal(inigrid.Pointmat)';
    Forecguess=guessobj.mobj.evaluateForec(inigrid.Pointmat)';
    transguess=guessobj.mobj.evaluateTrans(inigrid.Pointmat)';
    
end                  
            
% build approximating function
% create guess for solution variable functions
Pf=grid.LinearInterpFunction(inigrid,solmatguess);
% create guess for next-period functions
Vf=grid.LinearInterpFunction(inigrid,Vmatguess);
% and for state transition function
Tf=grid.LinearInterpFunction(inigrid,transguess);
% forecasting function for dealing with Rebates
Ff=grid.LinearInterpFunction(inigrid,Forecguess);

            
% create new object
% Here we construct mobj, this will be saved of the value functions.

mobj=instfun(expdef.params,endogenv,exogenv,Vf,Pf,Tf,Ff,basegrid,nd,wgt,stv.Add);

disp('-------------------------------------');
disp('Bounds of state space:');
disp(num2str(mobj.Vfct.SSGrid.StateBounds));


% ===========================================
% save initial object
% ===========================================

save(outfname,'stv','mobj');  

end