if usejava('desktop')
   clear; 
else
    ver;
end
close all;

if ~exist('no_par_processes','var')
    no_par_processes=16;
    open_parpool;
end


%--------------------------------------------------------------------------
% simulation setup
%--------------------------------------------------------------------------

% file with model
respath='./';
outpath='./Results/';

% Which experiment to run:

% 1. Baseline
experdef='20200904';
experinit = ['res_',experdef,'_bench'];
experfinal_list = {['res_',experdef,'_lowratefinal']};
translabel='_eta';
translabelfin='';
res_names = {['res_',experdef,'_trans01',translabel],...
             ['res_',experdef,'_trans02',translabel],...
             ['res_',experdef,'_trans03',translabel],...
             ['res_',experdef,'_trans04',translabel],...
             ['res_',experdef,'_trans05',translabel],...
             ['res_',experdef,'_trans06',translabel],...
             ['res_',experdef,'_trans07',translabel],...
             ['res_',experdef,'_trans08',translabel],...
             ['res_',experdef,'_lowratefinal',translabelfin]};


out_suffix='_simtrend';



% output table file
outfile=['ST_',experdef,translabel,out_suffix];


N_exper=numel(experfinal_list);

% Initial Economy Config
varlist={'simseries','statevec','indexmap','varnames'};
load(['sim_',experinit],varlist{:});
startmobjstr=load(experinit,'mobj');
startmobj=startmobjstr.mobj;

% compute Euler equation error?
compEEErr=1;


start_shock_sequence=[3,3,3,3,3,3,3,2];
start_ini=3;


shockcell.exst = 3;
shockcell.params = {};
shockcell.pts = {};


if ~isempty(start_shock_sequence)
    N_shock=length(start_shock_sequence);
else
    N_shock=0;
end


% number of periods and burn-in
N_runs=1000;
NT_sim=25;
NT_ini=0;


%% Set starting point
gv = @(x)indexmap.get(x);
% orig_statevec=statevec(2:end);

% N_vars=length(startvals);
% simulate for one period
NT_sim_ini=1;
%warning('NT_sim_ini = 0');
Nvarsim=startmobj.NSTEX + startmobj.NSTEN + startmobj.NSOL + startmobj.NV + ...
			startmobj.NADD + startmobj.NCOND + 1;
		
% Find states such that statevec = start_ini
% Note, statevec(2) corresponds to simseries(1)
idx_ini = find(statevec(2:end)==start_ini);

% Define initial values
startvals=mean(simseries(idx_ini,:));

% Simulate forward to get actual, not averaged, values for variables other
% than state variables
if NT_sim_ini>0
	firstrow=[start_ini, startvals(1:Nvarsim)];
%	startpt_vec=startvals([1,1+[gv('KB'),gv('LB'),gv('WI'),gv('BG')]]);
	startpt_vec=[start_ini,startvals([gv('whatM'),gv('eI')])];
	transprob=cumsum(startmobj.Exogenv.mtrans(start_ini,:));
	shock_prob=transprob(start_ini);
	if start_ini>1
		shock_prob_minus=transprob(start_ini-1);
	else
		shock_prob_minus=0;
	end
	rvar_next=(shock_prob+shock_prob_minus)/2;
	[simseries_ini,varnames_ini,~,~,~]=startmobj.simulate(NT_sim_ini,0,startpt_vec,compEEErr,rvar_next*ones(NT_sim_ini,1));

	[simseries_ini, varnames_ini] = startmobj.computeSimulationMoments(simseries_ini,varnames_ini,firstrow);
	startvals=simseries_ini(end,:);
end
N_vars=length(startvals);


% dY, dM, hY, mY, mM, hM
% 

% initial point in terms of quantities
% Because we don't have start-of-period quantities, use the previous period's
% end-of-period variables 
quantenstates = startvals([gv('dY'), gv('hY'), gv('mY'), gv('dM'), gv('hM'), gv('mM')]);
quantpoint = [start_ini, quantenstates];

% Prep initial guess for MIT state
guess_statevec=statevec(2:end);
guess_enstates = simseries(:, [ gv('whatM'), gv('eI') ] );

varnames_store = varnames;

simseries_median = cell(N_exper,1);
simseries_mean = cell(N_exper,1);
simseries_std = cell(N_exper,1);



for nexp=1:N_exper
    
    % final state experiment
    resfinal=experfinal_list{nexp};
    mobjfinal=getfield(load([respath,resfinal,'.mat'],'mobj'),'mobj');  
    if ~isfield(mobjfinal.Params,'rental')
        newparams=mobjfinal.Params;
        newparams.rental=false;
        mobjfinal=mobjfinal.augmentParams(newparams);
    end
    
    % trend path
    if ~isempty(res_names)
        N_steps = length(res_names);
        mobjlist=cell(N_steps,1);
        for nst=1:N_steps
            thismobj=getfield(load(res_names{nst},'mobj'),'mobj');
            if ~isfield(thismobj.Params,'rental')
                newparams=thismobj.Params;
                newparams.rental=false;
                thismobj=thismobj.augmentParams(newparams);
            end
            mobjlist{nst}=thismobj;
        end
    else % or just one-time shock?
        N_steps=1;
        mobjlist{1}=mobjfinal;
    end

    
    disp(['Experiment ',num2str(nexp),' of ',num2str(N_exper)]);
    tens_simseries = zeros(NT_sim+1,N_vars,N_runs);
    
    start_shock = shockcell.exst;
	idx_switch = 1 + find( guess_statevec(2:end)==start_shock & guess_statevec(1:end-1)==start_ini);
    guessenstates = mean(guess_enstates(idx_switch,:));

    simmit=zeros(N_steps,Nvarsim+1);
    % note: indices from startmobj
    quantindex=[gv('dY'), gv('hY'), gv('mY'), gv('dM'), gv('hM'), gv('mM')]+1;
    if mobjfinal.Params.rental
        quantexpindex=[1,3,4,6];
    else
        quantexpindex=[1:4,6];
    end

    % vectorize all params for MIT shocks,
    % since by the nature of MIT shocks parameters change over time
    % need to make all params (NT_sim+1) x 1 vectors
    % first period is initial period before shock hits
    paramNames = fieldnames(startmobj.Params);
    mitparams_sim = startmobj.Params;
    for p = paramNames'
        p = p{:};
        origval = startmobj.Params.(p);
        if ~isa( origval, 'function_handle')
            origval = origval(:)';
            mitparams_sim.(p) = repmat( origval, NT_sim+1, 1 );
        end
    end
    
    % run consecutive MIT shocks
    for nst=1:N_steps
        disp(['MIT shock period: ',num2str(nst)]);
        % compute equilibrium and transitions for period when MIT shock hits
        [transmat,simvec,mitparams] = computeFirstPeriod(mobjlist{nst},quantpoint,guessenstates,shockcell,0);
        simmit(nst,:)=simvec;
        if nst < N_steps
            % get initial point for next shock ready
            quantpoint=simvec(quantindex);
            quantpoint(quantexpindex)=exp(quantpoint(quantexpindex));
            quantpoint=[start_shock_sequence(nst),quantpoint];
            guessenstates=simvec([ gv('whatM'), gv('eI') ]+1);
            shockcell.exst = start_shock_sequence(nst);
        end
        % update changed parameters
        paramNames = fieldnames(mitparams);
        for p = paramNames'
            p = p{:};
            mitval = mitparams.(p);
            if ~isa( mitval, 'function_handle')
                mitval = mitval(:)';
                newparams=mitparams_sim.(p);
                newparams(nst+1,:)=mitval; % skip first entry
                mitparams_sim.(p) = newparams;
            end
        end
    end
    
    % compute entry of random number matrix that sets states
    % deterministically to start_shock_sequence
    if N_shock>N_steps-1
        start_shock_sequence=start_shock_sequence(N_steps:end);
        N_shockfinal=length(start_shock_sequence);
        rvar_sequence=zeros(N_shockfinal,1);
        this_mtrans=mobjfinal.Exogenv.mtrans;
        for s=1:N_shockfinal
            thisstate=start_shock_sequence(s);
            transprob=cumsum(this_mtrans(thisstate,:));
            shock_prob=transprob(start_shock_sequence(s));
            if start_shock_sequence(s)>1
                shock_prob_minus=transprob(start_shock_sequence(s)-1);
            else
                shock_prob_minus=0;
            end
            rvar_sequence(s)=(shock_prob+shock_prob_minus)/2;
        end
    else
        N_shockfinal=0;
    end

    % Create shock matrix
    NT_final=NT_sim-N_steps;
    rng(1);
    shmatfull = lhsdesign(N_runs,NT_final+1);
    
    fprintf([repmat('.',1,100) '\n\n']);

    parfor n=1:N_runs       
	%for n=1:N_runs
        %--------------------------------------------------------------------------
        % start simulation
        %--------------------------------------------------------------------------

        shmat = shmatfull(n,:)';
        if N_shockfinal>0
            shmat(1:N_shockfinal)=rvar_sequence;
        end
        simseries_all=zeros(NT_sim,Nvarsim+1);
        simseries_all(1:N_steps,:)=simmit;
        
        exnext=find(transprob-shmat(1)>0,1,'first');
        startpt_vec = [exnext,transmat(exnext,:)];
                             
        % remaining periods with final experiment
        [simseries,varnames]=mobjfinal.simulate(NT_final,0,startpt_vec,compEEErr,shmat(2:end));
        simseries_all(N_steps+1:end,:)=simseries;
        
        simseries_all = mobjfinal.computeSimulationMoments(simseries_all,varnames,firstrow,mitparams_sim);
        
        tens_simseries(:,:,n) = [startvals; simseries_all];
        %aaa = tens_simseries(:,:,1);
        if mod(n,N_runs/100)==0
            %disp([num2str(n),'/',num2str(N_runs),': ',num2str(round(1000*n/N_runs)/10),'% complete']);
            fprintf('\b|\n');
        end

    end
    
    fprintf('\n');
    varnames = varnames_store;
    nvars = length(varnames);

    % make HashMap with mapping of names to indices
    indexmap=java.util.HashMap;
    for i=1:nvars
        indexmap.put(varnames{i},i);
    end
%     varst=zeros(length(startpt_vec)-1,1);
%     for i=1:length(startpt_vec)-1
%         varst(i)=indexmap.get(mobj.En_names{i});
%     end

    %save(outfile,'tens_simseries','indexmap');

    simseries_median{nexp} = median(tens_simseries,3);
    simseries_mean{nexp} = mean(tens_simseries,3);
    simseries_std{nexp} = std(tens_simseries,[],3);

end

save(outfile,'simseries_mean','simseries_median','simseries_std','indexmap','NT_sim','N_shock','start_ini','start_shock_sequence');




function [transmat,simnext,params] = computeFirstPeriod(mobj,quantpoint,guessenstates,shockstruct,n)
	exst=shockstruct.exst;
	params = mobj.Params;
	paramDeviations = shockstruct.params;
	for ii=1:size(paramDeviations,1)
		params.(paramDeviations{ii,1}) = paramDeviations{ii,2};
	end
	
	exogenv = mobj.Exogenv;
	exogDeviations = shockstruct.pts;
	for ii=1:size(exogDeviations,1)
		colIdx = find( strcmp( exogDeviations{ii,1}, exogenv.exnames ) );
		exogenv.pts_perm(exst,colIdx) = exogDeviations{ii,2};
		exogenv.pts_all(exst,colIdx) = exogDeviations{ii,2};
	end

% 	% adjust initial guesses
% 	adjust = cell( mobj.NSOL, 1);
% 	adjust(:) = { @(x) x };
% 	% adjust initial guess for consumption
% 	adjust([6,11]) = { @(x) x - log(1 + params.demandShock * exp(x) ) };
	
	exnpt = mobj.Exogenv.exnpt;
	guesspoint = [exst, guessenstates];
	quantpoint(1) = exst;
	solguess=mobj.evaluatePol(guesspoint);
    Rebate_list=mobj.evaluateForec(guesspoint);
    rebguess = Rebate_list(exst);
% 	for ii=1:mobj.NSOL
% 		solguess(ii) = adjust{ii}(solguess(ii));
% 	end
	transguess=mobj.evaluateTrans(guesspoint);
	stateguess=guessenstates';
	guess=[solguess;transguess;stateguess;rebguess];
	objfun = @(x)computeMITShockState(quantpoint,x,mobj,params,exogenv);

	options=optimset('Display','off','TolX',1e-15,'TolFun',1e-12,...
			'MaxIter',100,'MaxFunEvals',100^2,'FinDiffType','central');
	[sol,fx,exfl] = fsolve(objfun,guess,options);
	if exfl<1
%		error('No Eqm Found in Run %d',n);
		disp('No Eqm Found');
	end
	[~,simnext, transmat]=objfun(sol);
end
