if ~exist('clear_flag', 'var'), clear_flag = 1; end

if usejava('desktop') && clear_flag
   clear; 
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
if ~exist('resfile','var')
    disp('resfile not defined. Opening default file instead.');
   resfile_list={'res_20200904_lowratefinal'}; 
end


for si=1:length(resfile_list)

resfile = resfile_list{si};    
    
load([respath,resfile,'.mat']);

% Update params
augmentParams=1;
expdef='experdef_20200904.m';
if augmentParams
    run(expdef);
    mobj=mobj.augmentParams(allexpers.bench_ini0.params);
end


% number of periods and burn-in
NT_sim=5000;
NT_ini=500;

% compute Euler equation error?
compEEErr=1;

% Winsorize
winsorize=0;
winsfile='';
cutoff=99.9;

% Force the creation of a sim_res file
force_output = 1;

% output table file
% the folder name will be Results
if ~exist('output_dir', 'var')
    output_dir = './Results/';
end

% create the folder
if ~exist(output_dir, 'dir')
    mkdir(output_dir)
end

% outpath='./Results_c/';
outstats_exog=[output_dir,'statsexog_',resfile,'.xls'];
errstats=[output_dir,'errstats_',resfile,'.xls'];
expdata=0;
outdata=[output_dir,'series_',resfile,'.csv'];
       
%--------------------------------------------------------------------------
% start simulation
%--------------------------------------------------------------------------

% set starting point

% Because the initial shock is the third one?
start_ex=3;
startpt=struct;
startpt.whatM=stv.State.whatM;
startpt.eI=stv.State.eI;
startpt=orderfields(startpt,mobj.En_names);
% Call the function inside of DSGEModel.m that converts the sttucture to cell.
% And the trasnforms the cell to a matrix.
startpt_vec=model.DSGEModel.structToVec(startpt)';
% Put all the state variables in the initialization in the vector.
startpt_vec=[start_ex,startpt_vec];

% aaa = struct2cell(startpt);
% bbb = cell2mat(aaa);

% simulate
% function inside OLGIntermediaryModel.m
[simseries,varnames,errmat,Wshtrans,SDFmat]=mobj.simulate(NT_sim+1,NT_ini,startpt_vec,compEEErr);
simseries_orig=simseries;
varnames_orig=varnames;
statevec = simseries(:,1);

% This function is in OLGIntermdiaryModel.m
[simseries, varnames] = mobj.computeSimulationMoments(simseries,varnames);
nvars = length(varnames);

% Create table object for easier access
% It will implicitely eliminate some repetated variables (mainly some
% prices).
simtable=array2table(simseries);
[~,ia,~]=unique(varnames);
simtable=simtable(:,ia);
simtable.Properties.VariableNames=varnames(ia);
dispnames=varnames(ia);

% make HashMap with mapping of names to indices
indexmap=java.util.HashMap;
for i=1:nvars
    indexmap.put(varnames{i},i);
end

% Check transition function errors
if compEEErr
    idx = sub2ind([NT_sim+1,mobj.Exogenv.exnpt],(1:NT_sim+1)',[statevec(2:end);1]);
    idx=idx(1:end-1);
    WhatMtrans=Wshtrans(:,1:mobj.Exogenv.exnpt);
    WhatM_err=simseries(:,indexmap.get('whatM')) - WhatMtrans(idx);
    errmat = [errmat, [WhatM_err;0]];
    eItrans=Wshtrans(:,mobj.Exogenv.exnpt+1:end);
    eI_err=simseries(:,indexmap.get('eI')) - eItrans(idx);
    errmat = [errmat, [eI_err;0]];
end

% WHAT IS THIS?
if winsorize
    if isempty(winsfile)
        prct=prctile(abs(errmat(2:end,end)),cutoff);
        bad_idx = any(abs(errmat(2:end,end)) > repmat(prct,NT_sim,1) , 2);
    else
       load([respath,winsfile,'.mat'],'bad_idx');
    end
    simseries = simseries(~bad_idx,:);
    errmat = errmat([true;~bad_idx],:);
    NT_sim = size(simseries,1) + 1;
    save([respath,resfile,'.mat'],'bad_idx','-append');
end

%--------------------------------------------------------------------------
% calculate stats
%--------------------------------------------------------------------------

% Recover the values of the state variables
varst=zeros(length(startpt_vec)-1,1);
for i=1:length(startpt_vec)-1
    varst(i)=indexmap.get(mobj.En_names{i});
end
            
% STATE VARIABLES means in stationary distribution
stvstat=mean(simseries(:,varst));


% calculate business cycle stats
statsout_exog=cell(4,1);
statsout_endog=cell(4,1);


mu_Y=exp(mobj.Params.mu_y);
% first for all periods, then separately for low and high sigma_omega states
% THERE ARE '4 EXPERIMENTS' HERE:
% 1 .taking the whole simulation
% 2. taking those simulations where the sigma_y is low and the mu_y is high
% 3. taking those simulations where the sigma_y is low and the mu_y is low
% 4. taking those simulations where the sigma_y is high and the mu_y is high

smpsel_exog={true(NT_sim-1,1), simseries(:,2)==mobj.Params.sig_epsY(1) & simseries(:,1) >=mu_Y, ...
                          simseries(:,2)==mobj.Params.sig_epsY(1)  & simseries(:,1) < mu_Y, ...
                          simseries(:,2)==mobj.Params.sig_epsY(2)  & simseries(:,1) < mu_Y};

gdp_idx = indexmap.get('C'); 

% Looping over the experiment proposed before
for j=1:numel(smpsel_exog)
    % taking only those simulated periods that are in accordance to the specification
    % of the experiment parameters.
    simtmp=simseries(smpsel_exog{j},:);
    
    % For each variable, we will calculate some statisitcs
    statstmp=zeros(nvars,7);
    
    %mean
    statstmp(:,1)=nanmean(simtmp)';
    %std
    statstmp(:,2)=nanstd(simtmp)';
    % contemp and first-order autocorrelations
    autocorrm=corrcoef([simtmp(2:end,:),simtmp(1:end-1,:)]);
    conm=autocorrm(1:nvars,1:nvars);
    lagm=autocorrm(nvars+1:end,1:nvars);
    % corr with shocks
    statstmp(:,3:4)=[conm(:,1),lagm(1,:)'];
    % corr with C
    statstmp(:,5:6)=[conm(:,gdp_idx),lagm(gdp_idx,:)'];
    % vector with fo autocorr
    statstmp(:,7)=diag(lagm);
    statsout_exog{j}=statstmp;
    
end
% The results are saved in the statstmp matrix
% saving them further in statsout_exog

%--------------------------------------------------------------------------
% output
%--------------------------------------------------------------------------

% overview output for eyeball check against analytic st.st. values

% make one big structure with steady-state values
stvbig=model.HelperCollection.combineStructs({stv.Sol,stv.State,stv.Add});

% output table
% make index vector
% displist is the list of the index list
% dispnames is the list linked to the index list before
[displist,dispnames]=model.HelperCollection.makeListFromNames(indexmap,dispnames);
% Notice that ndvars is less thatn nvars because there where some repeated
% variables.
ndvars=length(displist);

disp(' ');
disp('Simulation steady state');

% overview output 
fprintf('Frequency (exog subsamples): ');
% This first loop is just to show the numbers of simulated periods that
% satisify the restiction in each experiment before. For instace, the first
% should have 5,000 periods.
for j=1:numel(smpsel_exog)
    % select vars, cut the repated ones, saved the in a new place called tabout_exog
    tabout_exog{j}=statsout_exog{j}(displist,:);
    fprintf('%f\t',sum(smpsel_exog{j}));
end
fprintf('\n');
disp('-------------');

for s=1:ndvars
    if isfield(stvbig,dispnames{s})
        ststval=stvbig.(dispnames{s});
    else
        ststval=0;
    end
    if numel(dispnames{s}) > 7
        fprintf('%d\t%4s\t\t\t\t%f |',displist(s),dispnames{s},ststval);
    else
        fprintf('%d\t%4s\t\t\t\t\t%f |',displist(s),dispnames{s},ststval);
    end
    % disp('Exog subsamples')

    % use the results saved in tabout to show the mean and std.
    % The results will show 4 colums, each one represents one of the
    % experiments.
    for j=1:numel(smpsel_exog)
        fprintf('\t%f, %f |',tabout_exog{j}(s,1),tabout_exog{j}(s,2));
    end
    fprintf('\n');    
end

% Next show the EE
% The first 16 are from the 16 equations in the model
% the other 2 are for the 2 state variables.
errtab = [];
if compEEErr
    avg_err=mean(abs(errmat))';
    med_err=median(abs(errmat))';
    p75_err=prctile(abs(errmat),75)';
    p95_err=prctile(abs(errmat),95)';
    p99_err=prctile(abs(errmat),99)';
    p995_err=prctile(abs(errmat),99.5)';
    max_err=max(abs(errmat))';
    errtab=table(avg_err,med_err,p75_err,p95_err,p99_err,p995_err,max_err);
    errarr=table2array(errtab);
    disp(' ');
    disp('-----------------------------------------------');
    disp('Average and maximum Euler equation error');
    fprintf('Equ.no.\t\tAvg.\t\tMed.\t\tp75\t\t\tp95\t\t\tp99\t\t\tp99.5\t\tMax.\n');
    for s=1:length(avg_err)
        fprintf('%d\t\t\t%f\t%f\t%f\t%f\t%f\t%f\t%f\n',s,errarr(s,1),errarr(s,2),errarr(s,3), ...
            errarr(s,4),errarr(s,5),errarr(s,6),errarr(s,7));
    end
    
%     % plot EE error for these equations
%     plotEE_pol=[3,4];
%     plotEE_state=[0,0];
%     for i=1:length(plotEE_pol)
%         points=simseries(:,[6,7]);
%         errvals=abs(errmat(1:end-1,plotEE_pol(i)));
%         if plotEE_state(i)>0
%             itmp=(statvec==plotEE_state(i));
%             points=points(itmp,:);
%             errvals=errvals(itmp,:);
%         end
%         model.HelperCollection.scatterPoints2D(points,errvals);
%     end
    
end

% Check grid bounds: this will help us to get better grids.
state_range=4:5;
min_vec=min(simseries(:,state_range));
max_vec=max(simseries(:,state_range));
disp('State bounds:');
disp(mobj.Pfct.SSGrid.StateBounds(:,2:end));
disp('Simulation mins:');
disp(min_vec);
disp('Simulation max:');
disp(max_vec);


% write to file
values=struct2cell(mobj.Params);
paramout=cell2table(values,'RowNames',fieldnames(mobj.Params));
colnames={'mean','std','corrG','corrG_1','corrC','corrC_1','AC'};

% Save the results of the 4 experiments in a spreedsheet in excel.
% They will include all the stattisicts in colnames
for j=1:numel(smpsel_exog)
    tableout_exog=array2table(tabout_exog{j},'RowNames',dispnames,'VariableNames',colnames);
    writetable(tableout_exog,outstats_exog,'WriteRowNames',1,'FileType','spreadsheet','Sheet',j);
end

% Add the value of the parameters of each experiment
writetable(paramout,outstats_exog,'WriteRowNames',1,'FileType','spreadsheet','Sheet','params');

if compEEErr
% Save the EE errors.
    writetable(errtab,errstats,'FileType','spreadsheet');
end

% Saving everything in sim_res_test2.mat
if force_output
    params=mobj.Params;
    disp(['Saving simulation data to .mat file: ',['sim_',resfile,'.mat']]);
    save(['sim_',resfile,'.mat'],'simseries','displist','dispnames','errmat','tabout_exog','outstats_exog', ...
                'errstats','errtab','indexmap','NT_ini','NT_sim','smpsel_exog','statevec','statsout_exog','varnames','params');
end    


if expdata
    disp(' ');
    disp('Exporting simseries...');
    model.HelperCollection.tableExport(outdata,varnames,simseries);
end

% save model file with stationary state values
save([respath,resfile,'.mat'],'stvstat','-append');


end

