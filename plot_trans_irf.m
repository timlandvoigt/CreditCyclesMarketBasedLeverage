if ~exist('clear_flag', 'var'), clear_flag = 1; end

if usejava('desktop') && clear_flag
   clear; 
end

close all;
respath='./';
% outpath='./Results2/';
outpath='./Results/';
if ~exist('resfile','var')
    resfile='res_20200904_bench';
end
outfile=['GR_',resfile];
grayscale=0;
reportLevels=1;
batchMode=0;
relative_irf=1;

load([respath,resfile,'.mat']);
load([respath,'sim_',resfile,'.mat']);
load([respath,'GR_',resfile,'.mat']);

close all;

% Make Graphs
outpath=[outpath,outfile,'_'];
if relative_irf, outpath = [outpath, 'diff_']; end
tvec=0:NT_sim-1;

% file names for graphs (set to empty for no printing)
printfiles={[outpath,'IRF1'],[outpath,'IRF2'],[outpath,'IRF3'],[outpath,'IRF4']};        
% which variables
brsel1=[indexmap.get('Y'),indexmap.get('sig_epsY'),indexmap.get('P'),indexmap.get('Lrate'),...
        indexmap.get('totlev_m'), indexmap.get('totdep'), indexmap.get('Mspr'),indexmap.get('eI')];
brsel2=[indexmap.get('Defrate'),indexmap.get('MsprY'),indexmap.get('MsprM'),...
        indexmap.get('eI'),indexmap.get('totdep'),indexmap.get('R')];
brsel3=[indexmap.get('Lrate'),indexmap.get('totmdebt'),indexmap.get('xI'), ...
        indexmap.get('cY'),indexmap.get('cM'),indexmap.get('cO')];
brsel4=[indexmap.get('DefY'), indexmap.get('DefM'),indexmap.get('dY'), ...
    indexmap.get('dM'), indexmap.get('DWL'),indexmap.get('I')];


% How many shocks to plot (one less for relative)
if relative_irf
    N_shock_plot = N_shock - 1;
else
    N_shock_plot = N_shock;
end

brsel_all=[brsel1,brsel2,brsel3,brsel4];
nvar=length(brsel_all);   
brseries_gr=zeros(N_shock_plot, NT_sim, nvar);

for s=1:N_shock_plot
    if relative_irf
        if reportLevels==1
            brseries_gr(s,:,:) = simseries_diff_mean{s+1}(:,brsel_all);
        else
            brseries_gr(s,:,:) = 100*(simseries_diff_mean{s+1}(:,brsel_all) ./ repmat(simseries_diff_mean{s+1}(1,brsel_all),NT_sim,1)-1);
        end
    else
        if reportLevels==1
            brseries_gr(s,:,:) = simseries_mean{s}(:,brsel_all);
        else
            brseries_gr(s,:,:) = 100*(simseries_mean{s}(:,brsel_all) ./ repmat(simseries_mean{s}(1,brsel_all),NT_sim,1)-1);
        end
    end
end

colors={'b-o','r-o'};
if N_shock_plot==3
   colors = ['k-o',colors]; 
end

titles1={'Income','Housing Risk','House price','Loss rate','Leverage', 'Deposits','Mortg. Spread','Bank equity'}; 
titles2={'Default rate','Mortg. spread Y','Mortg. spread M','Bank equity','Deposits','Riskfree rate'};
titles3 = {'Loss rate','Mort Debt', 'Net dividend','Consumption Y','Consumption M','Consumption O'};
titles4 = {'Def rate Y', 'Def rate M', 'Deposits Y', ...
         'Deposits M', 'DWL', 'Equ issuance'};
     
if usejava('desktop')
    makeImpRes(brseries_gr(:,:,1:8),tvec,titles1,colors,[2,4],[],printfiles{1});   
end


if batchMode
    close all;
    clear;
end
