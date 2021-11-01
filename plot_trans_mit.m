respath='./';
outpath='./Results/';

% Which experiment to run:
resfiles = {'20200904'};
res_suffix='_MITloss';

plot_exper=[1];


labels=[];

colors_base={'k-o','b-o','r-o'};
colors={'k-o','b-o','r-o'};


printplots=true;
printfile = [outpath,'transitions_',resfiles{end},res_suffix];


grayscale=0;
reportLevels=1;
batchMode=0;



close all;

%% Make Plots

% titles{1} = {'Riskfree rate', 'House price', 'Mortg. debt','Leverage',...
%            'Default rate','Bank equ. ratio','Mortg. spread','Consumption'};
titles{1} = {'Loss rate', 'LTV', 'Bank equ.','Mortg. risk premium'};
       
titles{2} = {'Bank equ. ratio', 'Loss rate' ,'LTV', 'Deposits',...
           'Consumption Y','Consumption M', 'Risk premium M', 'Risk premium Y'};
titles{3} = {'Bank equ. ratio', 'Mortg. debt','Leverage','Deposits',...
           'Default rate','Mortg. spread','Exp. Excess Ret Y','Exp. Excess Ret M'};
    
titles{4} = {'Mort debt Y', 'Housing Y','Leverage Y',...
           'Consumption Y', 'Mortg. spread Y', 'Exp. Excess Ret Y'};

titles{5} = {'Mort debt M', 'Housing M','Leverage M',...
           'Consumption M', 'Mortg. spread M', 'Exp. Excess Ret M'};

titles{6} = {'Riskfree rate', 'Price-rent (FHFA)','Mort.debt/income (SCF)','Spread Y-M (Pr-Subpr)'};       

NT_sim = 26;
N_exper = length(plot_exper);


npl=length(titles);
plotformats=zeros(npl,1);
simseries_plots=cell(npl,1);
for p=1:npl
    plotformats(p)=length(titles{p});
    simseries_plots{p}=zeros(N_exper,NT_sim,plotformats(p));
end



for i=1:N_exper
	%tmp=load(['PT_',resfiles{i},res_suffix,'.mat']);
    tmp=load(['ST_',resfiles{i},res_suffix,'.mat']);
	tmp_simseries = tmp.simseries_mean{1};
    indexmap=tmp.indexmap;
	
    brsel{1} = [indexmap.get('Lrate'), indexmap.get('totlev_m'), indexmap.get('eI'), indexmap.get('expER_YMdiff')];
    
    brsel{2} = [indexmap.get('eIrat_m'), indexmap.get('Lrate'), indexmap.get('totlev_m'),indexmap.get('DI'), ...
        indexmap.get('cY'), indexmap.get('cM'), indexmap.get('expER_MMnet'), indexmap.get('expER_MYnet')];
    
    brsel{3} = [indexmap.get('eIrat_m'), indexmap.get('totmdebt'), indexmap.get('totlev_m'), indexmap.get('DI'), ...
        indexmap.get('Defrate'), indexmap.get('Mspr'), indexmap.get('expER_MY'), indexmap.get('expER_MM')];
    
    
    brsel{4} = [indexmap.get('mY'), indexmap.get('hY'), indexmap.get('lev_mY'), ...
        indexmap.get('cY'), indexmap.get('MsprY'), indexmap.get('expER_MY')];
    
    brsel{5} = [indexmap.get('mM'), indexmap.get('hM'), indexmap.get('lev_mM'), ...
        indexmap.get('cM'), indexmap.get('MsprM'), indexmap.get('expER_MM')];
    
     brsel{6} = [indexmap.get('rf'), indexmap.get('PRavg'), indexmap.get('totdti'), indexmap.get('MsprYM')];
    
    
    for p=1:length(brsel)
        thisplot=simseries_plots{p};
        thisplot(i,:,:)=tmp_simseries(:,brsel{p});
        simseries_plots{p}=thisplot;
    end

end

tvec=0:15;
       

% Plots for all the variables
for p=[2]
    thisprint=[];
    if printplots
        thisprint=[printfile,'_',num2str(p)];
    end
    plotdata=simseries_plots{p}(:,tvec+1,:);
    nplots=size(plotdata,3);
    makeImpRes(simseries_plots{p}(:,tvec+1,:),tvec,titles{p},colors,[2,4],labels,thisprint);
end

