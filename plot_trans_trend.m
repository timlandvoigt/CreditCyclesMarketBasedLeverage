respath='./';
outpath='./Results/';

% Which experiment to run:
%resfiles = {'20200904_eta'};
resfiles = {'20200904_eta','20200904_secur','20200904_belief'};
%resfiles = {'20200904_secur','20200904_secur_LTV','20200904_secur_DTI'};
res_suffix='_simtrend';



labels={'asset demand','+ dereg.', '+ credit risk'};
%labels={'no constr','LTV', 'DTI'};
%labels={'base'};

colors_base={'k-o','b-o','r-o'};
colors={'k-o','b-o','r-o'};

plot_exper=[1,2,3];
%plot_exper=[1];

printplots=true;
printfile = [outpath,'transitions_',resfiles{end},res_suffix];


grayscale=0;
reportLevels=1;
batchMode=0;

dattab=readtable('graph_data_SCF.csv');
data4=[dattab.rrt/100,dattab.PR,dattab.SCFdti,dattab.mspreaddiff100];

close all;

%% Make Plots

% titles{1} = {'Riskfree rate', 'House price', 'Mortg. debt','Leverage',...
%            'Default rate','Bank equ. ratio','Mortg. spread','Consumption'};
titles{1} = {'House price', 'LTV','Cons Y','Bank equ.','Loss rate','Mortg. risk premium'};
       
titles{2} = {'Bank equ. ratio', 'Default rate' ,'LTV', 'Mortg. spread',...
           'Consumption Y','Consumption M', 'Exp. Excess Ret Y', 'Exp. Excess Ret M'};
titles{3} = {'Bank equ. ratio', 'Mortg. debt','Leverage','Deposits',...
           'Default rate','Mortg. spread','Exp. Excess Ret Y','Exp. Excess Ret M'};
    
titles{4} = {'Mort debt Y', 'Housing Y','Leverage Y',...
           'Consumption Y', 'Mortg. spread Y', 'Exp. Excess Ret Y'};

titles{5} = {'Mort debt M', 'Housing M','Leverage M',...
           'Consumption M', 'Mortg. spread M', 'Exp. Excess Ret M'};

titles{6} = {'Riskfree rate', 'Price-rent (FHFA)','Mort.debt/income (SCF)','Spread Y-M (Subpr-Pr)'};       

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
	
    brsel{1} = [indexmap.get('P'), indexmap.get('totlev_m'), indexmap.get('cY'),...
                indexmap.get('eI'),indexmap.get('Lrate'), indexmap.get('expER_YMdiff')];
    
    brsel{2} = [indexmap.get('lev_mY'), indexmap.get('lev_mM'), indexmap.get('totlev_m'),indexmap.get('Mspr'), ...
        indexmap.get('cY'), indexmap.get('cM'), indexmap.get('expER_MY'), indexmap.get('expER_MM')];
    
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
for p=[1]
    thisprint=[];
    if printplots
        thisprint=[printfile,'_',num2str(p)];
    end
    plotdata=simseries_plots{p}(:,tvec+1,:);
    nplots=size(plotdata,3);
    makeImpRes(simseries_plots{p}(:,tvec+1,:),tvec,titles{p},colors,[2,3],labels,thisprint);
end
% makeImpRes(simseries_mean(:,tvec+1,brsel5),tvec,titles5,colors,[2,3],labels,[printfile,'_5']);
% makeImpRes(simseries_mean(:,tvec+1,brsel6),tvec,titles6,colors,[2,3],labels,[printfile,'_6']);

% Plots for the variables with data
printdat=[];
if printplots
    printdat=[printfile,'_data'];
end
series_model=simseries_plots{6}(:,tvec+1,:);
PRrat = series_model(:,:,2)./repmat(squeeze(series_model(:,1,2)),1,tvec(end)+1)*100;
series_model(:,:,2)=PRrat;
series_dat=[reshape(data4(tvec+1,:),[1,tvec(end)+1,4]); series_model];
makeImpRes(series_dat,tvec,titles{6},{'k--',colors_base{:}},[1,4],[{'data'},labels],printdat,[0.182257626463907,0.701389109065472,0.120929118773946,0.16283422459893]);

