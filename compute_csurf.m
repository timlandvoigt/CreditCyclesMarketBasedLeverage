if usejava('desktop')
   clear; 
end

close all;
respath='./';
% outpath='./Results2/';
outpath='./Results/';
if ~exist('resfile','var')
    resfile='res_20200904_bench';
    %resfile='res_20200904_lowratefinal_belief';
end
printfile=[outpath,resfile,'_csurf'];
outfile=['CS_',resfile,'.mat'];

load([respath,resfile,'.mat']);
load([respath,'sim_',resfile,'.mat']);

open_parpool;

% get average housing share, LTV, and LTW
hY=simseries(:,indexmap.get('hY'));
hM=simseries(:,indexmap.get('hM'));
mY=simseries(:,indexmap.get('mY'));
mM=simseries(:,indexmap.get('mM'));
levYH=simseries(:,indexmap.get('lev_mY'));
levMH=simseries(:,indexmap.get('lev_mM'));
wYYdef=simseries(:,indexmap.get('wYYdef'));
pY=simseries(:,indexmap.get('pY'));
pM=simseries(:,indexmap.get('pM'));
wMdef=simseries(:,indexmap.get('wMdef'));
P=simseries(:,indexmap.get('P'));
% levYW=mean(mY)/mean(wYYdef-pY);
% levMW=mean(mM)/mean(wMdef-pM);
levYW=mean(mY)/mean(wYYdef);
levMW=mean(mM)/mean(wMdef);
whatM=simseries(:,indexmap.get('whatM'));
eI=simseries(:,indexmap.get('eI'));

hYgrid=linspace(min(hY),max(hY),6);
hMgrid=linspace(min(hM),max(hM),6);
levYHgrid=linspace(mean(levYH)-0.2,mean(levYH)+0.2,10);
levMHgrid=linspace(mean(levMH)-0.1,mean(levMH)+0.1,6);
levYWgrid=linspace(levYW-0.05,levYW+0.05,8);
levMWgrid=linspace(levMW-0.05,levMW+0.05,8);

% credit surface variables
Ygrid=grid.TensorGrid({hYgrid,levYHgrid,levYWgrid});
Yvar=Ygrid.Pointmat;
Mgrid=grid.TensorGrid({hMgrid,levMHgrid,levMWgrid});
Mvar=Mgrid.Pointmat;

% histogram counts of stationary distribution
ungr=mobj.Pfct.SSGrid.Unigrids;
stategrid=[statevec(2:end),simseries(:,4:5)];
exst_edges=ungr{1};
whatM_edges=ungr{2};
whatM_edges=[whatM_edges(1),(whatM_edges(1:end-1)+whatM_edges(2:end))/2,whatM_edges(end)];
eI_edges=ungr{3};
eI_edges=[eI_edges(1),(eI_edges(1:end-1)+eI_edges(2:end))/2,eI_edges(end)];
dist=histcn(stategrid,exst_edges,whatM_edges,eI_edges);
prdist=dist/NT_sim;
indmat=grid.StateSpaceGrid.makeCombinations_rev(mobj.Pfct.SSGrid.Dimvec);
npt=mobj.Pfct.SSGrid.Npt;
prvec=zeros(npt,1);
for p=1:npt
    prvec(p)=prdist(indmat(p,1),indmat(p,2),indmat(p,3));
end

% compute surfaces
[QYmat,QMmat]=mobj.calcCSurf(Yvar,Mvar);

% average surfaces
levYWgridplot=levYWgrid*mean(wYYdef)/mean(wYYdef-pY);
levMWgridplot=levMWgrid*mean(wMdef)/mean(wMdef-pM);
levYWplot=mean(mY)/mean(wYYdef-pY);
levMWplot=mean(mM)/mean(wMdef-pM);
Ygridplot=grid.TensorGrid({hYgrid,levYHgrid,levYWgridplot});
Mgridplot=grid.TensorGrid({hMgrid,levMHgrid,levMWgridplot});
figs=plot_csurf(Ygridplot,Mgridplot,QYmat,QMmat,prvec,{[mean(hY),mean(levYH),levYWplot],[mean(hM),mean(levMH),levMWplot]},[],[],[],3);
zlab={'Y','M'};
for f=1:2
%     ax=get(figs{f},'CurrentAxes');
%     axes(ax);
%     hold(ax,'on');
%     %plot3(ax,xpt{f},ypt{f},zpt{f},'LineWidth',3,'Color','k');
%     arrow3(levpts0(f,:),levpts1(f,:),'k-2');
    print(figs{f},'-painters','-depsc',[printfile,zlab{f},'_avg.eps']);    
end

stategrid=mobj.Pfct.SSGrid.Pointmat;
save(outfile,'Ygrid','Mgrid','QYmat','QMmat','stategrid');

% % conditional on bank net worth
% point1=[2,mean(whatM),.025];
% dist=sum((stategrid-repmat(point1,size(stategrid,1),1)).^2,2);
% [~,distidx]=sort(dist);
% prvec1=zeros(npt,1);
% prvec1(distidx(1))=1;
% point2=[3,mean(whatM),.1];
% dist=sum((stategrid-repmat(point2,size(stategrid,1),1)).^2,2);
% [~,distidx]=sort(dist);
% prvec2=zeros(npt,1);
% prvec2(distidx(1))=1;
% prvec_ei=[prvec1,prvec2];
% 
% figs=plot_csurf(Ygrid,Mgrid,QYmat,QMmat,prvec_ei,{[mean(hY),levYW],[mean(hM),levMW]},{'r-','b-'},[],2);
% 






