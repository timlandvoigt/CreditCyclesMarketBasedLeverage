function [res,simnext,transmat]  = computeMITShockState( quantpoint , x, mobj, params, exogenv)
%[gv('dY'), gv('hY'), gv('mY'), gv('dM'), gv('hM'), gv('mM')]
exst = quantpoint(1);
dY = quantpoint(2);
hY = quantpoint(3);
mY = quantpoint(4);
dM = quantpoint(5);
hM = quantpoint(6);
mM = quantpoint(7);

NSOL=mobj.NSOL;
exnpt=exogenv.exnpt;

solvec = x(1:NSOL);
whatMtrans = x(NSOL+(1:exnpt));
eItrans = x(NSOL+exnpt+(1:exnpt));
whatM = x(NSOL+2*exnpt+1);
eI = x(NSOL+2*exnpt+2);
Reb_today = x(NSOL+2*exnpt+3);

point = [exst,whatM,eI];

%trans=[(1:exnpt)',KBtrans*ones(exnpt,1),LBtrans,WItrans,BGtrans*ones(exnpt,1)];
trans=[whatMtrans;eItrans];
%trans=reshape(trans,exnpt,4);

[nextst,outstr]=mobj.calcStateTransition(point,solvec,0,trans,params,exogenv);


G = params.mu_G;
delta_eta = params.delta_eta;
deltaH = params.deltaH;
lambdaY = params.lambdaY;
lambdaM = params.lambdaM;
piY = params.piY;
nu = params.nu;
xi = params.xi;
xiDWL = params.xiDWL;
if isfield(params,'MITshockloss')
    MITshockloss = params.MITshockloss;
else
    MITshockloss =0;
end

addvars=outstr.addvars;
Ystar_today=addvars.Y_star_today;
Rebate_today=addvars.Rebate_today;
xI=addvars.xI;
sigepsY=outstr.exstvec(2);
sigepsM=outstr.exstvec(3);
muepsY=-.5*sigepsY.^2;
muepsM=-.5*sigepsM.^2;
Ystar_today = Ystar_today - Rebate_today + Reb_today;

P=exp(solvec(2));
q=exp(solvec(3));
pY=exp(solvec(4));
pM=exp(solvec(5));

% young that stay young
wYYdef = dY./G + nu*Ystar_today + pY;
epsbarYY = (mY/G - lambdaY*wYYdef)./((1-deltaH)*P*hY);
[FepsYY,~,~,FepsminusYY]= OLGIntermediaryModel.integEps(epsbarYY,muepsY,sigepsY,[],[],[],[],[],1);

% wealth of the young that turn middle-aged
ystar_today=(1-nu)*Ystar_today;
wYMdef = dY./G + (delta_eta*ystar_today)/piY;
epsbarYM = (mY/G - lambdaY*wYMdef)./((1-deltaH)*P*hY);
[FepsYM,~,FepsplusYM,FepsminusYM]= OLGIntermediaryModel.integEps(epsbarYM,muepsY,sigepsY,[],[],[],[],[],1);
WYMdef = piY*wYMdef;
WYM = WYMdef.*(1-lambdaY.*FepsYM)...
    + piY*((1-deltaH)*FepsplusYM.*P*hY - (1-FepsYM)*mY/G);
% wealth of middle-aged
wMdef = dM./G + q+xI + ((1-delta_eta)*ystar_today+pM);
epsbarM = (mM/G - lambdaM*wMdef)./((1-deltaH)*P*hM);
[FepsM,~,FepsplusM,FepsminusM]= OLGIntermediaryModel.integEps(epsbarM,muepsM,sigepsM,[],[],[],[],[],1);
WM = WYM + (1-deltaH)*FepsplusM.*P*hM - (1-FepsM)*mM/G + (1-lambdaM*FepsM).*wMdef;
whatM_check = WM - (q+xI) - pM;

% intermediary wealth
DefY = (1-piY)*FepsYY + piY*FepsYM;
FepsminusY = (1-piY)*FepsminusYY + piY*FepsminusYM;
PayoffY = (1-DefY).*mY./G + (1-xi)*(1-deltaH)*FepsminusY.*P*hY;
PayoffM = (1-FepsM).*mM./G + (1-xi)*(1-deltaH)*FepsminusM.*P*hM;
eI_check= PayoffY + PayoffM - (dM+dY)./G - MITshockloss;

% Rebates
RebateY = (1-piY)*lambdaY*FepsYY.*wYYdef + lambdaY*FepsYM.*WYMdef + (xi-xiDWL)*(1-deltaH)*FepsminusY.*P*hY;
RebateM = (xi-xiDWL)*(1-deltaH)*FepsminusM.*P*hM + lambdaM*FepsM.*wMdef;

           
[fx,~,V]=mobj.calcEquations(point(1),nextst,solvec,outstr,1,0);
whatMtrans_check = V{2}(1:exnpt)';
eItrans_check = V{2}(exnpt+(1:exnpt))';

res = zeros(size(x));
res(1:NSOL) = fx;
res(NSOL+(1:exnpt)) = whatMtrans - whatMtrans_check;
res(NSOL+exnpt+(1:exnpt)) = eItrans - eItrans_check;
res(NSOL+2*exnpt+1) = whatM - whatM_check;
res(NSOL+2*exnpt+2) = eI - eI_check;
res(NSOL+2*exnpt+3) = Reb_today - (RebateY+RebateM);

if nargout>1
%	outstr.addvars.sig2B_zeta_xi_next = outstr.addvars.sig2B_zeta_xi_next(5,2);
	addvec=model.DSGEModel.structToVec(outstr.addvars)';
	valvec=V{1};
    [~,~,V]=mobj.calcEquations(point(1),nextst,solvec,outstr,2,0);
    condvec=model.DSGEModel.structToVec(V{1})';
% 	valvec( strcmp(mobj.V_names,'qB') ) = exp( solvec(strcmp(mobj.Sol_names,'qB')) );
% 	valvec( strcmp(mobj.V_names,'X') ) = solvec(strcmp(mobj.Sol_names,'X'));
% 	valvec( strcmp(mobj.V_names,'lamB') ) = solvec(strcmp(mobj.Sol_names,'lamB'));
% 	
% 	valvec( strcmp(mobj.V_names,'wbill') ) = outstr.addvars.Lscale * (...
% 		exp(solvec(strcmp(mobj.Sol_names,'wB'))) .* mobj.Params.Lbar(1) + ...
% 		exp(solvec(strcmp(mobj.Sol_names,'wS'))) .* mobj.Params.Lbar(2) );
	% write different categories of variables in one row
	simnext=[point(1),outstr.exstvec',point(2:end),solvec',valvec',addvec,Reb_today,condvec];
	
	transmat = zeros(exnpt,2);
	transmat(:,1) = whatMtrans;
	transmat(:,2) = eItrans;
	
end

end