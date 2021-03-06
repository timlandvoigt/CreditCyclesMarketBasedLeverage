classdef OLGIntermediaryModel < model.DSGEModel
    
    properties (SetAccess=protected)
        NSTEN % number of endogenous state variables: Vfct.Ndim-1 
        NSTEX % number of exogenous state variables
        NSOL % number of solution vars
        NV % number of forecasting variables
        NADD % number of additional endogenous variables
        NCOND % number of conditional expectations of period-ahead variables        
        Sol_names % NSOLx1 cell array with names of solution variables
                   % must be in order of functions of Pfct 
        V_names % NVx1 cell array with names of forecasting variables
                   % must be in order of functions of Vfct      
        Sol_baseguess           
        En_names % NSTENx1 cell array with names of endog state variables           
        Ex_names % NSTEXx1 cell array with names of exog state variables           
        Add_names % NADDx1 cell array with names of additional vars
        Cond_names %NCONDx1 cell array with names of conditional variables        
        Params % structure array with parameters
        Exogenv % structure array with fields mtrans, pts_perm, pts_all
                % specifying exog. Markov processes
        Vfct % ApproxFunction object for iteration-relevant functions
        Pfct % ApproxFunction object for solution jump variables
        Tfct % ApproxFunction object for transition of state variable(s) (optional)
        Ffct
        Basegrid % base grid 
        nodes % for Gaussian integration
        weights % for Gaussian integration
        Addit % Additional Variables
    end    
    
    
    methods
        % constructor
        function obj=OLGIntermediaryModel(params,endogenv,exogenv,vfct,pfct,tfct,ffct,basegrid,nodes,weights,addvar)
            % call superclass constructor
            obj=obj@model.DSGEModel(params,endogenv,exogenv,vfct,pfct,tfct,ffct);
            obj.Basegrid=basegrid; 
            obj.NCOND=length(endogenv.condnames);          
            obj.Cond_names=reshape(endogenv.condnames,1,obj.NCOND);
            obj.nodes=nodes;
            obj.weights=weights;
            % Add this to get the extra variables, needed to compute the
            % PTI constraints.
            obj.Addit = addvar;
        end
        
        function [nextst,outstr]=calcStateTransition(obj,point,solvec,mode,varargin)
            
            params=obj.Params;
            exogenv=obj.Exogenv;
            if ~isempty(varargin)
                if length(varargin)>1
                    params=varargin{2};
                    exogenv=varargin{3};
                end
            end
            % unpack relevant params
            theta=params.theta;
            nu=params.nu;
            tau=params.tau;
            Hbar=params.Hbar;
            deltaH=params.deltaH;
            phi=params.phi;
            gamma=params.gamma;
            chi=params.chi;
            piM=params.piM;
            G=params.mu_G;
            rental = params.rental;
            tax = params.tax;
            if isfield(params,'MITshockloss')
                MITshockloss = params.MITshockloss;
            else
                MITshockloss =0;
            end
            
            if rental
                kappaY = params.kappaY;
            end
            
            
            % compute next period's states
            if isempty(varargin)
                State_next=obj.evaluateTrans(point);
            else
                if mode>1
                    % simulation mode with different Tfct
                    thistfct=varargin{1};
                    State_next=thistfct.evaluateAt(point);
                else
                    % solution mode
                    State_next=varargin{1};
                end
            end            
            
            
            % extract state variables
            exst=point(1);
            whatM=point(2);
            eI=point(3);
            Y=exogenv.pts_perm(exst,1);
            sig_epsY=exogenv.pts_perm(exst,2);
            sig_epsM=exogenv.pts_perm(exst,3);
            
            
            
            % extract solution variables
            R=exp(solvec(1));
            P=exp(solvec(2));
            q=exp(solvec(3));
            pY=exp(solvec(4));
            pM=exp(solvec(5));
            dM=exp(solvec(6));
            dY=exp(solvec(7));
            cM=exp(solvec(8));
            cY=exp(solvec(9));
            
            if rental
            hY=solvec(10);
            else
            hY=exp(solvec(10));
            end
            
            Itilde=exp(solvec(13));
            RebMort=solvec(17);
            
            % Market clearing
            DI = dY + dM;
            hM = Hbar - hY;
            
            if rental
            sY = (Hbar)/(1+(cM/cY));
            sM =  Hbar - sY; 
            rhoM = theta*cM/((1-theta)*sM);
            rhoY = rhoM;
            else
            sM = hM;
            sY = hY;    
            rhoM = theta*cM/((1-theta)*sM);
            rhoY = theta*cY/((1-theta)*sY);
            end
            
            % Wealth distribution
            Yhat = Y;
            I = (1-Itilde)/chi;
            xI = tau*eI - I;
            WMO = whatM + q+xI + pM;
            Phi = phi^(1/gamma)/(1+phi^(1/gamma));
            
            % Values for the old agents
            WO = piM * WMO;
            BO = Phi * WO;
            cO = WO - BO;
            
            % Change: the taxation and the demand for assets affects the
            % wealth of the agents.
            whatY = (1-deltaH)*P + BO + RebMort + Yhat - whatM - eI - MITshockloss;   
            
            WY = whatY + pY;
            WM = WMO - WO;
            fullendogvars=[whatM,eI]; 
            
            % Budget Constraints
            QYHH = cY + rhoY*sY + (P-rhoY)*hY + dY/R + pY - WY;
            QMHH = cM + rhoM*sM + (P-rhoM)*hM + q + pM + dM/R - WM;
            SM = (P-rhoM)*hM + q + pM + dM/R - QMHH;  % implied savings 
            SY = (P-rhoY)*hY + pY + dY/R - QYHH;  % implied savings 
     
            % tax adjustment for mortgage rate
            QY = QYHH/(1+tax);
            QM = QMHH/(1+tax);
    
            % Find the Rebate:
            %point_state = [exst,mobj.whatM,eI];
            Rebate_list=obj.evaluateForec(point)';
            Rebate_today = Rebate_list(exst);
            Y_star_today = Y + BO + RebMort + Rebate_today;
           
            
            if mode>0
                % simulation, mode contains number of next period's state
                exst=mode;
                nst=exogenv.exnpt;
                cind=exogenv.pts_all(exst,end);
                whatMnext=State_next(exst,1);
                eInext=State_next(nst+exst,1);
                nextst=[cind,whatMnext,eInext];                
            else
                % solution mode, compute next period's state variables for all
                % possible Markov states
                cind=exogenv.pts_all(:,end);
                nst=exogenv.exnpt;
                whatMnext=State_next(1:nst,1);
                eInext=State_next(nst+(1:nst),1);
                % matrix of next period states
                nextst=[cind,whatMnext,eInext];                
            end
            
            addvars=struct('WY',WY,...
                'WM',WM,...
                'DI',DI,...
                'hM',hM,...
                'sM',sM,...
                'rhoM',rhoM,...
                'sY',sY,...
                'rhoY',rhoY,...
                'QYHH',QYHH,...
                'QMHH',QMHH,...
                'QY',QY,...
                'QM',QM,...
                'I',I,...
                'xI',xI,...
                'SM', SM, ...
                'SY', SY, ...
                'cO',cO,...
                'BO',BO,...
                'Y',Y,...
                'end_yieldM',obj.Addit.end_yieldM,...
                'end_yieldY',obj.Addit.end_yieldY,...
                'Rebate_today',Rebate_today,...
                'Y_star_today',Y_star_today);
            
           
            
            outstr=struct;
            outstr.addvars=addvars;
            outstr.exstvec=[Y;sig_epsY;sig_epsM];
            outstr.fullendogvars=fullendogvars;
            
        end
        
        
        function [fx,J,V]=calcEquations(obj,exst,nextst,solvec,instr,mode,checkTrans,varargin)
            J =[]; % compatibility with older version
            params=obj.Params;
            % unpack params
            beta=params.beta;
            betaO=params.betaO;
            gamma=params.gamma;
            theta=params.theta;
            psiY=params.psiY;
            psiM=params.psiM;
            chi=params.chi;
            tau=params.tau;
            nu=params.nu;
            lambdaY=params.lambdaY;
            lambdaM=params.lambdaM;
            xi=params.xi;
            xiDWL=params.xiDWL;
            piM=params.piM;
            piY=params.piY;
            ebarY=params.ebarY;
            ebarM=params.ebarM;
            deltaH=params.deltaH;
            Hbar=params.Hbar;
            G=params.mu_G;
            phi=params.phi;
            delta_eta=params.delta_eta;
            zeta=params.zeta;
            alpha=params.alpha;
            tax=params.tax;
            hardconstr_LTV = params.LTV;
            hardconstr_PTI = params.PTI;
            %hardconstr=params.hardconstr;
            rental = params.rental;
            
            if rental
            kappaY = params.kappaY;
            end
            
            %%%%%%%%%%%%%%%%%%%%%%
            % Parameters LTV,PTI %
            %%%%%%%%%%%%%%%%%%%%%%
            
            thetaM_LTV=params.thetaM_LTV;
            thetaY_LTV=params.thetaY_LTV;    
            % Change the value of theta PTI here!
            thetaM_PTI=params.thetaM_PTI;
            thetaY_PTI=params.thetaY_PTI;
 
           
            % extract endogeous variables 
            R=exp(solvec(1));
            P=exp(solvec(2));
            q=exp(solvec(3));
            pY=exp(solvec(4));
            pM=exp(solvec(5));            
            dM=exp(solvec(6));
            dY=exp(solvec(7));
            
            if rental
            hY=solvec(10);
            else
            hY=exp(solvec(10));
            end
            
            mY=exp(solvec(11));
            mM=exp(solvec(12));
            vM=exp(solvec(14));
            vY=exp(solvec(15));
            muI=solvec(16);
            RebMort=solvec(17);
            
            var_number = 17;
            
            if rental
            var_number =  var_number + 1;   
            muI_rentY=solvec(var_number);
            end
            
            if hardconstr_LTV || hardconstr_PTI
                var_number =  var_number + 1;
                if  hardconstr_LTV && hardconstr_PTI
                    lamM_LTV=solvec(var_number);
                    lamY_LTV=solvec(var_number+1);
                    lamM_PTI=solvec(var_number+2);
                    lamY_PTI=solvec(var_number+3);
                    var_number = var_number + 3;
                elseif hardconstr_LTV && not(hardconstr_PTI)
                    lamM_LTV=solvec(var_number);
                    lamY_LTV=solvec(var_number+1);
                    var_number = var_number + 1;
                elseif not(hardconstr_LTV) && hardconstr_PTI
                    lamM_PTI=solvec(var_number);
                    lamY_PTI=solvec(var_number+1);
                    var_number = var_number + 1;
                end
            end
            

            exnpt=obj.Exogenv.exnpt;
            if checkTrans
                WMtrans_check = solvec(var_number+(1:exnpt));
                eItrans_check = solvec(var_number+exnpt+(1:exnpt));
            end
            % allocate result
            fx=zeros(obj.NSOL+checkTrans*2*exnpt,1);            
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % Multiplier transformations %
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            % Financial Intermediaries
            muIplus=max(0,muI)^3;
            muIminus=max(0,-muI)^3; 
            
            % Non-Negativity const. housing of the young agents.
            if rental
            muIplus_rentY=max(0,muI_rentY)^3;
            muIminus_rentY=max(0,-muI_rentY)^3;
            end
            
            if hardconstr_LTV || hardconstr_PTI
                if  hardconstr_LTV && hardconstr_PTI
                    % LTV
                    lamMplus_LTV=max(0,lamM_LTV)^3;
                    lamMminus_LTV=max(0,-lamM_LTV)^3;
                    lamYplus_LTV=max(0,lamY_LTV)^3;
                    lamYminus_LTV=max(0,-lamY_LTV)^3;
                    % PTI
                    lamMplus_PTI=max(0,lamM_PTI)^3;
                    lamMminus_PTI=max(0,-lamM_PTI)^3;
                    lamYplus_PTI=max(0,lamY_PTI)^3;
                    lamYminus_PTI=max(0,-lamY_PTI)^3;
                    
                elseif hardconstr_LTV && not(hardconstr_PTI)
                    % LTV
                    lamMplus_LTV=max(0,lamM_LTV)^3;
                    lamMminus_LTV=max(0,-lamM_LTV)^3;
                    lamYplus_LTV=max(0,lamY_LTV)^3;
                    lamYminus_LTV=max(0,-lamY_LTV)^3;
                    
                    % PTI
                    lamMplus_PTI=0;
                    lamYplus_PTI=0;
                
                elseif not(hardconstr_LTV) && hardconstr_PTI
                    % LTV
                    lamMplus_LTV=0;
                    lamYplus_LTV=0;
                    
                    % PTI
                    lamMplus_PTI=max(0,lamM_PTI)^3;
                    lamMminus_PTI=max(0,-lamM_PTI)^3;
                    lamYplus_PTI=max(0,lamY_PTI)^3;
                    lamYminus_PTI=max(0,-lamY_PTI)^3;
                end
            else
                %LTV
                lamMplus_LTV=0;
                lamYplus_LTV=0;
                
                % PTI
                lamMplus_PTI=0;
                lamYplus_PTI=0;
            end
            
           
            
            % extract some other state-dependent variables
            envec=instr.addvars;
            DI=envec.DI;
            rhoM=envec.rhoM;
            rhoY=envec.rhoY;
            QYHH=envec.QYHH;
            QMHH=envec.QMHH;
            QY=envec.QY;
            QM=envec.QM;
            I=envec.I;
            xI=envec.xI;
            SM=envec.SM;
            SY=envec.SY;
            hM=envec.hM;
            eI=instr.fullendogvars(2);
            BO=envec.BO;
            WM=envec.WM;
            WY=envec.WY;
            Y_star_today = envec.Y_star_today;
            
            % Variables needed for the DTI constraints:
            ym = (1-nu)*Y_star_today;
            yy = nu*Y_star_today;
            
            % probabilities and states to compute expectation terms
            prnext=obj.Exogenv.mtrans(exst,:);  
            Ynext=obj.Exogenv.pts_perm(:,1);
            sigeps_nextY=obj.Exogenv.pts_perm(:,2);
            sigeps_nextM=obj.Exogenv.pts_perm(:,3);
            mueps_nextY=-.5*sigeps_nextY.^2;
            mueps_nextM=-.5*sigeps_nextM.^2;
            
            % projection evaluation
            if nargin>=8 && (~checkTrans)
                % 8 or more arguments were passed.
                % Means expectations have already been computed (fast
                % solution mode)
                Pol_next = varargin{1};
                if mode==3
                    % credit surface computation
                    Yvar = varargin{2};                    
                    Mvar = varargin{3};
                end
            else
                % Nothing passed in ((conditional moments and EE errors in
                % simulation) or consistency equations
                if checkTrans
                    nextst(:,2) = WMtrans_check;
                    nextst(:,3) = eItrans_check;
                end
                Pol_next = obj.evaluateVal(nextst)';
                if size(Pol_next,1)==1
                    prnext=1;
                end
                forec = obj.evaluateForec([exst,instr.fullendogvars]);
                Pol_next=[Pol_next,forec];
            end

          
            q_next = Pol_next(:,1);
            I_next = Pol_next(:,2);
            vM_next = Pol_next(:,3);
            vY_next = Pol_next(:,4);
            P_next = Pol_next(:,5);
            pY_next = Pol_next(:,6);
            pM_next = Pol_next(:,7);
            BO_next = Pol_next(:,8);
            R_next = Pol_next(:,9);
            RebMort_next = Pol_next(:,10);
            Reb_next = Pol_next(:,11);
            
            xI_next=tau*nextst(:,3)-I_next;
           
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % HERE I start using the equations %
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  

            % Mortgage payoff and default for young
            Ystar_next = Ynext + Reb_next + BO_next + RebMort_next;
            ystar_next = (1-nu)*Ystar_next;
            wYYdef_next = dY./G + nu*Ystar_next + pY_next;
            
            % Change: eliminate BO_next
            wYMdef_next = dY./G + (delta_eta*ystar_next)/piY;
            epsbarYY_next = (mY/G - lambdaY*wYYdef_next)./((1-deltaH)*P_next*hY);
            epsbarYM_next = (mY/G - lambdaY*wYMdef_next)./((1-deltaH)*P_next*hY);
            
            % First for the individuals that stay young.
            [Feps_nextYY,feps_nextYY,FepsplusYY,FepsminusYY,...
                RFepsVFYY,RFepsFOC1YY,RFepsFOC2YY]= OLGIntermediaryModel.integEps(epsbarYY_next,mueps_nextY,sigeps_nextY,gamma,...
                                                          G*(1-deltaH)*P_next*hY/SY,(G*wYYdef_next-mY)/SY,obj.nodes,obj.weights,0);          
            [Feps_nextYM,feps_nextYM,FepsplusYM,FepsminusYM,...
                RFepsVFYM,RFepsFOC1YM,RFepsFOC2YM]= OLGIntermediaryModel.integEps(epsbarYM_next,mueps_nextY,sigeps_nextY,gamma,...
                                                          G*(1-deltaH)*P_next*hY/SY,(G*wYMdef_next-mY)/SY,obj.nodes,obj.weights,0); 
            
            %%%%%%%%%%%%%%%%%%%%%%%
            % Payoffs Calculations%
            %%%%%%%%%%%%%%%%%%%%%%%
            
            % YOUNG
            Def_nextY = (1-piY)*Feps_nextYY + piY*Feps_nextYM;
            FepsminusY = (1-piY)*FepsminusYY + piY*FepsminusYM;
            Payoff_nextY = (1-Def_nextY).*mY./G + (1-deltaH)*P_next*hY.*(alpha + (1-xi)*FepsminusY);
            
            wMdef_next = dM./G + q_next+xI_next + ((1-delta_eta)*ystar_next+pM_next);
            epsbarM_next = (mM/G - lambdaM*wMdef_next)./((1-deltaH)*P_next*hM);
            
            [Feps_nextM,feps_nextM,FepsplusM,FepsminusM,...
                RFepsVFM,RFepsFOC1M,RFepsFOC2M]= OLGIntermediaryModel.integEps(epsbarM_next,mueps_nextM,sigeps_nextM,gamma,...
                                                          G*(1-deltaH)*P_next*hM/SM,(G*wMdef_next-mM)/SM,obj.nodes,obj.weights,0);
            
            Payoff_nextM = (1-Feps_nextM).*mM./G + (1-deltaH)*P_next*hM.*(alpha + (1-xi)*FepsminusM);
            
            %%%%%%%%%%%%%%%%%%
            % Financial Int. %
            %%%%%%%%%%%%%%%%%%

            eInext = Payoff_nextY + Payoff_nextM - DI./G;
            [min_eInext,minstate] = min(eInext);
            min_PayoffY=Payoff_nextY(minstate);
            min_PayoffM=Payoff_nextM(minstate);
             
            %%%%%%%%%%%%%%%%%%%%%%
            % Compute Middle SDF %
            %%%%%%%%%%%%%%%%%%%%%%
            
            Phi = (1/(1+ phi^(1/gamma)))^(1-gamma) + phi*(phi^(1/gamma)/(1+ phi^(1/gamma)))^(1-gamma) * betaO;
            RMdef_next = G*(1-lambdaM)*wMdef_next/SM;
            ThetaZ = ((1-theta)^(1-theta) * (theta/rhoM)^theta)^(1-gamma);
            AZ = psiM*(dM/SM)^(1-gamma) + ...
                beta*prnext*((RFepsVFM + Feps_nextM.*RMdef_next.^(1-gamma)) .* (piM*Phi + (1-piM)*vM_next));
            BZ = AZ^(1/gamma)/( AZ^(1/gamma) + ThetaZ^(1/gamma));

            MMnodef_m = beta * BZ^(-gamma)/vM * RFepsFOC1M .* (piM*Phi + (1-piM)*vM_next);
            MMnodef_h = beta * BZ^(-gamma)/vM * RFepsFOC2M .* (piM*Phi + (1-piM)*vM_next);
            MMdef = beta * BZ^(-gamma)/vM * RMdef_next.^(-gamma) .* (piM*Phi + (1-piM)*vM_next);

            SDFM = MMnodef_m + (1-lambdaM)*Feps_nextM.*MMdef;
            
            % With this we can compute the SDF of the financial
            % intermediares.
            SDFI = SDFM.*(tau + (1-tau)./(1-chi*I_next))*(1-chi*I);
                       
            
            %%%%%%%%%%%%%%%%%%%%%%
            % Derivatives of q_M %
            %%%%%%%%%%%%%%%%%%%%%%
            

            QMdiff_term = feps_nextM.*(mM - G*(1-xi)*(1-deltaH)*hM*P_next.*epsbarM_next);
            QMm_term = QMdiff_term./(G*(1-deltaH)*P_next*hM) + (1-deltaH)*P_next*hM/mM .*(alpha + (1-xi)*FepsminusM);
            QMh_term = G*(1-deltaH)*P_next.*(alpha + (1-xi)*FepsminusM + QMdiff_term.*epsbarM_next./(G*(1-deltaH)*P_next*hM));
            QMd_term = lambdaM*QMdiff_term./(G*(1-deltaH)*P_next*hM);
            QMb_term = QMd_term.*(q_next+xI_next)*G;
            QMeta_term = QMd_term.*((1-delta_eta)*ystar_next+pM_next)*G;
            
            QMh = (1+tax)*((1-ebarM)*muIplus*QMh_term(minstate) + prnext*(SDFI.*QMh_term))/(1+zeta);
            QMm = (1+tax)*((1-ebarM)*muIplus*QMm_term(minstate) + prnext*(SDFI.*QMm_term))/(1+zeta);
            QMd = (1+tax)*((1-ebarM)*muIplus*QMd_term(minstate) + prnext*(SDFI.*QMd_term))/(1+zeta);
            QMb = (1+tax)*((1-ebarM)*muIplus*QMb_term(minstate) + prnext*(SDFI.*QMb_term))/(1+zeta);            
            QMeta = (1+tax)*((1-ebarM)*muIplus*QMeta_term(minstate) + prnext*(SDFI.*QMeta_term))/(1+zeta);            
            
            %%%%%%%%%%%%%%%%%%%%%%
            % Compute Young SDF %
            %%%%%%%%%%%%%%%%%%%%%%
            
            RYYdef_next = G*(1-lambdaY)*wYYdef_next/SY;
            RYMdef_next = G*(1-lambdaY)*wYMdef_next/SY;
            ThetaY = ((1-theta)^(1-theta) * (theta/rhoY)^theta)^(1-gamma);
            AY = psiY*(dY/SY)^(1-gamma) ...
                 + beta*(1-piY)*prnext*((RFepsVFYY + Feps_nextYY.*RYYdef_next.^(1-gamma)) .* vY_next) ...
                 + beta*piY*prnext*((RFepsVFYM + Feps_nextYM.*RYMdef_next.^(1-gamma)) .* vM_next);
            BY = AY^(1/gamma)/( AY^(1/gamma) + ThetaY^(1/gamma));
 
            MYnodef_m = beta * BY^(-gamma)/vY * ((1-piY)*RFepsFOC1YY.*vY_next + piY*RFepsFOC1YM.*vM_next);
            MYnodef_h = beta * BY^(-gamma)/vY * ((1-piY)*RFepsFOC2YY.*vY_next + piY*RFepsFOC2YM.*vM_next);
            MYdef = beta * BY^(-gamma)/vY * ((1-piY)*Feps_nextYY.*RYYdef_next.^(-gamma).*vY_next + piY*Feps_nextYM.*RYMdef_next.^(-gamma).*vM_next );
            
            SDFY = MYnodef_m + (1-lambdaY)*MYdef;
            
            MYY = beta * BY^(-gamma)/vY * (1-piY)*(RFepsFOC1YY + (1-lambdaY)*Feps_nextYY.*RYYdef_next.^(-gamma)) .* vY_next;
            MYM = beta * BY^(-gamma)/vY * piY*(RFepsFOC1YM + (1-lambdaY)*Feps_nextYM.*RYMdef_next.^(-gamma)) .* vM_next;
            
            %%%%%%%%%%%%%%%%%%%%%%
            % Derivatives of q_Y %
            %%%%%%%%%%%%%%%%%%%%%%

            epshatYY = mY - G*(1-xi)*(1-deltaH)*hY*P_next.*epsbarYY_next;
            epshatYM = mY - G*(1-xi)*(1-deltaH)*hY*P_next.*epsbarYM_next;
            fhat_nextY = (piY*feps_nextYM.*epshatYM + (1-piY)*feps_nextYY.*epshatYY)./(G*(1-deltaH)*P_next*hY);
            

            QYm_term = fhat_nextY + (1-deltaH)*P_next*hY/mY .* (alpha + (1-xi)*FepsminusY);
            QYd_term = lambdaY*fhat_nextY;
            fhh_nextY = (piY*feps_nextYM.*epshatYM.*epsbarYM_next  ...
                         + (1-piY)*feps_nextYY.*epshatYY.*epsbarYY_next)./(G*(1-deltaH)*P_next*hY);         
            QYh_term = G*(1-deltaH)*P_next.*(alpha + FepsminusY*(1-xi)+ fhh_nextY);
            
            fhhh_nextY = (piY*((delta_eta*ystar_next)/piY).*feps_nextYM.*epshatYM ...
                    + (1-piY)*feps_nextYY.*(nu*Ystar_next + pY_next).*epshatYY)./(G*(1-deltaH)*P_next*hY);
            QYeta_term = G*lambdaY*fhhh_nextY;
            QYeta = (1+tax)*((1-ebarY)*muIplus*QYeta_term(minstate) + prnext*(SDFI.*QYeta_term))/(1+zeta);
            QYm = (1+tax)*((1-ebarY)*muIplus*QYm_term(minstate) + prnext*(SDFI.*QYm_term))/(1+zeta);
            QYh = (1+tax)*((1-ebarY)*muIplus*QYh_term(minstate) + prnext*(SDFI.*QYh_term))/(1+zeta);
            QYd = (1+tax)*((1-ebarY)*muIplus*QYd_term(minstate) + prnext*(SDFI.*QYd_term))/(1+zeta);
            
            
            %%%%%%%%%%%%%%%%%%%%%%
            % FOC's  Middle aged %
            %%%%%%%%%%%%%%%%%%%%%%
            
            % Effective returns on housing
            Rh_eff_next = G*(1-deltaH)*P_next/(P - rhoM - QMh);
            % Effective returns to mortgages 
            Rm_eff_next = 1/(QMHH/mM - QMm);
            % Effective returns to deposits. 
            Rd_eff = R/(1-R*QMd);
            % Effective returns to equity.
            Re_eff_next = G*(q_next+xI_next)/(q - QMb);
            % Effective returns to endowment.
            Reta_eff_next = G*((1-delta_eta)*ystar_next+pM_next)/(pM - QMeta);
            % This is the first element of equation (42)
            cyieldM = psiM*(dM*BZ/SM)^(-gamma)/vM;
            

            % The FOC written as no-excess return equations.
            fx(1) = prnext*(SDFM.*Re_eff_next - MMnodef_h.*Rh_eff_next) ...
                - lamMplus_LTV*thetaM_LTV*(BZ^(-gamma)/vM)*P/(P - rhoM - QMh);
            
            fx(2) = Rd_eff*cyieldM + prnext*(SDFM.*(Rd_eff - Re_eff_next));
            
            fx(3) = prnext*(MMnodef_h.*Rh_eff_next - MMnodef_m.*Rm_eff_next) ...
                +lamMplus_LTV*(BZ^(-gamma)/vM)*(thetaM_LTV*P/(P - rhoM - QMh) - Rm_eff_next)...
                -lamMplus_PTI*(BZ^(-gamma)/vM)*Rm_eff_next;
            
            payoffend_m = ym;
            
            fx(4) = prnext*(SDFM.*(Re_eff_next - Reta_eff_next))...
                -lamMplus_PTI*payoffend_m*thetaM_PTI*(BZ^(-gamma)/vM)*(1/(pM - QMeta));
            
            % Defintion of the savings for the middle aged.
            fx(5) = SM - BZ*WM;
                         
            % This is the defintion of v_m(Z) in the appenix.
            fx(6) = vM - ThetaZ*(1-BZ)^(1-gamma) - AZ*BZ^(1-gamma);
            
            %%%%%%%%%%%%%%%%%%%%%%%%
            % FOC's  Intermdiaries %
            %%%%%%%%%%%%%%%%%%%%%%%%
            
            fx(7)=1/R - muIplus/G - prnext*SDFI;
            fx(8)=(1+zeta)*QY - (1-ebarY)*muIplus*min_PayoffY - G*prnext*(SDFI.*Payoff_nextY);
            fx(9)=(1+zeta)*QM - (1-ebarM)*muIplus*min_PayoffM - G*prnext*(SDFI.*Payoff_nextM);
            
            %%%%%%%%%%%%%%%%%%%%%%%%
            % FOC's  Young agents  %
            %%%%%%%%%%%%%%%%%%%%%%%%
            
            % Effective returns on housing
            RYh_eff_next = G*(1-deltaH)*P_next/(P - rhoY - QYh);
            
            % Effective returns on mortgages.
            RYm_eff_next = 1/(QYHH/mY - QYm);
            
            % Effective returns on deposits (not an expectation!)
            RYd_eff = R/(1-R*QYd);
            
            % Effective returns of the endwoment asset if you don't age.
            RYYeta_eff_next = G*(nu*Ystar_next + pY_next)/(pY - QYeta);
            
            % Change: eliminate BO_next
            % Effective returns of the endwoment asset if you age.
            RYMeta_eff_next = G*((delta_eta*ystar_next)/piY)/(pY - QYeta);
            
            % Auxiliary variable: first part of equation (42)
            cyieldY = psiY*(dY*BY/SY)^(-gamma)/vY;
            
            % The FOC written as no-excess return equations.
            fx(10)= cyieldY*RYd_eff + prnext*(SDFY.*RYd_eff - MYnodef_m.*RYm_eff_next)...
                - lamYplus_LTV*(BY^(-gamma)/vY)*RYm_eff_next...
                - lamYplus_PTI*(BY^(-gamma)/vY)*RYm_eff_next;
            
            if rental
            fx(11)= prnext*(MYnodef_m.*RYm_eff_next - MYnodef_h.*RYh_eff_next) ...
                + lamYplus_LTV*(BY^(-gamma)/vY)*(RYm_eff_next - thetaY_LTV*P/(P - rhoY - QYh))...
                + lamYplus_PTI*(BY^(-gamma)/vY)*RYm_eff_next...
                - (muIplus_rentY+(kappaY*(hY)^(-gamma))) * (BY^(-gamma)/vY)*(1/(P - rhoY - QYh));
            else
            fx(11)= prnext*(MYnodef_m.*RYm_eff_next - MYnodef_h.*RYh_eff_next) ...
                + lamYplus_LTV*(BY^(-gamma)/vY)*(RYm_eff_next - thetaY_LTV*P/(P - rhoY - QYh))...
                + lamYplus_PTI*(BY^(-gamma)/vY)*RYm_eff_next;    
            end
                      
           payoffend_y = yy;

           if rental 
           fx(12)= prnext*(MYnodef_h.*RYh_eff_next - (MYY.*RYYeta_eff_next + MYM.*RYMeta_eff_next)) ...
                + lamYplus_LTV*(BY^(-gamma)/vY)*thetaY_LTV*P/(P - rhoY - QYh)...
                - lamYplus_PTI*(BY^(-gamma)/vY)*thetaY_PTI*(payoffend_y/(pY - QYeta))...
                + (muIplus_rentY+(kappaY*(hY)^(-gamma))) * (BY^(-gamma)/vY)*(1/(P - rhoY - QYh));
           else
           fx(12)= prnext*(MYnodef_h.*RYh_eff_next - (MYY.*RYYeta_eff_next + MYM.*RYMeta_eff_next)) ...
                + lamYplus_LTV*(BY^(-gamma)/vY)*thetaY_LTV*P/(P - rhoY - QYh)...
                - lamYplus_PTI*(BY^(-gamma)/vY)*thetaY_PTI*(payoffend_y/(pY - QYeta));    
           end
            
            % Defintion of the savings for the middle aged.
            fx(13)= SY - BY*WY;
            
            % This is the defintion of v_y(Z) in the appenix.
            fx(14)= vY - ThetaY*(1-BY)^(1-gamma) - AY*BY^(1-gamma);
            
            % CC for the Intermediary
            fx(15)=(1-ebarY)*min_PayoffY+(1-ebarM)*min_PayoffM - DI/G - muIminus; 
            
            % Budget Constraint of the Intermediary.
            fx(16)= (1+zeta)*(QY + QM) - ((1-tau)*eI + I - chi*I^2/2 + DI/R);
            
            % Rebates
            fx(17)= RebMort - ( zeta*(QY+QM) - (QMHH-QM) - (QYHH-QY) );
            
            eq_number = 17;
            
            if rental
            eq_number = eq_number+1;    
            fx(eq_number) = hY - muIminus_rentY;
            end
            
            if hardconstr_LTV || hardconstr_PTI
                eq_number = eq_number+1; 
                if  hardconstr_LTV && hardconstr_PTI
                    % LTV and PTI
                    fx(eq_number)=thetaM_LTV*P*hM - mM-lamMminus_LTV;
                    fx(eq_number+1)=thetaY_LTV*P*hY - mY-lamYminus_LTV;
                    fx(eq_number+2)=thetaM_PTI*payoffend_m - mM-lamMminus_PTI;
                    fx(eq_number+3)=thetaY_PTI*payoffend_y - mY-lamYminus_PTI;
                    eq_number = eq_number+3;
                elseif hardconstr_LTV && not(hardconstr_PTI)
                    % LTV and PTI
                    fx(eq_number)=thetaM_LTV*P*hM - mM-lamMminus_LTV;
                    fx(eq_number+1)=thetaY_LTV*P*hY - mY-lamYminus_LTV;
                    eq_number = eq_number+1;

                elseif not(hardconstr_LTV) && hardconstr_PTI
                    % LTV and PTI
                    
                    fx(eq_number)=thetaM_PTI*payoffend_m - mM-lamMminus_PTI;
                    fx(eq_number+1)=thetaY_PTI*payoffend_y - mY-lamYminus_PTI;
                    eq_number = eq_number+1;
                end
            end
                        
            %%%%%%%%%%%%%%%
            % Transitions % 
            %%%%%%%%%%%%%%%
            
            % Aggregate income of the young defaulters:
            WYYdef_next = (1-piY)*(dY./G + nu*Ystar_next + pY_next);
            % Change: eliminate BO_next
            WYMdef_next = piY*dY./G + (delta_eta*ystar_next);
            % Aggregate income of the middle aged defaulters:
            WMdef_next = dM./G + q_next+xI_next + (1-delta_eta)*ystar_next + pM_next;
            

            RebateY = (1-piY)*lambdaY*Feps_nextYY.*WYYdef_next + lambdaY*Feps_nextYM.*WYMdef_next + (xi-xiDWL)*(1-deltaH)*FepsminusY.*P_next*hY;
            RebateM = (xi-xiDWL)*(1-deltaH)*FepsminusM.*P_next*hM + lambdaM*Feps_nextM.*WMdef_next;
            RebateI = -alpha*P_next*Hbar;
            
            % Change: eliminate BO_next
            WYMnext = (delta_eta*ystar_next + piY*dY/G).*(1-lambdaY.*Feps_nextYM)...
                + piY*((1-deltaH)*FepsplusYM.*P_next*hY - (1-Feps_nextYM)*mY/G);  
            
            % Total wealth of the middle age and old age people:
            WMnext = WYMnext + (1-deltaH)*FepsplusM.*P_next*hM - (1-Feps_nextM)*mM/G + (1-lambdaM*Feps_nextM).*WMdef_next;
            
            % This is the value of the state variable:
            %Yhat_next = Ynext - lumptax_next;
            WMtrans = WMnext - (q_next+xI_next) - pM_next;
            eItrans = eInext;
            
            if checkTrans
                fx(eq_number+(1:exnpt)) = WMtrans - WMtrans_check;
                fx(eq_number+exnpt+(1:exnpt)) = eItrans - eItrans_check;
            end
                                
            
            V=cell(3,1);
            % marginal value functions
            
            if mode==1
                % Output new values for time iteration
                Vnext=zeros(obj.Vfct.Nof,1);
                Vnext(1)=q;
                Vnext(2)=I;
                Vnext(3)=vM;
                Vnext(4)=vY;
                Vnext(5)=P;
                Vnext(6)=pY;
                Vnext(7)=pM;
                Vnext(8)=BO;
                Vnext(9)=R;
                Vnext(10)=RebMort;
                
                V{1}=Vnext;
                % state transition
                V{2}=[WMtrans; eItrans]';
                % forecasting function for rebates
                V{3}=(RebateY+RebateM+RebateI)';
                
        
            elseif mode==2
                
                % Evaluation during simulation. Output conditional
                % variables
                
                % SDFs
                SDF.SDFM = SDFM;                             
                SDF.SDFI = SDFI;
                SDF.SDFY = SDFY;
               
                % Define returns 
                NetPayoff_nextM = Payoff_nextM - alpha.*(1-deltaH).*P_next*hM;
                NetPayoff_nextY = Payoff_nextY - alpha.*(1-deltaH).*P_next*hY;
                retMM = G*NetPayoff_nextM/QM;
                retMY = G*NetPayoff_nextY/QY;
                retH = G.*(1-deltaH)*P_next/(P-rhoM);
                retBext = G.*(q_next + xI_next)/q;
                retBin =(G*(Payoff_nextM+Payoff_nextY) - DI)./((1+zeta)*(QY+QM) - DI/R); 
                expR_MM = prnext * retMM;
                expR_MY = prnext * retMY;
                expR_H = prnext * retH;
                expR_Bext = prnext * retBext;
                expR_Bin = prnext * retBin;
                % Expected default rates                
                EDefY = prnext * ((1-piY)*Feps_nextYY + piY*Feps_nextYM);
                EDefM = prnext * Feps_nextM;
                
                
                condvars = struct('expR_MM',expR_MM, ...
                                  'expR_MY',expR_MY, ...
                                  'expR_H',expR_H,...
                                  'expR_Bext',expR_Bext,...
                                  'expR_Bin',expR_Bin,...
                                  'min_eInext',min_eInext,...
                                  'cyieldM',cyieldM,...
                                  'cyieldY',cyieldY,...
                                  'EDefY',EDefY,...
                                  'EDefM',EDefM);                        
                
%                Wtrans.WY = WYtrans;                          
                Wtrans.WM = WMtrans;                          
                Wtrans.eI = eItrans;
                V = {condvars,Wtrans,SDF};
                
            elseif mode==3
               % compute credit surface 
               % first young
               npty=size(Yvar,1);
               QYmat=zeros(npty,4);
               for yi=1:npty
                   hYi = Yvar(yi,1);
                   levYi = Yvar(yi,2);
                   mYi=P*hYi*levYi;
                   levYWi = Yvar(yi,3);  
                   wYYdef_next = mYi/levYWi;
                   wYMdef_next = mYi/levYWi;
                   epsbarYY_next = (mYi/G - lambdaY*wYYdef_next)./((1-deltaH)*P_next*hYi);
                   epsbarYM_next = (mYi/G - lambdaY*wYMdef_next)./((1-deltaH)*P_next*hYi);
                   [Feps_nextYY,~,~,FepsminusYY]= OLGIntermediaryModel.integEps(epsbarYY_next,mueps_nextY,sigeps_nextY,gamma,[],[],obj.nodes,obj.weights,1);
                   [Feps_nextYM,~,~,FepsminusYM]= OLGIntermediaryModel.integEps(epsbarYM_next,mueps_nextY,sigeps_nextY,gamma,[],[],obj.nodes,obj.weights,1);
                   Def_nextY = (1-piY)*Feps_nextYY + piY*Feps_nextYM;
                   FepsminusY = (1-piY)*FepsminusYY + piY*FepsminusYM;
                   Payoff_nextYi = (1-Def_nextY).*mYi./G + (1-xi)*(1-deltaH)*FepsminusY.*P_next*hYi;
                   % price this mortgage
                   QYmat(yi,1)=(1-ebarY)*muIplus*Payoff_nextYi(minstate) + G*prnext*(SDFI.*Payoff_nextYi);
                   QYmat(yi,2)=mYi;  % mortgage
                   QYmat(yi,3)=mYi/QYmat(yi,1); % spread
                   QYmat(yi,4)=mYi/QYmat(yi,1) - R; % spread
               end
               % then middle aged
               nptm=size(Mvar,1);
               QMmat=zeros(nptm,4);
               for mi=1:nptm 
                   hMi = Mvar(mi,1);
                   levMi = Mvar(mi,2);
                   mMi=P*hMi*levMi;
                   levMWi = Mvar(mi,3);    
                   wMdef_next = mMi/levMWi;
                   epsbarM_next = (mMi/G - lambdaM*wMdef_next)./((1-deltaH)*P_next*hMi);
                   [Feps_nextM,~,~,FepsminusM]= OLGIntermediaryModel.integEps(epsbarM_next,mueps_nextM,sigeps_nextM,gamma,[],[],obj.nodes,obj.weights,1);
                   Payoff_nextMi = (1-Feps_nextM).*mMi./G + (1-xi)*(1-deltaH)*FepsminusM.*P_next*hMi;
                   % price this mortgage
                   QMmat(mi,1)=(1-ebarM)*muIplus*Payoff_nextMi(minstate) + G*prnext*(SDFI.*Payoff_nextMi);  
                   QMmat(mi,2)=mMi;
                   QMmat(mi,3)=mMi/QMmat(mi,1);
                   QMmat(mi,4)=mMi/QMmat(mi,1) - R;
               end
               V={QYmat,QMmat};
            end            
        end
  
        
        function [errmat,solmat,condmat,Wshtrans,SDFmat]=calcEEError(obj,pointmat)
            % function to compute Euler equation error at points in state
            % space given by pointmat
            nst=size(obj.Exogenv.pts_all,1);
            
            errmat=zeros(size(pointmat,1),obj.Pfct.Nof);
            solmat=zeros(size(errmat));
            condmat=zeros(size(pointmat,1),obj.NCOND);
            SDFmat=zeros(size(pointmat,1),3*nst);
            Wshtrans=zeros(size(pointmat,1),obj.NSTEN*nst);
            
            evaluatePol = @(point)obj.evaluatePol(point);
            calcStateTransition = @(point,soltmp)obj.calcStateTransition(point,soltmp,0);
            calcEquations = @(exst,nextst,soltmp,outstr)obj.calcEquations(exst,nextst,soltmp,outstr,2,0);    
            params=obj.Params;
            % Should be parfor. Use for when debugging only
            parfor i=1:size(errmat,1)
%            for i=1:size(errmat,1)
                point=pointmat(i,:);
                soltmp=evaluatePol(point);
                % transition
                [nextst,outstr]=calcStateTransition(point,soltmp);
                % equations
                [fx,~,V]=calcEquations(point(1),nextst,soltmp,outstr);                                
                R=exp(soltmp(1));
                P=exp(soltmp(2));
                vM=exp(soltmp(14));
                vY=exp(soltmp(15));
                Reb=soltmp(17);
                QY=outstr.addvars.QY;
                QM=outstr.addvars.QM;
                DI=outstr.addvars.DI;
                SM=outstr.addvars.SM;
                SY=outstr.addvars.SY;
                
                if params.rental
                normvec=[R,R,R,R,SM,vM,1/R,QY,QM,R,R,R,SY,vY,DI,1,Reb,1];
                else
                normvec=[R,R,R,R,SM,vM,1/R,QY,QM,R,R,R,SY,vY,DI,1,Reb];
                end

                if params.PTI || params.LTV 
                    if  params.LTV && params.PTI
                        % LTV and PTI
                       normvec=[normvec,P,P,P,P];

                    elseif params.LTV && not(params.PTI)
                        % LTV and PTI
                        normvec=[normvec,P,P];
                    elseif not(params.LTV) && params.PTI
                        % LTV and PTI
                        normvec=[normvec,P,P];
                    end
                end

                condvars=V{1};
%                WYtrans=V{2}.WY;                
                WMtrans=V{2}.WM;                
                eItrans=V{2}.eI;
                SDFI=V{3}.SDFI;
                SDFM=V{3}.SDFM;
                SDFY=V{3}.SDFY;
                errmat(i,:)=fx'./normvec;
                solmat(i,:)=soltmp';
                condmat(i,:)=model.DSGEModel.structToVec(condvars)';
                Wshtrans(i,:) = [WMtrans', eItrans'];
                SDFmat(i,:) = [SDFM',SDFI',SDFY'];
            end
            
        end
        
        
        function [QYmat,QMmat,baseQ]=calcCSurf(obj,Yvar,Mvar)
            % function to compute credit surface at equilibrium
            % intermediary SDF
            npt=obj.Vfct.SSGrid.Npt;         
            ynpt=size(Yvar,1);
            mnpt=size(Mvar,1);
            NDIM=obj.Vfct.SSGrid.Ndim;
            exnpt=size(obj.Exogenv.pts_perm,1);            
            pointmat=obj.Vfct.SSGrid.Pointmat;
            
            QYmat=zeros(npt,ynpt,4);
            QMmat=zeros(npt,mnpt,4);
            baseQ=zeros(npt,2);

            % transitions
            transmat=obj.evaluateTrans(pointmat)';
            forecmat=obj.evaluateForec(pointmat)';
            transpts=reshape([repmat(1:exnpt,npt,1),transmat],npt*exnpt,NDIM);
            Vtrans = obj.evaluateVal(transpts)';
            % add forecasting points
            forecpts=reshape(forecmat,npt*exnpt,1); % stack columns representing different states
            Vtrans = [Vtrans, forecpts];
            % index matching for transitions
            refidx=kron((1:npt)',ones(1,exnpt));
            refidx=refidx(:);
                       
            evaluatePol = @(point)obj.evaluatePol(point);
            calcStateTransition = @(point,soltmp,trans)obj.calcStateTransition(point,soltmp,0,trans);
            calcEquations = @(exst,nextst,soltmp,outstr,vtrans)obj.calcEquations(exst,nextst,soltmp,outstr,3,0,vtrans,Yvar,Mvar);  
            
            % Should be parfor. Use for when debugging only
            parfor i=1:npt
            %for i=1:npt
                point=pointmat(i,:);
                trans=transmat(i,:)';
                vtrans=Vtrans(refidx==i,:);
                soltmp=evaluatePol(point);
                % transition
                [nextst,outstr]=calcStateTransition(point,soltmp,trans);
                % equations
                [~,~,V]=calcEquations(point(1),nextst,soltmp,outstr,vtrans);                                
                QYmat(i,:,:) = V{1};
                QMmat(i,:,:) = V{2};
                baseQ(i,:)=[outstr.addvars.QY,outstr.addvars.QM];
            end
            
        end
        
        % simulate model; overwrite method from superclass 
         function [simseries,varnames,errmat,Wshtrans,SDFmat,nextst]=simulate(obj,NT,NTini,inistvec,simerror,shmat,varargin)
            
             % to check you have the correct number of state variables.
             if length(inistvec)~=obj.Vfct.SSGrid.Ndim
                error('inistvec must be vector of length SSGrid.Ndim');
             end
            
            
            NTtot=NT+NTini;
            simseries=zeros(NTtot,1+obj.NSTEX+obj.NSTEN+obj.NSOL+obj.NV+(obj.NADD)+1);
            ffnext=zeros(obj.Exogenv.exnpt,1); % rebate forecasting function
            
            % if shock matrix wasn't passed in, create it
            preset_path = false;
            if nargin<6
                rng(10,'twister');
                shmat=lhsdesign(NTtot,1,'criterion','correlation'); % realizations of shocks for Markov state
            else 
                preset_path=isinteger(shmat);
            end        
            
            point=inistvec;
            % point = startpt_vec;
                       
            pointmat=zeros(NTtot,length(point));
            % NTtot = NT_sim + NT_ini;
            
            for t=1:NTtot
               pointmat(t,:)=point; 
               exst=point(1);
               rebate=ffnext(exst); % get rebate for this period here
               
                % next period's exog. state
               if preset_path
                    exnext=shmat(t);
               else
                    transprob=cumsum(obj.Exogenv.mtrans(exst,:));
                    % Finding the value of the exogenous shock for the next
                    % period.
                    exnext=find(transprob-shmat(t)>0,1,'first');
               end
               
               % transition to next period
               
               % As before, use interpolation to get the value of the
               % policy function, the value functions and the rebate at the
               % initial point.
               solvec=obj.evaluatePol(point);
               valvec=obj.evaluateVal(point)';     
               ffnext=obj.evaluateForec(point)';
               [nextst,outstr]=obj.calcStateTransition(point,solvec,exnext,varargin{:});
               
               if isempty(outstr.addvars.Rebate_today)
                    outstr.addvars.Rebate_today = 0;
               end
               
                if isempty(outstr.addvars.Y_star_today)
                    outstr.addvars.Y_star_today = 0;
               end
           
               
               addvec=model.DSGEModel.structToVec(outstr.addvars)';
               % write different categories of variables in one row
               % simseries=zeros(NTtot,1+obj.NSTEX+obj.NSTEN+obj.NSOL+obj.NV+obj.NADD+1);
               
               
               simnext=[point(1),outstr.exstvec',point(2:end),solvec',valvec,addvec,rebate];
               if length(simnext)~=size(simseries,2)
                    disp('problem');
               end
               simseries(t,:)=simnext;
               % Update the state for the next simulated period.
               point=nextst;
            end     
            
            simseries=simseries(NTini+1:end,:);
            varnames=[{'exst'}, obj.Ex_names, obj.En_names, obj.Sol_names, obj.V_names, obj.Add_names, {'Reb'}];
            
            
            errmat=[];
            Wshtrans=[];%zeros(size(obj.Exogenv.mtrans(simseries(:,1),:)));
            SDFmat=[];%zeros(size(obj.Exogenv.mtrans(simseries(:,1),:)));
            if simerror
                [errmat,~,condmat,Wshtrans,SDFmat]=obj.calcEEError(pointmat);
                errmat=errmat(NTini+1:end,:);
                condmat=condmat(NTini+1:end,:);
                Wshtrans=Wshtrans(NTini+1:end,:);
                SDFmat=SDFmat(NTini+1:end,:);
                simseries=[simseries,condmat];
                varnames=[varnames,obj.Cond_names];
            else
                simseries=[simseries,zeros(NT,obj.NCOND)];                
                % varnames=[varnames,strcat(obj.Cond_names,'_nan')];
                varnames=[varnames,strcat(obj.Cond_names)];
            end
         end
        
         
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % This polIter function is the one used in main_run_exper.m %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        function [mobj,failedPoints,dist,distT]=polIter(mobj,MAXIT,revisitFailed,printmode,damp,tol_avg,checkTrans)
            gridSt=mobj.Vfct.SSGrid.Pointmat; % use points from BaseGrid here
            NPT=mobj.Vfct.SSGrid.Npt;
            NDIM=mobj.Vfct.SSGrid.Ndim;
            exnpt=size(mobj.Exogenv.pts_perm,1);
            
            % initialize
            % QUESTION!!!
            resmat=mobj.evaluatePol(gridSt)';
            resmat_prev=resmat;
            
            % split up matrix of points for better output
            gr_points=cell(exnpt,1);
            gr_index=cell(exnpt,2);
            
            for i=1:exnpt
                grinlog=(gridSt(:,1)==i);
                grind=find(grinlog);
                gr_points{i}=gridSt(grinlog,:);
                gr_index{i,1}=grinlog;
                gr_index{i,2}=grind;
            end
            
            
            % value function
            
            % The function evaluateVal is in the script DSGEModel.m
         
            VF=mobj.evaluateVal(gridSt)';
            VFnext=zeros(size(VF));
            TF=mobj.evaluateTrans(gridSt)';
            TFnext=zeros(size(VF,1),mobj.Tfct.Nof);
            FF=mobj.evaluateForec(gridSt)';
            FFnext=zeros(size(VF,1),mobj.Ffct.Nof);
            resmat_old = resmat; % for dampening
                        
            % control flags
            iter=0;
                        
            disp(' ');
            disp('Starting main loop ...');
            disp(' ');
            while 1
                % counter
                iter=iter+1;
                
                % ===========================================
                % loop over state space
                % ===========================================
                
                % matrix for failed points
                failedPoints=[];
                
                % transitions
                transmat=mobj.evaluateTrans(gridSt)';
                forecmat=mobj.evaluateForec(gridSt)';
                
                failedPoints_trans_T=zeros(0,size(transmat,2));
                failedPoints_trans_I=[];
                failedPoints_trans_V=zeros(0,size(VF,2)+1);
                
                %resmat=(1-damp)*resmat+damp*resmat_old;
                % a rectangular grid to speed up VF interpolations solutions   
                
                % outer loop: all exogenous states
                for ei=1:exnpt
                    tmp_grid=gr_points{ei};
                    tmp_indlog=gr_index{ei,1};
                    tmp_index=gr_index{ei,2};

                    tmp_resmat=resmat(tmp_indlog,:);
                    tmp_resmat_prev=resmat_prev(tmp_indlog,:);
                    
                    tmp_transmat=transmat(tmp_indlog,:);
                    tmp_forecmat=forecmat(tmp_indlog,:);
                    tmp_NPT=size(tmp_index,1); % nb SS pts for exog state ei
                    
                    transpts=reshape([repmat(1:exnpt,tmp_NPT,1),tmp_transmat],tmp_NPT*exnpt,NDIM);
                    Vtrans = mobj.evaluateVal(transpts)';
                    % add forecasting points
                    forecpts=reshape(tmp_forecmat,tmp_NPT*exnpt,1); % stack columns representing different states
                    Vtrans = [Vtrans, forecpts];
                    
                    % index matching for transitions
                    refidx=kron((1:tmp_NPT)',ones(1,exnpt));
                    refidx=refidx(:);                    
                    
                    
                    disp(['State ',num2str(ei)]);
                   
                    % The function solvePointList is in DSGEModel.m
                    [tmp_resmat_new,tmp_VF,tmp_TF,tmp_FF,tmp_failed]=mobj.solvePointList(tmp_grid,tmp_resmat,tmp_transmat,...
                        refidx,Vtrans,tmp_resmat_prev,printmode,[],checkTrans);  
                    
                    failedPoints=[failedPoints; tmp_index(tmp_failed)];
                    
                    
                    if revisitFailed
                        failedPoints_trans_T=[failedPoints_trans_T; tmp_transmat(tmp_failed,:)];
                        refidx_failed = ismember(refidx,find(tmp_failed));
                        
                        failedPoints_trans_I=[failedPoints_trans_I; tmp_index(refidx(refidx_failed))];
                        failedPoints_trans_V=[failedPoints_trans_V; Vtrans(refidx_failed,:)];
                    end
                    
                    % To fill the matrix where we store the results.
                    % resmat, VFnext, TFnext, FFnext.
                    resmat_prev(tmp_indlog,:)=tmp_resmat;
                    resmat(tmp_indlog,:)=tmp_resmat_new;
                    VFnext(tmp_indlog,:)=tmp_VF;
                    TFnext(tmp_indlog,:)=tmp_TF;
                    FFnext(tmp_indlog,:)=tmp_FF;
                end  
                
                if (revisitFailed && ~isempty(failedPoints))
                    disp( '~~~~~~~~~~~~~~~~~~~');
                    disp(['Revisiting failed points: ',num2str(length(failedPoints)),' add. points ...']);
                    % try to solve at failed points
                    [new_resmat,new_VF,new_TF,new_FF,n_succ]=mobj.solvePointListFailed(gridSt,failedPoints,resmat,...
                        failedPoints_trans_T,failedPoints_trans_I,failedPoints_trans_V,...
                        1,printmode,[],checkTrans);
                    resmat(failedPoints,:)=new_resmat;
                    VFnext(failedPoints,:)=new_VF;
                    TFnext(failedPoints,:)=new_TF;
                    FFnext(failedPoints,:)=new_FF;
                    disp(['Revisiting solved ',num2str(n_succ),' points.']);
				end     
                
                
                % Update functions for next iteration (but slowly!)
                mobj=mobj.updateVfct((1-damp)*VFnext+damp*VF);
                mobj=mobj.updateTfct((1-damp)*TFnext+damp*TF);
                mobj=mobj.updateFfct((1-damp)*FFnext+damp*FF);
                
                % convergence criterion (based on points in BaseGrid)
                val_range=1:5;
                VF_val = VF(:,val_range);
                VFnext_val=VFnext(:,val_range);
                [dist,wh]=max(abs(VF_val(:)-VFnext_val(:)));
                [mean_dist,col]=max(abs(mean(VF_val-VFnext_val)));
                [distT,whT]=max(abs(TF(:)-TFnext(:)));      
                [distres,whres]=max(abs(resmat_old(:)-resmat(:))); 
                [wh_1,wh_2]=ind2sub(size(VFnext_val),wh);
                [whT_1,whT_2]=ind2sub(size(TFnext),whT);
                [whres1,whres2]=ind2sub(size(resmat),whres);
                

                disp(['-- Iteration: ',num2str(iter),', max distance: ',num2str(dist),' in ',char(mobj.V_names(val_range(1)-1+wh_2)), ...
                    ' at point ',num2str(wh_1),': ',num2str(mobj.Vfct.SSGrid.Pointmat(wh_1,:))]);

                disp(['-- Iteration: ',num2str(iter),', mean distance: ',num2str(mean_dist),' in ',char(mobj.V_names(val_range(1)-1+col))]);

                disp(['-- Iteration: ',num2str(iter),', max P distance: ',num2str(distres),' in var ', num2str(whres2), ...
                    ' at point ',num2str(whres1),': ',num2str(mobj.Vfct.SSGrid.Pointmat(whres1,:))]);
                
                disp(['-- Iteration: ',num2str(iter),', max T distance: ',num2str(distT),' in col ',num2str(whT_2), ...
                    ' at point ',num2str(whT_1),': ',num2str(mobj.Vfct.SSGrid.Pointmat(whT_1,:))]);
                
                disp(' ');

                % The criteria is the second one: the mean.
                if distT<tol_avg
                    disp('Converged.');
                    break;
                elseif iter>=MAXIT
                    disp('Max.iter. exceeded.');
                    break;
                end
                
                resmat_old=resmat;
                VF=VFnext;
                TF=TFnext; 
                FF=FFnext;
            end
            
            % resulting policy functions
            mobj=mobj.updatePfct(resmat);       
        end
        
        function [simseries, varnames] = computeSimulationMoments(obj, simseries, varnames, varargin)
                       
            %NT_sim = size(simseries,1);
            
            % make HashMap with mapping of names to indices
            indexmap=java.util.HashMap;
            for i=1:length(varnames)
                if isempty(indexmap.get(varnames{i}))
                    indexmap.put(varnames{i},i);
                end
            end
            
            %--------------------------------------------------------------------------
            % transformations of raw model output
            %--------------------------------------------------------------------------
            
            % list of indices
            if obj.Params.rental
               loglist=model.HelperCollection.makeListFromNames(indexmap,{'R','P','q','pY','pM','dM','dY','cM','cY','mY','mM','vM','vY','Itilde'});
            else
                loglist=model.HelperCollection.makeListFromNames(indexmap,{'R','P','q','pY','pM','dM','dY','cM','cY','hY','mY','mM','vM','vY','Itilde'});
            end
            multlist=model.HelperCollection.makeListFromNames(indexmap,{'muI'});    

            % conversion of log-values
            simseries(:,loglist)=exp(simseries(:,loglist));
            % conversion of multipliers
            simseries(:,multlist)=max(simseries(:,multlist),0).^(1/3);


			if nargin>3
				firstrow = varargin{1};
				simseries = [firstrow; simseries];
			end
			NT_sim = size(simseries,1);
			
			if nargin>4
				params=varargin{2};
			else
				params=obj.Params;
				paramNames = fieldnames(params);
				for p = paramNames'
					p = p{:};
					tmp = params.(p);
					if ~isa( tmp, 'function_handle')
						tmp = tmp(:)';
						params.(p) = repmat( tmp, NT_sim, 1 );
					end
				end
            end            
  
            G = params.mu_G;
            
            % Function to calculate log growth rates

            Y=simseries(:,indexmap.get('Y'));
            sig_epsY=simseries(:,indexmap.get('sig_epsY'));
            sig_epsM=simseries(:,indexmap.get('sig_epsM'));
            mu_epsY=-.5*sig_epsY.^2;
            mu_epsM=-.5*sig_epsM.^2;
            %loggr = @(var) diff(log(var)) + gY;            
            
            
            % ---------------------------------------------------------------------
            % borrower debt and default: young
            % ---------------------------------------------------------------------
            P= simseries(:,indexmap.get('P'));
            pY= simseries(:,indexmap.get('pY'));
            QY= simseries(:,indexmap.get('QY')); 
%             
%             if rental
%                 hY = log(simseries(:,indexmap.get('hY'))); 
%             else
            hY = simseries(:,indexmap.get('hY')); 
            sY = simseries(:,indexmap.get('sY')); 
            %end
    
            mY = simseries(:,indexmap.get('mY')); 
            dY = simseries(:,indexmap.get('dY')); 
            BO = simseries(:,indexmap.get('BO'));
          
            
            Reb = simseries(:,indexmap.get('Reb'));
            RebMort = simseries(:,indexmap.get('RebMort'));
            nu = params.nu;
            delta_eta = params.delta_eta;
            deltaH = params.deltaH;
            lambdaY = params.lambdaY;
            lambdaM = params.lambdaM;
            piY = params.piY;
            piM = params.piM;
            xi = params.xi;
            
           
            % Change: missing BO and adding the tax.
            Ystar = Y + Reb + BO + RebMort;
            ystar = (1-nu).*Ystar;
            wYYdef = dY(1:end-1)./G(1:end-1) + (nu(2:end).*Ystar(2:end) + pY(2:end));
            % Change: eliminate BO
            wYMdef = dY(1:end-1)./G(1:end-1) + (delta_eta(2:end).*ystar(2:end))./piY(2:end);
            epsbarYY = (mY(1:end-1)./G(1:end-1) - lambdaY(2:end).*wYYdef)./((1-deltaH(2:end)).*P(2:end).*hY(1:end-1));
            epsbarYM =(mY(1:end-1)./G(1:end-1) - lambdaY(2:end).*wYMdef)./((1-deltaH(2:end)).*P(2:end).*hY(1:end-1));
            [FepsYY,~,FepsplusYY,FepsminusYY] = OLGIntermediaryModel.integEps(epsbarYY,mu_epsY(2:end),sig_epsY(2:end),[],[],[],[],[],1);            
            [FepsYM,~,FepsplusYM,FepsminusYM] = OLGIntermediaryModel.integEps(epsbarYY,mu_epsY(2:end),sig_epsY(2:end),[],[],[],[],[],1);            
            DefY = (1-piY(2:end)).*FepsYY + piY(2:end).*FepsYM;
            FepsY = (1-piY(2:end)).*(1-FepsYY) + piY(2:end).*(1-FepsYM);
            condminusY = piY(2:end).*FepsminusYM./FepsYM + (1-piY(2:end)).*FepsminusYY./FepsYY;
            condplusYY = FepsplusYY./(1-FepsYY);
            condplusYM = FepsplusYM./(1-FepsYM);
            

            wYYnodef = (1-deltaH(2:end)).*condplusYY.*P(2:end).*hY(1:end-1) - mY(1:end-1)./G(1:end-1) + dY(1:end-1)./G(1:end-1) + nu(2:end).*Y(2:end) + pY(2:end);
            % Change: eliminate BO
            wYMnodef = (1-deltaH(2:end)).*condplusYM.*P(2:end).*hY(1:end-1) - mY(1:end-1)./G(1:end-1) + dY(1:end-1)./G(1:end-1) + delta_eta(2:end).*ystar(2:end);
            FepsminusY = (1-piY(2:end)).*FepsminusYY + piY(2:end).*FepsminusYM;
            NetPayoffY = FepsY.*mY(1:end-1)./G(1:end-1) + (1-deltaH(2:end)).*(1-xi(2:end)).*FepsminusY.*P(2:end).*hY(1:end-1);  
            
            lev_mY = mY./(P.*hY);            
            lev_QY = QY./(P.*hY);
            LrateY = 1 - G(1:end-1).*NetPayoffY./mY(1:end-1);
            LGDY = 1 - G(1:end-1).*(1-deltaH(2:end)).*(1-xi(2:end)).*condminusY.*P(2:end).*hY(1:end-1)./mY(1:end-1);                       

            % HO rate
            HOrateY = sY./hY;
            
            % deposits/income
            dYrate = dY./nu;
            
            % ---------------------------------------------------------------------
            % borrower debt and default: middle
            % ---------------------------------------------------------------------
            QM= simseries(:,indexmap.get('QM')); 
            hM = simseries(:,indexmap.get('hM')); 
            mM = simseries(:,indexmap.get('mM')); 
            dM = simseries(:,indexmap.get('dM')); 
            q = simseries(:,indexmap.get('q')); 
            xI = simseries(:,indexmap.get('xI')); 
            pM= simseries(:,indexmap.get('pM'));
            
            wMdef = dM(1:end-1)./G(1:end-1) + q(2:end)+xI(2:end) + pM(2:end) + (1-delta_eta(2:end)).*ystar(2:end);
            
            epsbarM = (mM(1:end-1)./G(1:end-1) - lambdaM(2:end).*wMdef)./((1-deltaH(2:end)).*P(2:end).*hM(1:end-1));
            [FepsM,~,FepsplusM,FepsminusM] = OLGIntermediaryModel.integEps(epsbarM,mu_epsM(2:end),sig_epsM(2:end),[],[],[],[],[],1);            
            DefM = FepsM;
            condminusM = FepsminusM./FepsM;
            condplusM = FepsplusM./(1-FepsM);
            wMnodef = (1-deltaH(2:end)).*condplusM.*P(2:end).*hM(1:end-1) - mM(1:end-1)./G(1:end-1) + wMdef;
            NetPayoffM = (1-FepsM).*mM(1:end-1)./G(1:end-1) + (1-deltaH(2:end)).*(1-xi(2:end)).*P(2:end).*FepsminusM.*hM(1:end-1);
            lev_mM = mM./(P.*hM);
            lev_QM = QM./(P.*hM);            
            LrateM = 1 - G(1:end-1).*NetPayoffM./mM(1:end-1);
            LGDM = 1 - G(1:end-1).*(1-deltaH(2:end)).*(1-xi(2:end)).*condminusM.*P(2:end).*hM(1:end-1)./mM(1:end-1);             
                        
            totmdebt = mY+mM;
            totdep = dY+dM;
            Defrate = DefY.*mY(1:end-1)./totmdebt(1:end-1) + DefM.*mM(1:end-1)./totmdebt(1:end-1);
            Lrate = mY(1:end-1).*LrateY./totmdebt(1:end-1) + mM(1:end-1).*LrateM./totmdebt(1:end-1);
            totlev_Q = (QM + QY)./P;
            totlev_m = (mM + mY)./P;
            totdti = totmdebt./Ystar;
            
            yY = nu.*Ystar;
            yM = (1-nu).*Ystar;
            DTI_mM = mM./yM;
            DTI_QM = QM./yM;
            DTI_mY = mY./yY;
            DTI_QY = QY./yY; 
                        
%            WMrat = (simseries(:,indexmap.get('WM'))-pM)./ystar;
            WMrat = (simseries(:,indexmap.get('WM')))./ystar;
            
            % ratios to income
            totmdebt_Y=totmdebt./Y;
            totdep_Y = simseries(:,indexmap.get('DI'))./Y;
            
            
            % ---------------------------------------------------------------------
            % interest rates and returns
            % ---------------------------------------------------------------------

            R = simseries(:,indexmap.get('R'));
            bind_muI= simseries(:,indexmap.get('muI'))>0;
            MrateY = mY./QY - 1;
            MrateM = mM./QM - 1;
            rf = R - 1;
            % spreads fro the morgage interest rate.
            MsprY = MrateY - rf;
            MsprM = MrateM - rf;
            % average spread
            Mspr = mY./totmdebt.*MsprY + mM./totmdebt.*MsprM;
            % risky spread
            MsprYM = MrateY - MrateM;
            % tax subsidy
            QMHH= simseries(:,indexmap.get('QMHH')); 
            QYHH= simseries(:,indexmap.get('QYHH')); 
            taxsubs = mY./totmdebt.* (mY./QY - mY./QYHH) + mM./totmdebt.* (mM./QM - mM./QMHH);
            qM = QM./mM;
            tax = params.tax;
            taxsubsrate = tax.*qM.^2./(1-qM);
            
            % price rent ratios
            rhoY= simseries(:,indexmap.get('rhoY'));
            rhoM= simseries(:,indexmap.get('rhoM'));            
            PRY = P./rhoY;
            PRM = P./rhoM;
            PRavg = hY.*PRY + hM.*PRM;            
            
            expR_Bext= simseries(:,indexmap.get('expR_Bext'));
            expR_Bin = simseries(:,indexmap.get('expR_Bin'));
            expR_MM= simseries(:,indexmap.get('expR_MM'));
            expR_MY= simseries(:,indexmap.get('expR_MY'));
            
            if ~isempty(expR_Bext)
                expER_Bext=expR_Bext - R;
                expER_Bin=expR_Bin - R;
                expER_MM=expR_MM - R;
                expER_MY=expR_MY - R;
                expER_YMdiff = expER_MY - expER_MM;
                expER_MYnet = expER_MY - params.zeta;
                expER_MMnet = expER_MM - params.zeta;
            else
                expER_Bext=nan(NT_sim,1);
                expER_Bin=nan(NT_sim,1);
                expER_MM=nan(NT_sim,1);
                expER_MY=nan(NT_sim,1);
                expER_YMdiff=nan(NT_sim,1);
                expER_MYnet = nan(NT_sim,1);
                expER_MMnet = nan(NT_sim,1);
            end                

            eI = simseries(:,indexmap.get('eI'));
            eIrat_Q = eI./(QY+QM);
            eIrat_m = eI./(mY+mM);
            eIrat_Q_lag = eI(2:end)./(QY(1:end-1)+QM(1:end-1));
            eIrat_m_lag= eI(2:end)./(mY(1:end-1)+mM(1:end-1));
            payratI = xI./eI;
            min_eI = simseries(:,indexmap.get('min_eInext'));
            eqbufferI = (eI - min_eI)./eI;
            I = simseries(:,indexmap.get('I'));
            netiss = I./eI;
            netissEOP = I./(I+eI);
           
            % rate subsidy
            alpha = params.alpha;
            ratesubsY = alpha(2:end).*P(2:end).*hY(1:end-1)./QY(1:end-1);
            ratesubsM = alpha(2:end).*P(2:end).*hM(1:end-1)./QM(1:end-1);
            ratesubs = mY(1:end-1)./totmdebt(1:end-1).*ratesubsY + mM(1:end-1)./totmdebt(1:end-1).*ratesubsM;
            
            
            % ---------------------------------------------------------------------
            % Consumption and welfare
            % ---------------------------------------------------------------------                      
            
            cY=simseries(:,indexmap.get('cY'));
            cM=simseries(:,indexmap.get('cM'));
            cO=simseries(:,indexmap.get('cO'));
            loggr = @(x)(log(x(2:end)) - log(x(1:end-1)));
            cY_gr = loggr(cY);
            cM_gr = loggr(cM);
            cO_gr = loggr(cO);
            
            C = cY + cM + cO;
            consMbyY = cM./cY;
            consMObyY = (cM+cO)./cY;
            houseMbyY = hM./hY;
            equI = params.chi/2 .* I.^2;
            DWL = params.xiDWL(2:end).*(1-deltaH(2:end)).*(P(2:end).*FepsminusM.*hM(1:end-1) + FepsminusY.*P(2:end).*hY(1:end-1));
            totcost = equI(2:end)+DWL;
            Hmaint = deltaH.*P.*params.Hbar;    
            
            Ycheck = C + equI + Hmaint -  Y;
                                   
            simseries=[simseries(:,2:end),... state vars
                                      rf,MrateY,MrateM,MsprY,MsprM,Mspr,MsprYM, ... rates
                                      expER_Bext,expER_Bin,expER_MY,expER_MM,expER_YMdiff,expER_MYnet,expER_MMnet, ... returns
                                      lev_mY,lev_QY,lev_mM,lev_QM, totlev_m, totlev_Q,... borrower leverage
                                      eIrat_Q, eIrat_m, bind_muI, totmdebt, totdti, totdep, payratI, eqbufferI, ... intermediary
                                      DTI_mM, DTI_QM,DTI_mY,DTI_QY,PRY,PRM,PRavg,HOrateY,dYrate,WMrat,taxsubs,taxsubsrate,...
                                      C, equI, Hmaint, Ycheck,...
                                      totmdebt_Y,totdep_Y,netiss,netissEOP,consMbyY,consMObyY,houseMbyY];

            simseries=[simseries(2:end,:),wYYnodef, wYYdef, wYMnodef, wYMdef, ...
                                         wMnodef, wMdef, ...
                                         epsbarYY, epsbarYM, DefY, NetPayoffY, LrateY, LGDY,... % young default
                                         epsbarM, DefM, NetPayoffM, LrateM, LGDM,ratesubs,... % middle default
                                         cY_gr, cM_gr, cO_gr,... % consumption growth
                                         Defrate, Lrate, DWL, totcost,...
                                         eIrat_Q_lag, eIrat_m_lag]; 

            varnames_add={'rf','MrateY','MrateM','MsprY','MsprM','Mspr','MsprYM',...
                      'expER_Bext','expER_Bin','expER_MY','expER_MM','expER_YMdiff','expER_MYnet','expER_MMnet',...
                      'lev_mY','lev_QY','lev_mM','lev_QM','totlev_m','totlev_Q', ...
                      'eIrat_Q', 'eIrat_m','bind_muI','totmdebt','totdti','totdep','payratI', 'eqbufferI',...
                      'DTI_mM', 'DTI_QM','DTI_mY','DTI_QY','PRY','PRM','PRavg','HOrateY','dYrate','WMrat','taxsubs','taxsubsrate',...
                      'C','equI','Hmaint','Ycheck',...
                      'totmdebt_Y','totdep_Y','netiss','netissEOP','consMbyY','consMObyY','houseMbyY',...
                      'wYYnodef', 'wYYdef', 'wYMnodef', 'wYMdef', 'wMnodef', 'wMdef', ...
                      'epsbarYY', 'epsbarYM', 'DefY' ,'PayoffY', 'LrateY', 'LGDY',...
                      'epsbarM', 'DefM','PayoffM', 'LrateM', 'LGDM', 'ratesubs',... 
                      'cY_gr', 'cM_gr', 'cO_gr',...
                      'Defrate', 'Lrate', 'DWL', 'totcost',...
                      'eIrat_Q_lag','eIrat_m_lag'};
            varnames=[varnames(2:end), varnames_add];
            
        end
        
        
    end %of object methods
        
    
    %==============================================================================
    methods (Static)
        % static class-specific methods
        
        % adapative intgration for testing
        function f=integEpsFun(x,A,B,m,s,gamma)
           f=zeros(1,3);
           f(1)=(A*x+B).^(1-gamma)*lognpdf(x,m,s);
           f(2)=(A*x+B).^(-gamma)*lognpdf(x,m,s);
           f(3)=x*(A*x+B).^(-gamma)*lognpdf(x,m,s);
        end
        
        % Gauss-Legendre for speed
        function f=integEpsGL(epsbar,A,B,m,s,gamma,nodes,weights)
           f=zeros(1,3);
           up=12*s;
           nd=nodes*(up-epsbar)/2 + (up+epsbar)/2;
           lnpdfv=lognpdf(nd,m,s);
           f(1)= sum(weights.* (A*nd+B).^(1-gamma) .* lnpdfv)*(up-epsbar)/2;
           f(2)= sum(weights.* (A*nd+B).^(-gamma) .* lnpdfv)*(up-epsbar)/2;
           f(3)= sum(weights.* nd.*(A*nd+B).^(-gamma) .* lnpdfv)*(up-epsbar)/2;            
        end
        
        function [Feps,feps,Fepsplus,Fepsminus,RepsVF,RepsFOC1,RepsFOC2] = integEps(epsbar,mueps,sigeps,gamma,Avec,Bvec,nodes,weights,simmode)            
            % default rates and truncated expectations
            epsbar(epsbar<0)=0;
            Feps=logncdf(epsbar,mueps,sigeps);
            feps=lognpdf(epsbar,mueps,sigeps);
            sigeps2 = sigeps.^2;
            Fepsminus = exp(mueps + sigeps2/2) .* normcdf((log(epsbar) - mueps - sigeps2)./sigeps);
            Fepsplus = exp(mueps + sigeps2/2) .* normcdf((mueps + sigeps2 - log(epsbar))./sigeps);
            % integrals for SDFs
            nst=length(epsbar);
            RepsVF=zeros(nst,1);
            RepsFOC1=zeros(nst,1);
            RepsFOC2=zeros(nst,1);
            if ~simmode
                for n=1:nst
                                     % fun=@(x)OLGIntermediaryModel.integEpsFun(x,Avec(n),Bvec(n),mueps(n),sigeps(n),gamma);
                                    % int=integral(fun,epsbar(n),Inf,'ArrayValued',true,'AbsTol',1e-3,'RelTol',1e-3);
                    int=OLGIntermediaryModel.integEpsGL(epsbar(n),Avec(n),Bvec(n),mueps(n),sigeps(n),gamma,nodes,weights);
                    RepsVF(n)=int(1);
                    RepsFOC1(n)=int(2);
                    RepsFOC2(n)=int(3);
                end
            end
        end
         
 

%% Steady state


        function [fx,stvals]=compStSt(sol,params,print,nodes,weights)
          % unpack relevant parameters
          beta=params.beta;
          betaO=params.betaO;
          gamma=params.gamma;
          theta=params.theta;
          psiY=params.psiY;
          psiM=params.psiM;
          chi=params.chi;
          nu=params.nu;
          tau=params.tau;
          ebarY=params.ebarY;
          ebarM=params.ebarM;
          lambdaY=params.lambdaY;
          lambdaM=params.lambdaM;
          sig_epsY=params.sig_epsY(1);
          mu_epsY=params.mu_epsY(1);
          sig_epsM=params.sig_epsM(1);
          mu_epsM=params.mu_epsM(1);
          xi=params.xi;
          xiDWL=params.xiDWL;
          Hbar=params.Hbar;
          G=params.mu_G;
          g=params.mu_g;
          piM=params.piM;
          piY=params.piY;
          deltaH=params.deltaH;
          phi=params.phi;
          delta_eta=params.delta_eta;
          zeta=params.zeta;
          alpha=params.alpha;
          tax=params.tax;
          % Foreign demand for domestic assets
          % Hard constraints
          hardconstr_LTV = params.LTV;
          hardconstr_PTI = params.PTI;
          % Rental market
          rental = params.rental;
            
          if rental
            kappaY = params.kappaY;
          end
          

          % solution vars
          P=exp(sol(1));
          r=sol(2);
          q=exp(sol(3));
          pY=exp(sol(4));
          pM=exp(sol(5));
          cY=exp(sol(6));
          if rental
          hY=sol(7);
          else 
          hY=exp(sol(7));
          end
          mY=exp(sol(8));
          cM=exp(sol(9));
          dM=exp(sol(10));
          dY=exp(sol(11));
          mM=exp(sol(12));
          xI=sol(13);
          BO=exp(sol(14));
          vM=exp(sol(15));
          vY=exp(sol(16));
          Reb=exp(sol(17));
          
          % Multiplier of the Intermediary's capital req. constraint.
          muI=sol(18);
          muIplus=max(0,muI)^3;
          muIminus=max(0,-muI)^3;
          
          % Multiplier of the non-neg. constraint of housing for the young.
          if rental
          muI_rentY=sol(19);
          muIplus_rentY=max(0,muI_rentY)^3;
          muIminus_rentY=max(0,-muI_rentY)^3;
          end
          
          
          %sY = exp(sol(20));

          %%%%%%%%%%%%%%%%%%%
          % Market clearing %
          %%%%%%%%%%%%%%%%%%%
          % Market clearing for the deposits.
          DI = dY + dM;
          % Market clearing for the housing capital.
          hM = Hbar - hY;
          % Market clearing for the rental market, using implictely the 
          % solution from the HH problem of both generations.
          if rental
          sY = (Hbar)/(1+(cM/cY));
          sM = Hbar - sY;
          else
          % Market clearing for the inidividual rental markets
          sM = hM;
          sY = hY;    
          end
          
         
          %%%%%%%%%%%%%%%%%%%%%%
          % Young HH Equations %
          %%%%%%%%%%%%%%%%%%%%%% 
          
          % Aggregate Income:
          % Adding bequests and Rebates
          Ystar = 1 + Reb + BO;
          ystar = (1-nu)*Ystar;
          wdefYY = nu*Ystar + pY + dY/G; 
          epsbarYY = (mY/G - lambdaY*wdefYY)/((1-deltaH)*P*hY);
          
          [FepsYY,fepsYY,FepsplusYY,FepsminusYY,...
              RepsVFYY,RepsFOC1YY,RepsFOC2YY] = OLGIntermediaryModel.integEps(epsbarYY,mu_epsY,sig_epsY,gamma,(1-deltaH)*P*hY,wdefYY-mY/G,nodes,weights,0);
          wdefYM = (delta_eta*ystar)/piY + dY/G;
          wdefYM_agg = delta_eta*ystar + piY*dY/G;
          epsbarYM = (mY/G - lambdaY*wdefYM)/((1-deltaH)*P*hY);
          
          [FepsYM,fepsYM,FepsplusYM,FepsminusYM,...
              RepsVFYM,RepsFOC1YM,RepsFOC2YM] = OLGIntermediaryModel.integEps(epsbarYM,mu_epsM,sig_epsM,gamma,(1-deltaH)*P*hY,wdefYM-mY/G,nodes,weights,0);
          
          DefY = (1-piY)*FepsYY + piY*FepsYM;
          FepsminusY = (1-piY)*FepsminusYY + piY*FepsminusYM;
          RebateY = (1-piY)*FepsYY*lambdaY*wdefYY + FepsYM*lambdaY*wdefYM_agg + (xi-xiDWL)*FepsminusY*(1-deltaH)*P*hY;     
          wnodefYY_agg = FepsplusYY*(1-deltaH)*P*hY + (1-FepsYY)*(nu*Ystar + pY - mY/G + dY/G);
          
          % Total wealth of young agents:
          WY = piY*(nu*Ystar + pY) + (1-piY)*(wnodefYY_agg + FepsYY*(1-lambdaY)*wdefYY);  
          % Total wealth of young agents tha age: (written as in the paper)
          WYM = piY*(FepsplusYM*(1-deltaH)*P*hY - (1-FepsYM)*mY/G) + wdefYM_agg*(1-lambdaY*FepsYM);       

          rhoY = theta*cY/((1-theta)*sY);
          QYHH = cY + rhoY*sY + (P-rhoY)*hY + dY/(1+r) + pY - WY;                    
          QY = QYHH/(1+tax);
         
          %%%%%%%%%%%%%%%%%%%%%%%%%%%%
          % Middle Aged HH Equations %
          %%%%%%%%%%%%%%%%%%%%%%%%%%%%

          wdefM = q+xI + dM/G + pM + (1-delta_eta)*ystar;
          epsbarM = (mM/G - lambdaM*wdefM)/((1-deltaH)*P*hM);

          [FepsM,fepsM,FepsplusM,FepsminusM,...
              RepsVFM,RepsFOC1M,RepsFOC2M] = OLGIntermediaryModel.integEps(epsbarM,mu_epsM,sig_epsM,gamma,(1-deltaH)*P*hM,wdefM-mM/G,nodes,weights,0);
          
          wdefM_agg = (q+xI) + dM/G + (pM + (1-delta_eta)*ystar);

          RebateM = (xi-xiDWL)*FepsminusM*(1-deltaH)*P*hM + FepsM*lambdaM*wdefM_agg;
  
          wnodefM_agg = FepsplusM*(1-deltaH)*P*hM + (1-FepsM)*(q+xI + pM + (1-delta_eta)*ystar +  dM/G - mM/G);
          
          % Total wealth of middle and old agents 
          WMO = WYM + wnodefM_agg + FepsM*(1-lambdaM)*wdefM_agg;
          WO = piM*WMO;
          WM = WMO - WO;
          
          % Now, both generations face the same price
          if rental
          rhoM = rhoY;
          else
          % With separate rental markets
          rhoM = theta*cM/((1-theta)*sM); 
          end
          
          QMHH = cM + rhoM*sM + (P-rhoM)*hM + q + pM + (dM)/(1+r) - WM;
          QM = QMHH/(1+tax);

          cO = WO/(1 + phi^(1/gamma));
          Phi = (1/(1+ phi^(1/gamma)))^(1-gamma) + phi*(phi^(1/gamma)/(1+ phi^(1/gamma)))^(1-gamma) * betaO;
          
          %%%%%%%%%%%%%%%%%%%%%%%%%%%
          % Intermediary  Equations %
          %%%%%%%%%%%%%%%%%%%%%%%%%%%
          
          PayoffY = ((1-piY)*(1-FepsYY) + piY*(1-FepsYM))*mY/G  ...
                       + P*hY*(1-deltaH)*( alpha + ((1-piY)*FepsminusYY + piY*FepsminusYM)*(1-xi) );
          PayoffM = (1-FepsM)*mM/G + P*hM*(1-deltaH)*( alpha + FepsminusM*(1-xi) );
          NetPayoffY = PayoffY - P*hY*(1-deltaH)*alpha;
          NetPayoffM = PayoffM - P*hM*(1-deltaH)*alpha;
          
          eI = PayoffY + PayoffM - DI/G;
          % Using the defintion of the intermediary's equity (tau*eI-I).
          RebMort_today = zeta*(QY+QM) - (QMHH-QM) - (QYHH-QY);
          RebateI = RebMort_today - (1-deltaH)*alpha*P*Hbar;
          I = tau*eI - xI;
                    
          %%%%%%%%%%%%%%%%%%%%%%%%%%%%
          % Middle Age agents part 2 % 
          %%%%%%%%%%%%%%%%%%%%%%%%%%%%
          
          SM = (P-rhoM)*hM + q + pM + (dM)/(1+r) - QMHH;  % implied savings
          R_h = G*(1-deltaH)*P/(P-rhoM);  % equilibrium returns
          R_e = G*(q+xI)/q;
          dMhat = (dM)/SM;  
          RepsVFM = RepsVFM*(G/SM)^(1-gamma);
          RMdef = G*(1-lambdaM)*wdefM/SM;
          AZ = psiM*(dMhat)^(1-gamma) + beta*(piM*Phi + (1-piM)*vM)*(RepsVFM + FepsM*RMdef^(1-gamma));
          ThetaZ = ((1-theta)^(1-theta) * (theta/rhoM)^theta )^(1-gamma);
          SMshare = AZ^(1/gamma) / (AZ^(1/gamma) + ThetaZ^(1/gamma));  
          vMcheck = AZ * ( AZ^(1/gamma) / (AZ^(1/gamma) + ThetaZ^(1/gamma)) )^(1-gamma) ...
              + ThetaZ * ( ThetaZ^(1/gamma) / (AZ^(1/gamma) + ThetaZ^(1/gamma)) )^(1-gamma);
          
          %%%%%%%%%%%%%%%%%%%%%%%%%%%
          % Middle Age agents SDF's % 
          %%%%%%%%%%%%%%%%%%%%%%%%%%%
          
          RepsFOC1M = RepsFOC1M*(G/SM)^(-gamma);
          RepsFOC2M = RepsFOC2M*(G/SM)^(-gamma);
          

          MMnodef_m = beta * (SMshare)^(-gamma)*RepsFOC1M * (piM*Phi/vM + 1-piM);
          MMnodef_h = beta * (SMshare)^(-gamma)*RepsFOC2M * (piM*Phi/vM + 1-piM);
          MMdef = beta * (SMshare*RMdef)^(-gamma) * (piM*Phi/vM + 1-piM);
          
          SDFM = MMnodef_m + FepsM*(1-lambdaM)*MMdef;
          MI = SDFM * (tau*(1-chi*I) + (1-tau));
          
          %%%%%%%%%%%%%%%%%%%%%%
          % Derivatives of q_M %
          %%%%%%%%%%%%%%%%%%%%%%

          epshatM = mM/G - (1-xi)*(1-deltaH)*P*hM*epsbarM;
          QM_m = (1+tax)*(muIplus*(1-ebarM)  + MI)/(1+zeta)*(fepsM*epshatM/(G*(1-deltaH)*P*hM) + (alpha + (1-xi)*FepsminusM)*(1-deltaH)*P*hM/mM );
          QM_h = (1+tax)*(muIplus*(1-ebarM)  + MI)/(1+zeta)*G*(1-deltaH)*P*(alpha + (1-xi)*FepsminusM+ fepsM*epshatM*epsbarM/(G*(1-deltaH)*P*hM));
          QM_d = (1+tax)*(muIplus*(1-ebarM)  + MI)/(1+zeta)*fepsM*epshatM*lambdaM/(G*(1-deltaH)*P*hM);
          QM_b = (1+tax)*(muIplus*(1-ebarM)  + MI)/(1+zeta)*G*fepsM*epshatM*(q+xI)*lambdaM/(G*(1-deltaH)*P*hM);
          QM_eta = (1+tax)*(muIplus*(1-ebarM)  + MI)/(1+zeta)*G*fepsM*epshatM*((1-delta_eta)*ystar+pM)*lambdaM/(G*(1-deltaH)*P*hM);
          
          %%%%%%%%%%%%%%%%%%%%%%%%%%%
          % Young Age agents part 2 % 
          %%%%%%%%%%%%%%%%%%%%%%%%%%%% 

          SY = (P-rhoY)*hY + pY + dY/(1+r) - QYHH;  % implied savings
          R_hY = G*(1-deltaH)*P/(P-rhoY);  % equilibrium returns
          dYhat = dY/SY;
          RYYdef = G*(1-lambdaY)*wdefYY/SY;
          RepsVFYY = RepsVFYY*(G/SY)^(1-gamma);
          RYMdef = G*(1-lambdaY)*wdefYM/SY;
          RepsVFYM = RepsVFYM*(G/SY)^(1-gamma);
          AY = psiY*(dYhat)^(1-gamma) ...
                + beta*(1-piY)*vY*(RepsVFYY + FepsYY*RYYdef^(1-gamma)) ...
                + beta*piY*vM*(RepsVFYM + FepsYM*RYMdef^(1-gamma));
          ThetaY = ( (1-theta)^(1-theta) * (theta/rhoY)^theta )^(1-gamma);
          SYshare = AY^(1/gamma) / (AY^(1/gamma) + ThetaY^(1/gamma));
          vYcheck = AY * ( AY^(1/gamma) / (AY^(1/gamma) + ThetaY^(1/gamma)) )^(1-gamma) ...
              + ThetaY * ( ThetaY^(1/gamma) / (AY^(1/gamma) + ThetaY^(1/gamma)) )^(1-gamma);
          
          
          %%%%%%%%%%%%%%%%%%%%%%%%%%
          % Young Age agents SDF's % 
          %%%%%%%%%%%%%%%%%%%%%%%%%%       
          

          RepsFOC1YY = RepsFOC1YY*(G/SY)^(-gamma);
          RepsFOC2YY = RepsFOC2YY*(G/SY)^(-gamma);
          RepsFOC1YM = RepsFOC1YM*(G/SY)^(-gamma);
          RepsFOC2YM = RepsFOC2YM*(G/SY)^(-gamma);
          
          MYnodef_m = beta * ((1-piY)*(SYshare)^(-gamma)*RepsFOC1YY + piY*vM/vY*(SYshare)^(-gamma)*RepsFOC1YM );
          MYnodef_h = beta * ((1-piY)*(SYshare)^(-gamma)*RepsFOC2YY + piY*vM/vY*(SYshare)^(-gamma)*RepsFOC2YM );
          MYdef = beta * ((1-piY)*FepsYY*(SYshare*RYYdef)^(-gamma) + piY*FepsYM*vM/vY*(SYshare*RYMdef)^(-gamma) );
          SDFY = MYnodef_m + (1-lambdaY)*MYdef;
          MYY = beta * (1-piY)*(RepsFOC1YY*(SYshare)^(-gamma) + FepsYY*(1-lambdaY)*(SYshare*RYYdef)^(-gamma));
          MYM = beta * piY*vM/vY*(RepsFOC1YM*(SYshare)^(-gamma) + FepsYM*(1-lambdaY)*(SYshare*RYMdef)^(-gamma));
          
          %%%%%%%%%%%%%%%%%%%%%%
          % Derivatives of q_Y %
          %%%%%%%%%%%%%%%%%%%%%%
          
          epshatYY = mY/G - (1-xi)*(1-deltaH)*P*hY*epsbarYY;
          epshatYM = mY/G - (1-xi)*(1-deltaH)*P*hY*epsbarYM;
          FhatY = piY*(1-FepsYM) + (1-piY)*(1-FepsYY);
          fhatY = (piY*fepsYM*epshatYM + (1-piY)*fepsYY*epshatYY)/(G*(1-deltaH)*P*hY);
          QY_m = (1+tax)*(muIplus*(1-ebarY)  + MI)/(1+zeta)*(fhatY + (alpha + FepsminusY*(1-xi))*(1-deltaH)*P*hY/mY);
          QY_d = (1+tax)*(muIplus*(1-ebarY)  + MI)/(1+zeta)*lambdaY*fhatY;
          fhhY = (piY*fepsYM*epshatYM*epsbarYM + (1-piY)*fepsYY*epshatYY*epsbarYY)/(G*(1-deltaH)*P*hY);
          QY_h = (1+tax)*(muIplus*(1-ebarY)  + MI)/(1+zeta)*G*(1-deltaH)*P*(alpha + FepsminusY*(1-xi)+ fhhY);
          fhhhY = (piY*fepsYM*epshatYM*((delta_eta*ystar)/piY) + (1-piY)*fepsYY*epshatYY*(nu*Ystar + pY))/(G*(1-deltaH)*P*hY);
          QY_eta = (1+tax)*(muIplus*(1-ebarY)  + MI)/(1+zeta)*G*lambdaY*fhhhY;

          %%%%%%%%%%%%%%%%%%%%%%%%%%
          % First-order conditions %
          %%%%%%%%%%%%%%%%%%%%%%%%%%
          
          % Middle-aged %
         
          % Effective returns to housing
          Rh_eff = G*(1-deltaH)*P/(P - rhoM - QM_h);
          % Effective returns to mortgages
          Rm_eff = 1/(QMHH/mM - QM_m);
          % Effective returns to deposits
          Rd_eff = (1+r)/(1-(1+r)*QM_d);
          % Effective returns to bank the specific generation asset: equity 
          Re_eff = G*(q+xI)/(q - QM_b);
          % Effective returns of the specific generation asset for the middle age
          % agents for endowments.
          RMMeta_eff = G*((1-delta_eta)*ystar + pM)/(pM - QM_eta);
          cyieldM = psiM*(dMhat*SMshare)^(-gamma)/vM;
                        
          
          % FOC: no excess return between deposits and equity (EQ 32)  
          fx(1)= cyieldM* Rd_eff + SDFM*(Rd_eff - Re_eff);
          % FOC: no excess return between the equity and the housing
          fx(2)= SDFM*Re_eff - MMnodef_h*Rh_eff; % PF choice forbank equity
          % Equation (33), this combines de FOC of housing and mortgages.
          fx(3)= MMnodef_h*Rh_eff - MMnodef_m*Rm_eff; % housing
          % This euqation relates both specific assets for the middle aged
          % generation: equity and the endowment asset.
          fx(4)= Re_eff - RMMeta_eff; 
          
          
          fx(5)= SM - SMshare*WM;
          fx(6)= vM - vMcheck; % recursion for VF
          
          % Young-age %
          
          % Effective returns of housing.
          RYh_eff = G*(1-deltaH)*P/(P - rhoY - QY_h);
          % Effective returns of mortgages.
          RYm_eff = 1/(QYHH/mY - QY_m);
          % Effective returns of deposits.                    
          RYd_eff = (1+r)/(1-(1+r)*QY_d);
          % Effective returns of the specific generation asset when you
          % don't age.
          RYYeta_eff = G*(nu*Ystar + pY)/(pY - QY_eta);
          % Effective returns of the specific generation asset when you
          % age.          
          RYMeta_eff = G*((delta_eta*ystar)/piY)/(pY - QY_eta);
          cyieldY = psiY*(dYhat*SYshare)^(-gamma)/vY;
          
          % Equation (42), this combines de FOC of deposits and mortgages.
          fx(7)= cyieldY*RYd_eff + SDFY*RYd_eff - MYnodef_m*RYm_eff;
          % Equation (43), this combines de FOC of housing and mortgages.
          if rental
          fx(8)= MYnodef_m*RYm_eff - MYnodef_h*RYh_eff - (muIplus_rentY+(kappaY*(hY)^(-gamma))) * (SYshare^(-gamma)/vY)*(1/(P - rhoY - QY_h));
          else
          fx(8)= MYnodef_m*RYm_eff - MYnodef_h*RYh_eff;
          end
          
          % Equation (44), this combines de FOC of housing andt the asset
          % that is specific to generation y.
          if rental
          fx(9)= MYnodef_h*RYh_eff - (MYY*RYYeta_eff + MYM*RYMeta_eff) + (muIplus_rentY+(kappaY*(hY)^(-gamma))) * (SYshare^(-gamma)/vY)*(1/(P - rhoY - QY_h));
          else
          fx(9)= MYnodef_h*RYh_eff - (MYY*RYYeta_eff + MYM*RYMeta_eff);
          end
          
          fx(10)= SY - SYshare*WY;
          fx(11)= vY - vYcheck;
          
          % Intermediaries %
          % Equation (34)
          fx(12)= (1+zeta)*QY - (muIplus*(1-ebarY) + MI)*PayoffY*G;
          fx(13)= (1+zeta)*QM - (muIplus*(1-ebarM)  + MI)*PayoffM*G;
          fx(14)= 1/(1+r) - muIplus  - MI;
          fx(15)= (1-ebarY)*PayoffY*G + (1-ebarM)*PayoffM*G - DI - muIminus;
          fx(16)= (1+zeta)*(QY+QM) - ((1-tau)*eI + I - chi*I^2/2 + DI/(1+r));
          
          
          % Market clearing equation %
          % Budget of the old agent
          fx(17)= BO - (WO-cO);
          % That the rebates sum to the total one used in the defition of
          % Y_star
          fx(18)= Reb - RebateM - RebateY - RebateI;
          
          % Non-negativity constraint housing for the young
          if rental
          fx(19) = hY -muIminus_rentY;
          end
           
          %%%%%%%%%%%%%%%%%%%
          % Extra equations %
          %%%%%%%%%%%%%%%%%%%
          
          % For the Rebate:
          Rebate_today = Reb;
          Y_star_today = Ystar;

          % check that goods market adds up
          C = cY + cM + cO;
          YDWL = 1 - xiDWL*(FepsminusM*(1-deltaH)*P*hM + FepsminusY*(1-deltaH)*P*hY);
          Ycheck = C + chi*I^2/2 + deltaH*P*Hbar  -  YDWL ;
          Wcheck = WY + WMO + eI - (q+xI) - pY - pM - (1-deltaH)*P*Hbar - BO - RebMort_today;
          
          % returns etc
          R_mY = G*PayoffY/QY;
          R_mM = G*PayoffM/QM;
          R_ine = (G*(PayoffY+PayoffM) - DI)/(QM+QY - DI/(1+r));
          lrateY = 1 - G*NetPayoffY/mY;
          lrateM = 1 - G*NetPayoffM/mM;

          % wealth states for dynamic code
          whatM = WMO-(q+xI)-pM;
          whatY = WY-pY;
          
          % Endowment yields for the PTI constraints
          end_yieldM = (1-delta_eta)*(1-nu)*Ystar;
          end_yieldY = ((1-piY)*nu+piY*((delta_eta)/piY)*(1-nu))*Ystar;
          
          % mortgage rates
          mrateY = mY/(QY)-1;
          mrateM = mM/(QM)-1;
          
          if print
              % print steady state values
              disp(' ');
              disp('Analytic steady state');
              disp('--- Prices ---');
              disp(['r: ',num2str(r)]);
              disp(['P: ',num2str(P)]);
              disp(['pY: ',num2str(pY)]);
              disp(['pM: ',num2str(pM)]);
              disp(['MrateY: ',num2str(mrateY)]);
              disp(['MrateM: ',num2str(mrateM)]);
              disp(['rhoY ',num2str(rhoY)]);
              disp(['rhoM: ',num2str(rhoM)]);
              disp(['cyieldM: ',num2str(cyieldM)]);
              disp(['cyieldY: ',num2str(cyieldY)]);
              disp('--- Housing ---');
              disp(['hY: ',num2str(hY)]);
              disp(['hY p.c.: ',num2str(piY*hY)]);
              disp(['hM: ',num2str(hM)]);
              disp(['hM p.c.: ',num2str(piM*hM)]);
              disp('--- Returns ---');
              disp(['R_e: ',num2str(R_e)]);
              disp(['R_mM: ',num2str(R_mM)]);
              disp(['R_mY: ',num2str(R_mY)]);
              disp(['R_hM: ',num2str(R_h)]);
              disp(['R_hY: ',num2str(R_hY)]);
              disp(['R_ine: ',num2str(R_ine)]);
              disp('--- Debt ---');
              disp(['mY: ',num2str(mY)]);
              disp(['mM: ',num2str(mM)]);
              disp(['LTVY: ',num2str(mY/(P*hY))]);
              disp(['LTVQY: ',num2str(QY/(P*hY))]);
              disp(['LTVM: ',num2str(mM/(P*hM))]);
              disp(['LTVQM: ',num2str(QM/(P*hM))]);
              disp(['dY: ',num2str(dY)]);
              disp(['dY/inc: ',num2str(dY/nu)]);
              disp(['dM: ',num2str(dM)]);
              disp(['epsbarYY: ',num2str(epsbarYY)]);
              disp(['epsbarYM: ',num2str(epsbarYM)]);
              disp(['DefY: ',num2str(DefY)]);
              disp(['epsbarM: ',num2str(epsbarM)]);
              disp(['DefM: ',num2str(FepsM)]);
              disp(['LrateM: ',num2str(lrateM)]);
              disp(['LrateY: ',num2str(lrateY)]);              
              disp('--- Intermediary ---');
              disp(['eI: ',num2str(eI)]);
              disp(['eIratm: ',num2str(eI/(mM+mY))]);
              disp(['I: ',num2str(I)]);
              disp(['muI: ',num2str(muIplus)]);
              disp(['xI: ',num2str(xI)]);
              disp(['xIrat: ',num2str(xI/eI)]);
              disp(['q: ',num2str(q)]);
              disp('--- Wealth and consumption ---');
              disp(['YDWL: ',num2str(YDWL)]);
              disp(['Rebate: ',num2str(Reb)]);
              disp(['cY/C: ',num2str(cY/C)]);
              disp(['cM/C: ',num2str(cM/C)]);
              disp(['cO/C: ',num2str(cO/C)]);
              disp(['wO: ',num2str(WO)]);
              disp(['BO: ',num2str(BO)]);
              disp(['WM: ',num2str(WM)]);
              disp(['WO/WM: ',num2str(WO/WM)]);
              disp(['WTIM: ',num2str((WM)/(1-nu))]);
              disp(['whatM: ',num2str(whatM)]);
              disp(['vM: ',num2str(vM)]);
              disp(['vY: ',num2str(vY)]);
              disp(['WY: ',num2str(WY)]);
              disp(['WTIY: ',num2str((whatY+nu)/nu)]);
              disp(['whatY: ',num2str(whatY)]);
              disp('--- Endowment Yields ---');
              disp(['end_yieldM: ',num2str(end_yieldM)]);
              disp(['end_yieldY: ',num2str(end_yieldY)]);
              disp('--- Checks ---');
              disp(['Ycheck: ',num2str(Ycheck)]);
%             disp(['SDFMcheck: ',num2str(MMnodef*FepsM-1/Rm_eff)]);
%              disp(['SDFMcheck: ',num2str(SDFM-1/RMMeta_eff)]);
              disp(['SDFMcheck: ',num2str(SDFM-1/Re_eff)]);
              %disp(['SDFYcheck: ',num2str(((kappaY/hY)) * (SYshare^(-gamma)/vY)*(1/(P - rhoY - QY_h)) + MYnodef_h * RYh_eff-1)]);
              disp(['Wcheck: ',num2str(Wcheck)]);
              disp('=============================================');
              disp('Calibration')
              disp('=============================================');
              disp(['r: ',num2str(r)]);
              disp(['orig. spr.: ',num2str(mrateM-r)]);
              disp(['Y-M spr.: ',num2str(mrateY-mrateM)]);
              disp(['r: ',num2str(r)]);  
              disp(['Net rebate: ',num2str(RebateI)]);                            
              disp('----------------------------------');
              disp(['H/Y: ',num2str(P)]);
              disp(['dY/inc: ',num2str(dY/nu)]);
              disp(['WM/Y: ',num2str(WM-2/3*pM)]);
              disp(['WTIM: ',num2str((WM-2/3*pM)/(1-nu))]);
              disp(['LTVY: ',num2str(mY/(P*hY))]);
              disp(['LTVM: ',num2str(mM/(P*hM))]);
              disp('----------------------------------');
              disp(['Def.agg: ',num2str(mY/(mY+mM)*DefY + mM/(mY+mM)*FepsM)]);
              disp(['Lrate.agg: ',num2str(mY/(mY+mM)*lrateY + mM/(mY+mM)*lrateM)]);
              disp('----------------------------------');
              disp(['cY/C: ',num2str(cY/C)]);
              disp(['cO/C: ',num2str(cO/C)]);
              disp('----------------------------------');
              disp(['cyieldM: ',num2str(cyieldM)]);        
              disp('----------------------------------');
              disp(['xI rate: ',num2str(xI/eI)]);
              disp('----------------------------------');              
              disp(['H.O. rate: ',num2str(hY/sY)]);
          end         

          if rental
          Sol=struct('R',1+r,...
              'P',P,...
              'q',q,...
              'pY',pY,...
              'pM',pM,...
              'dM',dM,...
              'dY',dY,...
              'cM',cM,...
              'cY',cY,...
              'hY',hY,...
              'mY',mY,...
              'mM',mM,...
              'Itilde',1-chi*I,...
              'vM',vM,...
              'vY',vY,...
              'muI',muIplus,...
              'RebMort',RebMort_today,...
              'muI_rentY',muI_rentY);
          else
          Sol=struct('R',1+r,...
              'P',P,...
              'q',q,...
              'pY',pY,...
              'pM',pM,...
              'dM',dM,...
              'dY',dY,...
              'cM',cM,...
              'cY',cY,...
              'hY',hY,...
              'mY',mY,...
              'mM',mM,...
              'Itilde',1-chi*I,...
              'vM',vM,...
              'vY',vY,...
              'muI',muIplus,...
              'RebMort',RebMort_today);
         end
                
            if hardconstr_LTV || hardconstr_PTI
                if  hardconstr_LTV && hardconstr_PTI
                    % LTV
                    Sol.lamM_LTV=-0.5;
                    Sol.lamY_LTV=-0.5;
                    Sol.lamM_PTI=-0.5;
                    Sol.lamY_PTI=-0.5;

                elseif hardconstr_LTV && not(hardconstr_PTI)
                    % LTV
                    Sol.lamM_LTV=-0.5;
                    Sol.lamY_LTV=-0.5;
                
                elseif not(hardconstr_LTV) && hardconstr_PTI
                    % LTV
                    Sol.lamM_PTI=-0.5;
                    Sol.lamY_PTI=-0.5;
                end
            
            end
          
          
          V=struct('q',q,...
              'I',I,...
              'vM',vM,...
              'vY',vY,...
              'P',P,...
              'pY',pY,...
              'pM',pM,...
              'BO',BO,...
              'R',1+r,...
              'RebMort',RebMort_today);

          
          Add=struct('WY',WY,...
                'WM',WM,...
                'DI',DI,...
                'hM',hM,...
                'sM',sM,...
                'rhoM',rhoM,...
                'sY',sY,...
                'rhoY',rhoY,...
                'QYHH',QYHH,...
                'QMHH',QMHH,...
                'QY',QY,...
                'QM',QM,...
                'I',I,...
                'xI',xI,...
                'SM', SM, ...
                'SY',SY,...
                'cO',cO,...
                'BO',BO,...
                'gY',g,...
                'end_yieldY',end_yieldY,...
                'end_yieldM',end_yieldM,...
                'Rebate_today',Rebate_today,...
                'Y_star_today',Y_star_today);
            
          
          
          State=struct('whatM',whatM,...
              'eI',eI,...
              'Y',1);
          
          
          stvals=struct('Sol',Sol,...
              'V',V,...
              'Add',Add,...
              'State',State,...
              'Rebate',Reb);
          
        end
                      
        
        function [solguessvec,Vguessvec,V_names]=assignGuess(stv,hardconstr_LTV,hardconstr_PTI,rental)

           solguess=struct('R',log(stv.Sol.R),...
                'P',log(stv.Sol.P),...
                'q',log(stv.Sol.q),...
                'pY',log(stv.Sol.pY),...
                'pM',log(stv.Sol.pM),...
                'dM',log(stv.Sol.dM),...
                'dY',log(stv.Sol.dY),...
                'cM',log(stv.Sol.cM),...
                'cY',log(stv.Sol.cY),...
                'hY',log(stv.Sol.hY),...
                'mY',log(stv.Sol.mY),...
                'mM',log(stv.Sol.mM),...
                'Itilde',log(stv.Sol.Itilde),...
                'vM',log(stv.Sol.vM),...
                'vY',log(stv.Sol.vY),...
                'muI',0.5,...
                'RebMort',stv.Sol.RebMort);
                
            if rental
                solguess.muI_rentY=1;
                solguess.hY=stv.Sol.hY;
            end
            
            if hardconstr_LTV || hardconstr_PTI
                if  hardconstr_LTV && hardconstr_PTI
                    % LTV
                    solguess.lamM_LTV=-0.5;
                    solguess.lamY_LTV=-0.5;
                    solguess.lamM_PTI=-0.5;
                    solguess.lamY_PTI=-0.5;

                elseif hardconstr_LTV && not(hardconstr_PTI)
                    % LTV
                    solguess.lamM_LTV=-0.5;
                    solguess.lamY_LTV=-0.5;
                
                elseif not(hardconstr_LTV) && hardconstr_PTI
                    % LTV
                    solguess.lamM_PTI=-0.5;
                    solguess.lamY_PTI=-0.5;
                end
            
            end

            solguessvec=model.DSGEModel.structToVec(solguess);
                     
            Vguess=struct('q',stv.V.q,...
                'I',stv.V.I,...
                'vM',stv.V.vM,...
                'vY',stv.V.vY,...
                'P',stv.V.P,...
                'pY',stv.V.pY,...
                'pM',stv.V.pM,...
                'BO',stv.V.BO,...
                'R',stv.V.R,...
                'RebI',stv.V.RebMort);
            
            Vguessvec=model.DSGEModel.structToVec(Vguess);
            V_names=fieldnames(Vguess);
            
        end
                
        
        function [x,fx,exit,i]=tryOtherGuesses(fhand,gvec,options)
            % list of guesses
            
            %  Toggle Lagrange multipliers as follows:
            %
            %                              
            %   lamB   +   -   -   +
            %   muR    -   -   -   - 
            %   lamR   -   +   -   +
            %   lamS   -   -   -   - 
           

            gindex={16,16};
            gvals={-0.5,0.5};            
            
            for i=1:length(gindex)
                newguess=gvec;
                newguess(gindex{i})=gvals{i};
                [x,fx,exit]=fsolve(fhand,newguess,options);
%                if exit>=0
                if exit>0
                    break;
                end                                
            end
        end
                     
        
    end % of static methods
    
    
end
