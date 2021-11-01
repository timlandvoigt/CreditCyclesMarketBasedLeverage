classdef MixedFunction < grid.ApproxFunction

   
    properties (SetAccess=protected)
        % inherited properties (abstract in superclass)
        SSGrid
        Nof
        Vals      
        % mixed specific properties
        Nofcheb
        Nofspl
        Nterm
        Chebfct
        Splfct
        Spldegvec
        Splextrap
        Inifit
    end    
    
    methods
        function mf=MixedFunction(mfgrid,vals,spldegvec,extrap) 
            mf.SSGrid=mfgrid;
            % determine number of points and functions
            [npt,nof]=size(vals);
            if npt~=mf.SSGrid.Npt
                error('Value matrix must have dimensions (Spline.Npt*Cheby.Npt x Nof)');
            end
            mf.Nof=nof;      
            mf.Vals=vals;
            mf.Spldegvec=spldegvec;
            mf.Splextrap=extrap;
            % number of functions for cheby function    
            mf.Nofcheb=mf.Nof*mf.SSGrid.Tensorgrid.Npt;
            % initial call for fitTo
            mf.Inifit=true;
            mf=fitTo(mf,vals);
            mf.Inifit=false;
        end
        
        function mf=set.SSGrid(mf,mfgrid)
            % check that grid object is of class MixedGrid
            if ~isa(mfgrid,'grid.MixedGrid')
                error('StateSpaceGrid must be a MixedGrid');
            end
            mf.SSGrid=mfgrid;            
        end
        
        function mf=fitTo(mf,points)
            npc=mf.SSGrid.Chebygrid.Npt;
            nps=mf.SSGrid.Tensorgrid.Npt;
            % reshape values for cheby function
            valscheb=reshape(points,npc,mf.Nofcheb);
            if mf.Inifit
                % create cheby function
                mf.Chebfct=grid.ChebyFunction(mf.SSGrid.Chebygrid,valscheb);
                % number of functions for spline function
                mf.Nterm=size(mf.Chebfct.Coefs,1);
                mf.Nofspl=mf.Nof*mf.Nterm;
            else
                % or rerun fit
                mf.Chebfct=mf.Chebfct.fitTo(valscheb);
            end
            % reshape coefficient matrix as values for spline function
            coefsspl=reshape(mf.Chebfct.Coefs,mf.Nterm,nps,mf.Nof);
            % the coef matrix has dim (nps x nof*nterm), where the columns
            % are in blocks of terms per function; this is for fast
            % evaluation
            coefsspl=reshape(permute(coefsspl,[1,3,2]),mf.Nterm*mf.Nof,nps)';
            if mf.Inifit
                % create spline function
                mf.Splfct=grid.SplineFunction(mf.SSGrid.Tensorgrid,coefsspl,mf.Spldegvec,mf.Splextrap);          
            else
                % or rerun fit
                mf.Splfct=mf.Splfct.fitTo(coefsspl);
            end
        end
        
        function vals=evaluateAt(mf,points)
            nds=mf.SSGrid.Tensorgrid.Ndim;
            [np,ndim]=size(points);
            if ndim~=mf.SSGrid.Ndim
                error('Point matrix must have dimensions (#points x Spline.Npt*Cheby.Npt)');
            end
            % first evaluate spline function, polcoefs has dim (nof*nterm x np),
            % where rows are in blocks of terms per function
            polcoefs=mf.Splfct.evaluateAt(points(:,1:nds));
            % get polynomial terms at points, terms has dimension (nterm x np)
            terms=mf.Chebfct.computeTerms(points(:,nds+1:end));
            % replicate terms to get dim of polcoefs
            terms=repmat(terms,mf.Nof,1);
            % coefficients times terms by point and function
            valbig=polcoefs.*terms;
            % reshape to sum up terms, result has dim (nof x np)
            valbig=reshape(valbig,mf.Nterm,mf.Nof,np);
            vals=squeeze(sum(valbig,1));
            if isrow(vals)
                vals=vals';
            end
        end
        
    end
    
    
end