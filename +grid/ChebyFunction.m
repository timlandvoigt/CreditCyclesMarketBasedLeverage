classdef ChebyFunction < grid.ApproxFunction
    
    properties (SetAccess=protected)
        % inherited properties (abstract in superclass)
        SSGrid
        Nof
        Vals        
        % chebychev specific properties
        Coefs
        Maxpow      
    end
    
    methods
        % constructor
        function cf=ChebyFunction(ssgrid,vals)
           cf.SSGrid=ssgrid;
           cf.Maxpow=max(cf.SSGrid.Powers)';
           cf=fitTo(cf,vals);
        end
        
        function cf=set.SSGrid(cf,ssg)
            % check that grid object is of class ChebyGrid
            if ~isa(ssg,'grid.ChebyGrid')
                error('StateSpaceGrid must be a ChebyGrid');
            end
            cf.SSGrid=ssg;
        end
        
        function cf=fitTo(cf,vals)
            [npt,nof]=size(vals);
            if npt~=cf.SSGrid.Npt
                error('Value matrix must have dimensions (Npt x Nof)');
            end
            cf.Nof=nof;
            cf.Vals=vals;
            terms=cf.SSGrid.Terms;
            % get coefficients for all functions
            % this interpolates if collocation, otherwise least-squares fit
            cf.Coefs=terms\vals;
        end
        
        function vals=evaluateAt(cf,points)
            ndim=size(points,2);
            if ndim~=cf.SSGrid.Ndim
                error('Point matrix must have dimensions (#points x Ndim)');
            end
            % transform points to [-1,1] interval
            cheby_points=grid.ChebyGrid.SSToCheby(points,cf.SSGrid.StateBounds);
            % force into bounds
            cheby_points(cheby_points>1)=1;
            cheby_points(cheby_points<-1)=-1;
            % call to mex evaluation routine
            vals=evalchebyC(cf.Coefs,cheby_points',cf.SSGrid.Powers',cf.Maxpow,[]);
        end
        
        function terms=computeTerms(cf,points)
            ndim=size(points,2);
            if ndim~=cf.SSGrid.Ndim
                error('Point matrix must have dimensions (#points x Ndim)');
            end
            % transform points to [-1,1] interval
            cheby_points=grid.ChebyGrid.SSToCheby(points,cf.SSGrid.StateBounds);
            % force into bounds
            cheby_points(cheby_points>1)=1;
            cheby_points(cheby_points<-1)=-1;
            % only pass coefs for first function (coefs are irrelevant for
            % this method)
            [~,terms]=evalchebyC(cf.Coefs(:,1),cheby_points',cf.SSGrid.Powers',cf.Maxpow,[]);
        end
       
    end
    
    
end