classdef CompleteChebyGrid < grid.TensorGrid & grid.ChebyGrid
   
    properties (SetAccess=protected)
        % inherited from ChebyGrid (abstract)
        Powers
        Terms
        ChebyPointmat
    end    
    
  methods
        % constructor
        function ssg=CompleteChebyGrid(stateBounds,dimvec,compdeg)
            % quick check on inputs
            nr=size(stateBounds,1);
            if nr~=2
                error('StateBounds must be a 2xNdim matrix');
            end
            % build grids for each dimension
            ndim=length(dimvec);
            unigrids=cell(ndim,1);
            cheb_unigrids=cell(ndim,1);
            for i=1:ndim
                cheb_unigrids{i}=grid.ChebyGrid.getChebyZeros(dimvec(i));
                unigrids{i}=grid.ChebyGrid.chebyToSS(cheb_unigrids{i},stateBounds(:,i));
            end
            % call TensorGrid constructor with these
            ssg=ssg@grid.TensorGrid(unigrids);
            % set remaining properties
            ssg.ChebyPointmat=grid.ChebyGrid.SSToCheby(ssg.Pointmat,ssg.StateBounds);
            % power matrix and terms
            powmat=grid.StateSpaceGrid.makeCombinations(dimvec)-1;
            compsum=sum(powmat,2);
            ssg.Powers=powmat(compsum<=compdeg,:);
            ssg.Terms=grid.ChebyGrid.evalcheby_precomp(ssg.ChebyPointmat,ssg.Powers);
            ssg.Type='CompleteChebyGrid';
        end    
  end
    
end