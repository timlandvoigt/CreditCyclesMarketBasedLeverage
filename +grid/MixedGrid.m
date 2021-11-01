classdef MixedGrid < grid.StateSpaceGrid
% mixed approximation with spline and chebychev polynomial basis
% to create an object of this class, two grids need to be specified
% 1: a TensorGrid as basis for the spline part
% 2: a ChebyGrid as basis for the polynomial part
% The overall grid is the tensor product of points of both grids. 


   properties (SetAccess=protected)
        % inherited properties (abstract)
        Ndim
        Npt
        Pointmat 
        Type
        % mixed grid specific properties
        Chebygrid
        Tensorgrid
   end    
    
   methods
       % constructor
       function mssg=MixedGrid(tsgrid,chgrid) 
          % set partial grids (check in setter methods)
          mssg.Tensorgrid=tsgrid;
          mssg.Chebygrid=chgrid;
          % combined dimensions and points
          mssg.Ndim=tsgrid.Ndim+chgrid.Ndim;
          mssg.Npt=tsgrid.Npt*chgrid.Npt;
          mssg.StateBounds=[tsgrid.StateBounds,chgrid.StateBounds];
          % make combined Pointmat using forward ordering
          indmat=grid.StateSpaceGrid.makeCombinations([tsgrid.Npt,chgrid.Npt]);
          mssg.Pointmat=[tsgrid.Pointmat(indmat(:,1),:),chgrid.Pointmat(indmat(:,2),:)];
          mssg.Type='MixedGrid';
       end

      function mssg=set.Tensorgrid(mssg,tsgrid)
            % check that grid object is of class TensorGrid
            if ~isa(tsgrid,'grid.TensorGrid')
                error('First argument must be a TensorGrid');
            end
            mssg.Tensorgrid=tsgrid;           
      end
       
       function mssg=set.Chebygrid(mssg,chgrid)
            % check that grid object is of class ChebyGrid
            if ~isa(chgrid,'grid.ChebyGrid')
                error('Second argument must be a ChebyGrid');
            end
            mssg.Chebygrid=chgrid;           
       end
       
   end
   
   
end