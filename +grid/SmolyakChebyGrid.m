classdef SmolyakChebyGrid < grid.ChebyGrid & grid.StateSpaceGrid
   
    properties (SetAccess=protected)
        % inherited from StateSpaceGrid (abstract)
        Ndim
        Npt
        Pointmat 
        Type
        % inherited from ChebyGrid (abstract)
        Powers
        Terms
        ChebyPointmat
    end    
    
  methods
        % constructor
        function ssg=SmolyakChebyGrid(stateBounds,smDeg)
            ssg.StateBounds=stateBounds;
            ssg.Ndim=size(ssg.StateBounds,2);
            % generate point and powers matrix
            [points,powers]=grid.ChebyGrid.smolyakMe(ssg.Ndim,smDeg);
            ssg.Npt=size(points,1);
            ssg.ChebyPointmat=points;
            ssg.Powers=powers;
            ssg.Pointmat=grid.ChebyGrid.chebyToSS(points,ssg.StateBounds);
            ssg.Terms=grid.ChebyGrid.evalcheby_precomp(ssg.ChebyPointmat,ssg.Powers);
            ssg.Type='SmolyakChebyGrid';
        end
  end
    
    
end