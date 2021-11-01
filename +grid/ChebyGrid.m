classdef ChebyGrid < grid.StateSpaceGrid
    
    properties (Abstract, SetAccess=protected)
        Powers
        Terms
        ChebyPointmat
    end   
    
%--------------------------------------------------------------------------
% Static methods (miscellaenous functions for chebychev polynomials)
%--------------------------------------------------------------------------
    methods (Static)
        
        % function for zeros of cheb polyn of degree D
        function nodes=getChebyZeros(D)   
            nodes = zeros(D,1);
            for i=1:D
                nodes(i) = -cos((2*i-1)*pi/(2*D));
            end            
        end
        
        % recursive evaluation routine 
        function eval=chebyraw(x,n)          
            np=length(x);
            
            t1 = ones(np,1);
            if n==0
                eval=t1;
                return
            end
            t2 = x;
            eval = t2;
            for i=2:n
                t3 = 2*x.*t2 - t1;
                eval = t3;
                t1=t2;
                t2=t3;
            end           
        end
        
        % precompute terms without coefs and summing up
        function f=evalcheby_precomp(x,powers)           
            [np,dim]=size(x);
            nter=size(powers,1);
            f=zeros(np,nter);
            for t=1:nter
                temp=ones(np,1);
                for d=1:dim
                    temp=temp.*grid.ChebyGrid.chebyraw(x(:,d),powers(t,d));
                end
                f(:,t) = temp;
            end        
        end
        
        % linear transformation [-1,1] to actual state space
        function pout=chebyToSS(pin,stBounds)
            npt=size(pin,1);
            pout=0.5*(pin.*(ones(npt,1)*(stBounds(2,:)-stBounds(1,:))) ...
                    + ones(npt,1)*(stBounds(2,:) + stBounds(1,:)));
        end
        
        % linear transformation actual state space to [-1,1]       
        function pout=SSToCheby(pin,stBounds)
            npt=size(pin,1);
            pout=(2*pin-ones(npt,1)*(stBounds(2,:)+stBounds(1,:)))...
                ./(ones(npt,1)*(stBounds(2,:)-stBounds(1,:)));
        end
        
        % complete combinations for Smolyak's rule
        function A=completeCombinations(d,c)         
            % initialize
            basevec=(1:c-1)';
            A0=basevec;
            
            for i=2:d
                ml=size(A0,1);
                attach=kron(basevec,ones(ml,1));
                A1=[repmat(A0,c-1,1),attach];
                compsum=(sum(A1,2)<=c);
                A0=A1(compsum,:);
            end
            
            A=A0; 
        end
        
        % Smolyak's rule
        function [points,powmat]=smolyakMe(d,q)
            
            % degree vector
            degvec=2.^((2:q-d+1)-1) + 1;
            degvec=[1,degvec];
            
            % nested chebychev extrema up to max degree
            % incrementally sorted for indexing
            chex=-cos(pi/2);
            for k=2:length(degvec)
                deg=degvec(k);
                tempex=zeros(1,deg);
                for i=1:deg
                    tempex(i)=-cos( pi*(i-1)/(deg-1) );
                end
                chex=[chex,setdiff(tempex,chex)];
            end
            
            % complete combinations up to degree q for dim d
            cc=grid.ChebyGrid.completeCombinations(d,q);
            % discard rows with total degree < q-d+1
            cc=cc(sum(cc,2)>=q-d+1,:);
            scc=size(cc,1);
            
            % now union of actual points
            % first indices
            indexmat=ones(1,d);
            for i=1:scc
                tmp=grid.StateSpaceGrid.makeCombinations(degvec(cc(i,:)));
                % now add new points/terms
                tf=~ismember(tmp,indexmat,'rows');
                indexmat=[indexmat;tmp(tf,:)];
            end
            % actual points
            points=chex(indexmat);
            points=sortrows(points,1:d);
            
            % power matrix
            powmat=indexmat-1;            
        end
        
    end
%--------------------------------------------------------------------------
% End static methods 
%--------------------------------------------------------------------------
    
end