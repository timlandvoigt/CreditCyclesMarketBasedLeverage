classdef HelperCollection
   
    properties
       dummy         
    end
    
    methods
        function obj=HelperCollection()
            obj.dummy=1;
            disp('No need to create an instance of this class.');
        end
    end
    
    methods (Static)
        
        function bigstr=combineStructs(strArray)
            % strArray must be cell array of structs
            
            bigstr=struct;
            Ninstr=length(strArray);
            for i=1:Ninstr
                strtmp=strArray{i};
                fntmp=fieldnames(strtmp);
                for j=1:length(fntmp)
                    thisfield=fntmp{j};
                    bigstr.(thisfield)=strtmp.(thisfield);
                end
            end
        end
        
        function [indlist,effnames]=makeListFromNames(hashmap,names)            
            n=length(names);
            indlist=zeros(n,1);
            for i=1:n
                tmp=hashmap.get(names{i});
                if ~isempty(tmp)
                    indlist(i)=tmp;
                else
                    indlist(i)=0;
                end
            end
            effnames=names(indlist~=0);
            indlist=indlist(indlist~=0);            
        end
        
        function dum=tableExport(filename, varnames, data)
            
            dum=[];
            
            ncol=length(varnames);
            fid=fopen(filename,'w');
            for c=1:ncol-1
                fprintf(fid,'%s,',varnames{c});
            end
            fprintf(fid,'%s\n',varnames{ncol});
            fclose(fid);
            
            dlmwrite(filename,data,'-append');
            
        end
              
        function dum=scatterPoints2D(X,Y)
            dum=[];
            
            bY=max(Y);
            scale=100/bY;
            figure; hold on;
            scatter(X(:,1),X(:,2),Y*scale);            
        end
        
        function [x,w]=lgwt(N,a,b)
            
            % lgwt.m
            %
            % This script is for computing definite integrals using Legendre-Gauss
            % Quadrature. Computes the Legendre-Gauss nodes and weights  on an interval
            % [a,b] with truncation order N
            %
            % Suppose you have a continuous function f(x) which is defined on [a,b]
            % which you can evaluate at any x in [a,b]. Simply evaluate it at all of
            % the values contained in the x vector to obtain a vector f. Then compute
            % the definite integral using sum(f.*w);
            %
            % Written by Greg von Winckel - 02/25/2004
            N=N-1;
            N1=N+1; N2=N+2;
            
            xu=linspace(-1,1,N1)';
            
            % Initial guess
            y=cos((2*(0:N)'+1)*pi/(2*N+2))+(0.27/N1)*sin(pi*xu*N/N2);
            
            % Legendre-Gauss Vandermonde Matrix
            L=zeros(N1,N2);
            
            % Derivative of LGVM
            Lp=zeros(N1,N2);
            
            % Compute the zeros of the N+1 Legendre Polynomial
            % using the recursion relation and the Newton-Raphson method
            
            y0=2;
            
            % Iterate until new points are uniformly within epsilon of old points
            while max(abs(y-y0))>eps
                
                
                L(:,1)=1;
                Lp(:,1)=0;
                
                L(:,2)=y;
                Lp(:,2)=1;
                
                for k=2:N1
                    L(:,k+1)=( (2*k-1)*y.*L(:,k)-(k-1)*L(:,k-1) )/k;
                end
                
                Lp=(N2)*( L(:,N1)-y.*L(:,N2) )./(1-y.^2);
                
                y0=y;
                y=y0-L(:,N2)./Lp;
                
            end
            
            % Linear map from[-1,1] to [a,b]
            x=(a*(1-y)+b*(1+y))/2;
            
            % Compute the weights
            w=(b-a)./((1-y.^2).*Lp.^2)*(N2/N1)^2;            
        end
    end
end