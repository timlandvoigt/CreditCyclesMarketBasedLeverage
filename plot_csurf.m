function figlist=plot_csurf(Ygrid,Mgrid,QYmat,QMmat,prvec,val_other,linespec,printfile,figsin,mode,varargin)

if nargin>10
    legendstr=varargin{1};
end
if nargin>11
    levPts=varargin{2};
end


zlab={'Y','M'};

if mode==3
    RYsurfvals=squeeze(QYmat(:,:,4));
    RMsurfvals=squeeze(QMmat(:,:,4));
    
    RYsurfvals=100*prvec'*RYsurfvals;
    RYsurf=grid.LinearInterpFunction(Ygrid,RYsurfvals');
    
    RMsurfvals=100*prvec'*RMsurfvals;
    RMsurf=grid.LinearInterpFunction(Mgrid,RMsurfvals');    
        
    % Make Graphs
    surfct={RYsurf,RMsurf};
        
    figlist=cell(2,1);
    for f=1:length(surfct)
        surf=surfct{f}.plot3D(1,[2,3],val_other{f}(1),0);
        ax=surf.Parent;
        ax.XLabel.String=['LTV' ,zlab{f}];
        ax.YLabel.String=['LTW', zlab{f}];
        ax.ZLabel.String=['Spread ',zlab{f}];
        lev=surfct{f}.evaluateAt(val_other{f});
        levpt0=[val_other{f}(2:3),lev*2];
        levpt1=[val_other{f}(2:3),lev];
        hold(ax,'on');
        arrow3(levpt0,levpt1,'k-2');
        fig=ax.Parent;
        fig.PaperUnits='centimeters';
        fig.PaperOrientation='landscape';
        fig.PaperPosition=[2 2 25 16];
        if ~isempty(printfile)
            print(fig,'-dpdf','-r0',[printfile,zlab{f},'.pdf']);
        end
        figlist{f}=fig;
    end
    
elseif mode==2
%     RYsurfvals=squeeze(QYmat(:,:,3))-1;
%     RMsurfvals=squeeze(QMmat(:,:,3))-1;
    RYsurfvals=squeeze(QYmat(:,:,4));
    RMsurfvals=squeeze(QMmat(:,:,4));
        
    npl=size(prvec,2); % multiple columns --> plot for different states
    surfct=cell(2,npl);
    for p=1:npl
        prthis=squeeze(prvec(:,p));
        RYs=prthis'*RYsurfvals;
        surfct{1,p}=grid.LinearInterpFunction(Ygrid,RYs');
        
        RMs=prthis'*RMsurfvals;
        surfct{2,p}=grid.LinearInterpFunction(Mgrid,RMs');        
    end
    
    figlist=cell(2,1);
    figlist{2}=[];
    for f=1:1
        if ~isempty(figsin)
            thisfig=figsin{f};
            ax=get(thisfig,'CurrentAxes');
        else
            ax=[];
        end
        for p=1:npl
            surf=surfct{f,p}.plot2D(1,2,val_other{f},linespec(p),0,ax);
            surf{1}.LineWidth=2;
            if isempty(ax)
                ax=surf{1}.Parent;
            end
            levpoint=[val_other{f}(1),levPts(f,p),val_other{f}(2)];
            plot(ax,levPts(f,p),surfct{f,p}.evaluateAt(levpoint),[linespec{p},'o'],'LineWidth',2);
            ax.XLabel.String=['LTV ',zlab{f}];
            ax.YLabel.String=['Adjusted Mortgage Spread ',zlab{f}];
            ax.XLim=[0.55,0.75];
            ax.YLim=[0.0,0.06];
            ax.FontSize=14;
        end
        fig=ax.Parent;
        fig.PaperUnits='centimeters';
        fig.PaperOrientation='landscape';
        fig.PaperPosition=[2 2 25 16];
        if ~isempty(legendstr)
            legend(legendstr);
        end
        if ~isempty(printfile)
            print(fig,'-dpdf','-r0',[printfile,zlab{f},'.pdf']);
        end
        figlist{f}=fig;
    end
end    



end