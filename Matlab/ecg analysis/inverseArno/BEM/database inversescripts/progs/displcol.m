% displcol.m
% script of triplot.m
% input and process new column
% 20130826

coltemp=get(sltri,'value')
column=round(coltemp)

set(sltri,'value',column)

pause(0.5)

[azim,elev]=view;

clf
if exist('contourcolor'),
    triplot_contour,
else
    triplot
end



if ~exist('displmode'), displmode=0; end

if displmode==1;
    if exist('ECGG'),
        % for newinit:
        figure(2);
        clf
        sigplot(ECG,' ',LAY,.7,'b',0)
        dep=TIMS(:,column);
        gets;
        PHI=AA*S(:,1:tmaxi);
        rd=norm(ECG-PHI,'fro')/norm(PHI,'fro');
        sigplot(PHI,' ',LAY,.7,'r',0),
        % display rd
        ui35=uicontrol('style','text');
        uibox35=[.7 .9 .12 .05];
        set(ui35,'units','norm','position',uibox35,'string',...
            [' rd=' num2str(rd)]);
    end
end

if displmode==2,
    % for ecganal
    set(marktim,'xdata',[0 0],'ydata',[ty(column) ty(column)], ...
        'zdata',[scal*PHI(leadnos(node),column+tbeg-1)+z0-.2 ...
        scal*PHI(leadnos(node),column+tbeg-1)+z0+.2],'color','r');
end

set(ui17,'string',num2str(column))

if exist('ui4'),
    set(ui4,'string',sprintf('%3.3f',fun(node)));
end

if exist('ELPOS'),
    if ~isempty(ELPOS)
        show_elecs
    end
end

if exist('EDGENODES'),
    hold on
    plot3(EDGENODES(:,1),EDGENODES(:,2),EDGENODES(:,3),'k','linewidth',1.5);
end

if exist('keep_edges'),
    if keep_edges==1,
        hold on
        plot3(VER(edge1,1),VER(edge1,2),VER(edge1,3),'k-','linewidth',1.5)
        plot3(VER(edge2,1),VER(edge2,2),VER(edge2,3),'k-','linewidth',1.5)
        plot3(VER(edge3,1),VER(edge3,2),VER(edge3,3),'k-','linewidth',1.5)
    end  
end


