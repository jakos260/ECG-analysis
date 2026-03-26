% displcol.m
% script of triplot.m
% input and process new column
% NB: changes made may need to be replicated in setnode.m, which
% implements the effect of node as found by script nodesel.m
% 20140929


%if exist('column'), colprev=column; end
coltemp=get(sltri,'value');
column=round(coltemp);

set(sltri,'value',column);

pause(0.5)

[azim,elev]=view;

clf
if exist('contourcolor')
%     fun=VALS(:,column);
%     set(hs,'FaceVertexCData',fun)
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

if displmode==3,
    figure(3)
    set(marksigs,'Xdata',x(column),'Ydata',14*std_sigs(column)/max(std_sigs))
    figure(4)
end


% if gcf==10,
%     figure(9)
%     shiftmark=column-column_ref;
%     set(plt1,'xdata',imark+shiftmark,'ydata',SHOW(single_lead,column))
%     set(plt2,'xdata',imark+shiftmark,'ydata',sstd(column))
%     figure(10);
% end


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
    plot3(EDGENODES(:,1),EDGENODES(:,2),EDGENODES(:,3),'k','linewidth',1);
end

if exist('keep_edges'),
    if keep_edges==1,
        hold on
        if exist('edge1'), plot3(VER(edge1,1),VER(edge1,2),VER(edge1,3),'k-','linewidth',1),end
        if exist('edge2'), plot3(VER(edge2,1),VER(edge2,2),VER(edge2,3),'k-','linewidth',1), end
        if exist('edge3'), plot3(VER(edge3,1),VER(edge3,2),VER(edge3,3),'k-','linewidth',1), end
        if exist('edge4'), plot3(VER(edge4,1),VER(edge4,2),VER(edge4,3),'k-','linewidth',1) ,end
        if exist('elg_markers'),
            nelgmarkers=size(elg_markers,1);
            plot3(VER(elg_markers,1), VER(elg_markers,2),VER(elg_markers,3),'*k')
            for i=1:nelgmarkers,
                text(fac*VER(elg_markers(i),1), fac*VER(elg_markers(i),2),fac*VER(elg_markers(i),3),num2str(i),'fontsize',12)
            end
        end
        
    end
end

if exist('figsig','var'),
    % NB: replicate any changes in the script below in: setnode.m
    fignow=gcf;
    if exist('figsig')
        if ~isempty(figsig)
            if ~isempty(SIGS1),
                figure(figsig)
                clf
                plot(SIGS1(node,:))
                titel=['signal at node= ' num2str(node)];
                title(titel);
                xlabel('time  [ms]')
                ylabel('node voltage  [mV]')
                hold on
                plot([0 size(SIGS1,2)],[0 0],'g:')
                plot(column,SIGS1(node,column),'r*')
                if exist('SIGS2','var'),
                    if ~isempty(SIGS2),
                        plot(SIGS2(node,:),'r')
                        titel=[titel ' red: its estimate'];
                        title(titel);
                    end
                end
                
                %         if exist('dep'),
                %             if ~isempty(dep),
                %                 plot(dep(node),SIGS1(node,round(dep(node))),'k*')
                %                 plot(rep(node),SIGS1(node,round(rep(node))),'k*')
                %                 title(['SIG node ' num2str(node) ' dep: ' num2str(dep(node)) ...
                %                     ' rep: ' num2str(rep(node))]);
                %                 if exist('SRC'),
                %                     plot(SRC(node,1),SIGS1(node,round(SRC(node,1))),'m*')
                %                     plot(SRC(node,2),SIGS1(node,round(SRC(node,2))),'m*')
                %                 end
                %             end
                %         endend
            end
            
        end
        figure(fignow)
    end
end



