% file setnode.m
% display the node in triplot as, e.g., identified by nodesel

% changes made may need to be replicated in displcol
% A. van Oosterom, 2016_05_10

if ~exist('iss'),
    trisnode=find(ITRI(:,1)==node|ITRI(:,2)==node|ITRI(:,3)==node);
    if ~isempty(trisnode)
        iss=trisnode(1);
        is=1;
        SELECT(is,2)=1;
    end
end

if exist('is')
    if exist('ui4'), set(ui4,'string',sprintf('%3.3f',fun(node))); end
end
if exist('ui5') set(ui5,'string',num2str(node)); end

if exist('ui14')& exist('iss'), set(ui14,'string',num2str(iss)); end

if exist('is'),
    if SELECT(is,2)>= 0 ,
        if exist('ui5'), set(ui5,'foreground','b');end
        if exist('ui14'),set(ui14,'foreground','b');end
    else
        set(ui5,'foreground','r');
        if exist('ui14'),set(ui14,'foreground','r');end
    end
end

if exist('ui13') & exist('ui14'), set([ui13;ui14],'vis','on'), end
if exist('hcbarpos'),set(colorbar,'position',hcbarpos); end

if exist('ht'), set(ht,'pos',VER(node,:));end

if exist('figsig','var'),
    % NB: replicate any changes in the script below in:  displcol.m
    if ~empty(figsig),
        fignow=gcf;
        
        if ~isempty(SIGS1),
            figure(figsig)
            clf
            plot(SIGS1(node,:))
            titel=['SIG1 at node= ' num2str(node)];
            title(titel);
            xlabel('time  [ms]')
            
            ylabel('node voltage  [mV]')
            hold on
            plot([0 size(SIGS1,2)],[0 0],'g:')
            plot(column,SIGS1(node,column),'r*')
            
            if exist('SIGS2','var'),
                if ~isempty(SIGS2),
                    plot(SIGS2(node,:),'r')
                    titel=[titel ' red: SIG2'];
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
            %         end
            
        end
        
    end
    figure(fignow)
end





