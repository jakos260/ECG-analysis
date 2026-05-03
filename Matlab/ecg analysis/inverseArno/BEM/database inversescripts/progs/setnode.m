% file setnode.m
% display the node in triplot as, e.g., identified by nodesel

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
