
global nodes

L   =   get(gca,'currentpoint');
ver=get(findobj(gca,'Type','patch'),'Vertices');
itri=get(findobj(gca,'Type','patch'),'faces');
if iscell(ver)
	ver=ver{end};
	itri=itri{end};
end
% intersections:
SELECT=linetris(ver,itri,L(1,:),L(2,:));
% only the facing ones are relavant
% SELECT=SELECT(SELECT(:,2)>0,:);

[small,is]=min(SELECT(:,5));
veris=itri(SELECT(is,1),:);
[small,ii]=min(sum((ver(veris,:)-(ones(1,3)'*(L(1,:)+SELECT(is,5)*(L(2,:)-L(1,:))))).^2'));
node=veris(ii);
delete(findobj('Tag','dot'));
line(ver(node,1),ver(node,2),ver(node,3),'Markersize',8,'Marker','o','Color','w','Linewidth',3,'Tag','dot');
triag = SELECT(is,1);

disp(['vertex ' num2str(node) '   of triangle   ' num2str(triag)])
% if isempty(SELECT(SELECT(:,2)>0,:))
% 	title(['node: ' num2str(node) ' of triangle ' num2str(triag)],'color','r');
% else
% 	title(['node: ' num2str(node) ' of triangle ' num2str(triag)],'color','b');
% end

if ~exist('nodes')
    nodes=[];
end
nodes = [nodes node];