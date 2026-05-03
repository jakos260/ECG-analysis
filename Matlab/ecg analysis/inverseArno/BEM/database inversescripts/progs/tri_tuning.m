% tri_tuning.m
% tune triangulated surface by zooming in on area around a node
% define:
%        focus: the node around which the geometry is inspected
%        order: default 1; order is the neighbor edgedistance inspected
%                         
       
% if any change is made:
% run tricheck to test the result before storing it !!
% A. van Oosterom 20090616


if ~exist('focus'), focus=1; end
if ~exist('order'), order=1; end

ntri=size(ITRI,1);
nver=size(VER,1);

ED=edge_dist(ITRI,focus);
ned=size(ED,1);
ed=ED(ED(:,2)<=ones(ned,1)*order,1);

ITRI=ITRI(ismember(ITRI(:,1),ed)& ismember(ITRI(:,2),ed) & ismember(ITRI(:,3),ed),:);

VER=VER-ones(nver,1)*VER(focus,:);

newboxsize=max(norm3d(VER(ed,:)));



figure(1)
clf

node=focus;
VALS=VER;
triplot

hold on


figure(2)
clf
crossec


