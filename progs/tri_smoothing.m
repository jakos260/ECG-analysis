% tri_smoothing.m
% ver= tri_smoothing(VER,ITRI,node,mode,alpha);
% smoothing local mesh around node by adjusting its position
% mode==1: move node by a factor alpha toward center of its direct neighbours; repeat if
%          desired
% mode=2: project node on the local spherical cap estimated by that from
%         the node's direct neighbours
% mode=3: project node on the local spherical cap estimated by that from
%         node's current position AND that of its direct neighbours

% date: 20101105

function ver=tri_smoothing(VER,ITRI,node,mode)

if nargin<4, mode=1; end

bver=buren(ITRI,node);

if mode==1,
   if nargin<5, alpha=0.5;end
   ver=alpha*VER(node,:)+(1-alpha)*mean(VER(bver,:));
end

if mode==2,
   USE=VER(bver,:);
   [r,cc]=fit_sphere(USE);
   ver=VER(node,:);
   ver=cc + r/norm3d(vernode-cc)*(vernode-cc);
end

if mode==3,
   USE=VER([node; bver],:);
   [r,cc]=fit_sphere(USE);
   vernode=VER(node,:);
   ver=cc + r/norm3d(vernode-cc)*(vernode-cc);
   
end






