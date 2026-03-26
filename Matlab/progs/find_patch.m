% find_patch.m
% [patchnodes,patchtris]=find_patch(ITRI,border,seed)
% Finds all node indices of a patch of a triangulated surface outlined by all nodes
% enclosing the patch (including those of the border)
% as well as the indices of the triangles constituting the patch.
% The border may consist of multiple, non-intersecting CLOSED, contours
% seed: one internal node of the desired patch (which contains at least one
% internal node)
% for frequent use of large ADJ and DIST matrices: define them as global

% A. van Oosterom; 2013_10_18

function [patchnodes,patchtris]=find_patch(ITRI,border,seed) % patchtris = patch triangle indices (column vector);

% force: column vectors
if size(border,1)<size(border,2), border=border'; end
if size(seed,1)<size(seed,2), seed=seed'; end

border=unique(border);
patchnodes=[border;seed];
ADJ=graphdist(ITRI);
nn=size(ADJ,1);
nb=size(border,1);
ADJ(border,:)=zeros(nb,nn);
ADJ(:,border)=zeros(nn,nb);
DIST=graphdist(ADJ);
DIST(:,seed);
intnodes=find(DIST(:,seed)>zeros(nn,1));
patchnodes=unique([intnodes; patchnodes;]);
    
% find triangle indices of the path
ntri=size(ITRI,1);
patchtris=find(sum(ismember(ITRI(:,[1 2 3]),patchnodes),2)==3*ones(ntri,1));



