function varargout=calcAnisotropicDist(ADJ,ITRI,anisotropyRatio)

% Calculate the distance matrix based on a (given) adjacency matrix in 
% which the connections are adapted. the connections though the wall are
% multiplied by the anisotropyRatio, ususally to make them longer due to
% the fact that the velocity through the wall is slower. The
% anisotropyRatio representes therefore the ratio between the axial and
% transverse myocardial cell velocity.

% Peter van Dam 11-03-2008


disp('Busy!!!!!')
if size(ITRI,1)==size(ITRI,2) || size(ITRI,2)~=3
	error('second parameter should be the indexes of the triangles')
elseif size(ADJ,1)~=size(ADJ,2) 
	error('first parameter should be the adjacency matrix')
end
buur=graphdist(ITRI);
ADJW=ADJ;
ADJW(buur==0)=ADJW(buur==0)*anisotropyRatio;
varargout{1}=graphdist(ADJW);
if nargout==2
	varargout{2}=ADJW;
end
