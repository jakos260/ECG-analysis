function B = forudlinf(VER,ITRI,OBS)
% forudlinf.m
%
% solve the forward problem for a closed UDL source
% in an infinite homogeneous medium
%
% B = forudlinf(VER,ITRI,OBS)
% VER, ITRI: triangulated surface
% OBS: observation points (Nobs x 3-matrix)
% B: Nobs x Nver-matrix of distributed solid angles
% 

if nargin<3
	error('Not enough input arguments.');
end

if min(ITRI(:))<1 | max(ITRI(:))>size(VER,1)
	error('Index in ITRI out of bounds');
end

s = 1/(2*pi);
B = s * forudlinf_i(VER',ITRI',OBS');
