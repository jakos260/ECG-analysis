% rotation_matrix.m
% ROT=rotation_matrix(alpha,u);
% compute matrix for rotation about axis u; passing through the origin
% alpha desired rotation; u: row vector (=rotation axis)
% specify alpha in radians
% 
% A. van Oosterom 2014_08_19


function ROT=rotation_matrix(alpha,u)
if size(u,2)==1, u=u'; end
u=u/norm(u);
U=u'*u;
N=[0 -u(3) u(2); u(3) 0 -u(1); -u(2) u(1) 0;];
ROT=cos(alpha)*eye(3)+sin(alpha)*N+(1-cos(alpha))*U;





