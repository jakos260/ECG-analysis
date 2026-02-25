% script: timediffmat.m
% create matrix for performing differentiatiation with respect
% to time (row elements)
% method: mean of forward and backward step differences
% useage, applied to PHI: S=PHI*D
 
function D=timediffmat(n)
D=zeros(n,n);
D(1,1)=-.5;
D(2,1)=.5;
for j=2:n-1,
for i=1:n,
D(i,j)=0.;
end
D(j-1,j)=-.5;
D(j+1,j)=.5;
end
D(n-1,n)=-.5;
D(n,n)=.5;

