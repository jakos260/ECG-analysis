% loadsymmat.m
% load symmetric matrix
% 20050926; A. van Oosterom
function [M,extra]=loadsymmat(file,extra)

if nargin>1,
   all=loadmat(file,extra);
else,
    all=loadmat(file);extra=[];
end

[k,l]=size(all);
if l~=1, 'error: file should contain a single column only', stop, end
n=(-1+sqrt(1+8*k))/2;


if n~=round(n), 'error: non symmetric matrix', stop, end 
M=zeros(n);
ibeg=0;; iend=0;
k=n*(n+1)/2;
for j=1:n,
    ibeg=iend+1; iend=iend+n-j+1;
    M(j:n,j)=all(ibeg:iend);
end
M=M+tril(M,-1)';
