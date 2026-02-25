% savesymmat.m
% save symmetrical matrix
% 20050926; A. van Oosterom
function savesymmat(file,M,extra)
[m,n]=size(M);
if m~=n, 'error: non square matrix', pause, end
if sum(sum(M~=M'))>0, 'error: non symmetric matrix', pause, end 
k=n*(n+1)/2;
all=[];
for j=1:n,
    all=[all;M(j:n,j)];
end

if nargin<3, 
    savemat(file,all),
else,
    savemat(file,all,extra)
end
    