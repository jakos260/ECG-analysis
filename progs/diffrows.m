function D=diffrows(M, mode)
% function D=diffrows(M,mode)
% differentiate rows of matrix M
% accurate for local parabolae
% mode=2  (default) assume parabolic continuation beyond both ends
% mode=1            assume constant value beyond both ends, i.e. :
%                   y(0)=y(1) and y(n+1)=y(n)
% mode=0            assign zero gradient to both ends
% 20060905

if nargin==1 | mode>2 | mode<0,
    mode=2;
end


n=size(M,2);

if mode==2,
   % end effect: parabolic continuation beyond both ends;
   %M=[3*M(:,1)-3*M(:,2)+M(:,3)  M  3*M(:,n)-3*M(:,n-1)+M(:,n-2);];
   M=[3*M(:,1)-3*M(:,2)+M(:,3)  M  3*M(:,n)-3*M(:,n-1)+M(:,n-2);];
end

if mode==1,
   % end effect: function constant; beyond both ends
   M=[M(:,1) M  M(:,n)];
end

if mode==0,
    % zero gradient assigned to both ends
    M=[M(:,2) M  M(:,n-1)];  
end


D=(M(:,3:n+2)-M(:,1:n))/2;


