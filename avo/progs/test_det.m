% test_det.m
R=rand(3);

%R=eye(3)

cr=cross(R(2,:),R(3,:));
cross(R(1,:),cr)

R(2,:)*(R(1,:)*R(3,:)')-R(3,:)*(R(1,:)*R(2,:)')



cr=cross(R(1,:),R(2,:));
cross(cr,R(3,:))
-R(1,:)*(R(2,:)*R(3,:)')+R(2,:)*(R(1,:)*R(3,:)')

 


% 
% E=[R(2,:)-R(1,:);
%    R(3,:)-R(2,:); 
%    R(1,:)-R(3,:);];
% 
% n=cross(E(1,:),E(2,:));
% 
% % n=n/norm(n);
% 
block = det3d(R(1,:),R(2,:),n)


cr=cross(R(1,:),R(2,:));
n*cr'


% 
% 
% 
% 
% cr=cross(R(1,:),R(2,:));
% 
% cr2=cross(R(2,:),cr)
% 
% R(1,:)*(-cr1+cr2)'
