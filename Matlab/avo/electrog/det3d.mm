% det3d.m
% function block=det3d(R1,R2,R3)
% block product i.e. determinant of three 3d vectors: R(:,1) R(:,2) R(:,3)
  function block=det3d(R1,R2,R3)

block=R1(:,1).*(R2(:,2).*R3(:,3)-R2(:,3).*R3(:,2)) - ...
      R1(:,2).*(R2(:,1).*R3(:,3)-R2(:,3).*R3(:,1)) + ...
      R1(:,3).*(R2(:,1).*R3(:,2)-R2(:,2).*R3(:,1));

