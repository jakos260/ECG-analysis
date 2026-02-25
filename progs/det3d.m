% det3d.m
% function block=det3d(R1,R2,R3)
% each row: block products: twice the volumes of the paralellepipids spanned from the origin to the row vectors R1,R2,R3
% i.e. the determinants, of the set of three 3d-vectors: R1(:,1:3) R2(:,1:3) R3(:,1:3)  
% For size(R)=3,3, det(R) gives the same answer  for a single tetraheder, including sign;
% A. van Oosterom

  function block=det3d(R1,R2,R3)

block=R1(:,1).*(R2(:,2).*R3(:,3)-R2(:,3).*R3(:,2)) - ...
      R1(:,2).*(R2(:,1).*R3(:,3)-R2(:,3).*R3(:,1)) + ...
      R1(:,3).*(R2(:,1).*R3(:,2)-R2(:,2).*R3(:,1));

