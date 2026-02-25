% make_square.m
% function [VER, ITRI]=make_square(nx,ny);
% generates triangular mesh for unit square with edges from origin to 
% [1 0 0] [0 1 0] with nodes nx and ny, respectively
% hence: its normal is oriented along z-axis

% 2011_11_29; A. van Oosterom

function [VER, ITRI]=make_square(nx,ny);
nver=nx*ny;
VER=zeros(nver,3);
ntri=2*(nx-1)*(ny-1);
ITRI=zeros(ntri,3);

iver=0;
for j=1:nx,
    for i=1:ny,
        iver=iver+1;
        VER(iver,:)=[(j-1)/(nx-1) (i-1)/(ny-1) 0];
    end
end

itri=0;
for i=1:ny-1,
   for j=1:nx-1;
       itri=itri+1;
       ITRI(itri,:)=[i+(j-1)*ny  i+j*ny  i+1+(j-1)*ny  ];
       itri=itri+1;
       ITRI(itri,:)=[i+j*ny i+1+j*ny  i+1+(j-1)*ny     ];
       
   end
end



