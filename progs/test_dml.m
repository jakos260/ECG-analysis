% test_dml
% A. van Oosterom;



clear

%check dml outcome at center of a rectangle

% specify vertices of the rectangle; sides a and b

% all results expressed in units(4*pi)   %%

a=1; b=2;

VER=[-a/2 -b/2 0;
     -a/2  b/2 0;
      a/2  b/2 0;
      a/2 -b/2 0;];

ITRI=[1 2 3; 
      2 4 3];

obs=[0 0 0];

phi_center=dml(VER,ITRI,obs);

phi_center=sum(sum(phi_center))
c=sqrt(a^2+b^2);


compare=a*log((c+b)/(c-b))+ b*log((c+a)/(c-a))


if b==a,
  
   compare=2*a*log( (sqrt(2)+1)/(sqrt(2)-1))
   % = 3.5255*a   ??
   2*a*log(3+2*sqrt(2))
    
end





