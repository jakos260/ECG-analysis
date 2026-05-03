% rhoek.m
 %[sa,intri]=rhoek(r1,r2,r3)
% compute solid angle sa subtended by a triangle,
% defined by (row)vectors r1,r2,r3 specifying its vertices
% as viewed from the origin

% if, when viewed from the origin, r1,r2,r3 are seen in clockwise order, the sign is
% negative
% use solida if multiple solid angles need to be computed 

% if the origin (observation point) is NOT in the plane of the triangle: intri=0
% else  if the origin is external to the triangle                      : intri=1
%       if the origin lies on an edge of the  triangle                 : intri=2
%       if the origin is an internal point of the triangle             : intri=3

% basis: A. van Oosterom and J. Strackee; Trans IEEE BME-30/2 Feb 1983; 125-126

       function [sa,intri]=rhoek(r1,r2,r3)
       R=[r1;r2;r3];
       block=det(R);
       n1=norm(r1);
       n2=norm(r2);
       n3=norm(r3);
       dot12=r1*r2';
       dot23=r2*r3';
       dot13=r1*r3';
       denom=n1*n2*n3+dot12*n3+dot23*n1+dot13*n2;
       if abs(block) < 1.e-10,
             sa=0;
	     intri=2;
             if denom < -1.e-10, intri=3;
             end
	     if denom > 1.e-10, intri=1;
             end
       else
       sa=-2*atan2(block,denom);
       intri=0;
       end
       