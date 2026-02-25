

function [newVER ROT]=doRot(VER,azel,alpha)


u = azel(:)/norm(azel);

alph = alpha*pi/180;
cosa = cos(alph);
sina = sin(alph);
vera = 1 - cosa;
x = u(1);
y = u(2);
z = u(3);
ROT= [cosa+x^2*vera x*y*vera-z*sina x*z*vera+y*sina; ...
       x*y*vera+z*sina cosa+y^2*vera y*z*vera-x*sina; ...
       x*z*vera-y*sina y*z*vera+x*sina cosa+z^2*vera]';


 x = VER(:,1);
 y = VER(:,2);
 z = VER(:,3);
   
origin=[0 0 0];   
   
[m,n] = size(x);
newxyz = [x(:)-origin(1), y(:)-origin(2), z(:)-origin(3)];
newxyz = newxyz*ROT;
newx = origin(1) + reshape(newxyz(:,1),m,n);
newy = origin(2) + reshape(newxyz(:,2),m,n);
newz = origin(3) + reshape(newxyz(:,3),m,n);


newVER=[newx newy newz];