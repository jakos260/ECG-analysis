% d/dphi of polar coordinates a(2)=theta, a(3)=phi
function der=ddphi(a)
der(1)=-sin(a(2))*sin(a(3));
der(2)=sin(a(2))*cos(a(3));
der(3)=0;

