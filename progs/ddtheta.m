% d/dpol with respect to polar angle theta at a(2)=theta, a(3)=phi
function der=ddtheta(a)
der(1)=cos(a(2))*cos(a(3));
der(2)=cos(a(2))*sin(a(3));
der(3)=-sin(a(2));

