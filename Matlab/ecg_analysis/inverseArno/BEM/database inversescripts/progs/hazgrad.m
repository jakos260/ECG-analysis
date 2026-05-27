function G=hazgrad(par,tvals)
tvals=tvals';
f1= (1./(1.+par(1)./100*(par(2)-tvals)));
f2= f1.^par(3);
f3=-f1.^(-2);
G(1,:)= par(3).*f1.^(par(3)-1).*f3.*(par(2)-tvals)/100;
G(2,:)= par(3).*f1.^(par(3)-1).*f3.*par(1);
G(3,:)= f2.*log(f1);

