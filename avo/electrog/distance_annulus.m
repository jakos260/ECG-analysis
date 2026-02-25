% distance_annulus.m
% 20051105
% computes shortest distance between two ponts, pnt1 and pnt2, while traveling through an annulus
% with radii 0<a<b
function distance=distance_annulus(pnt1,pnt2,a,b)

%if norm(pnt1)<a-eps | norm(pnt1)>b+eps, 'invalid data', pause, end
%if norm(pnt2)<a-eps | norm(pnt2)>b+eps, 'invalid data', pause, end

if size(pnt1,2)==2, pnt1=[pnt1 0]; end
if size(pnt2,2)==2, pnt2=[pnt2 0]; end

if norm(pnt1-pnt2) <= eps, distance=0; return, end

% test if line intersects inner circle
[dist,lambda]=dis2line([0 0 0],pnt1,pnt2);

if dist<a && lambda<1 && lambda>0,
    %'line segment intersects inner circle' 
    % find points of contact of the tangents drawn to the inner circle;
    % first from pnt1
    
    if abs(pnt1(2))>=abs(pnt1(1)),
       argroot=max(0,(2*a^2*pnt1(1))^2 - 4 * (pnt1*pnt1') * a^2*(a^2-pnt1(2)^2));
       r11(1)=(2*a^2*pnt1(1)+sqrt(argroot))/(2*pnt1*pnt1');
       r11(2)=(a^2-r11(1)*pnt1(1))/pnt1(2);
       r12(1)=(2*a^2*pnt1(1)-sqrt((2*a^2*pnt1(1))^2 - 4 * (pnt1*pnt1') * a^2*(a^2-pnt1(2)^2)  ))/(2*pnt1*pnt1');
       r12(2)=(a^2-r12(1)*pnt1(1))/pnt1(2);
    else
       argroot=max(0,(2*a^2*pnt1(2))^2 - 4 * (pnt1*pnt1') * a^2*(a^2-pnt1(1)^2));
       r11(2)=(2*a^2*pnt1(2)+sqrt(argroot))/(2*pnt1*pnt1');
       r11(1)=(a^2-r11(2)*pnt1(2))/pnt1(1);
       r12(2)=(2*a^2*pnt1(2)-sqrt((2*a^2*pnt1(2))^2 - 4 * (pnt1*pnt1') * a^2*(a^2-pnt1(1)^2)  ))/(2*pnt1*pnt1');
       r12(1)=(a^2-r12(2)*pnt1(2))/pnt1(1);
   end
         
    P=[r11;r12];
    p1=[P(1,1:2) 0];
    p2=[P(2,1:2) 0];
    c1=cross(pnt1-p1,pnt2-p1);
    c2=cross(pnt1-p2,pnt2-p2);
    if abs(c2(3))<abs(c1(3)), p1=p2; end 
        
    % then from pnt2
    if abs(pnt2(2))>=abs(pnt2(1)),
       argroot=max(0,(2*a^2*pnt2(1))^2 - 4 * (pnt2*pnt2') * a^2*(a^2-pnt2(2)^2));
       r11(1)=(2*a^2*pnt2(1)+sqrt(argroot ))/(2*pnt2*pnt2');
       r11(2)=(a^2-r11(1)*pnt2(1))/pnt2(2);
       r12(1)=(2*a^2*pnt2(1)-sqrt((2*a^2*pnt2(1))^2 - 4 * (pnt2*pnt2') * a^2*(a^2-pnt2(2)^2)  ))/(2*pnt2*pnt2');
       r12(2)=(a^2-r12(1)*pnt2(1))/pnt2(2);
   else,
       argroot=max(0,(2*a^2*pnt2(2))^2 - 4 * (pnt2*pnt2') * a^2*(a^2-pnt2(1)^2));
       r11(2)=(2*a^2*pnt2(2)+sqrt(argroot ))/(2*pnt2*pnt2');
       r11(1)=(a^2-r11(2)*pnt2(2))/pnt2(1);
       r12(2)=(2*a^2*pnt2(2)-sqrt((2*a^2*pnt2(2))^2 - 4 * (pnt2*pnt2') * a^2*(a^2-pnt2(1)^2)  ))/(2*pnt2*pnt2');
       r12(1)=(a^2-r12(2)*pnt2(2))/pnt2(1);
   end
    
    
    P=[r11; r12;];
    q1=[P(1,1:2) 0];
    q2=[P(2,1:2) 0];
    c1=cross(pnt2-q1,pnt1-q1);
    c2=cross(pnt2-q2,pnt1-q2);
    if abs(c2(3))<abs(c1(3)), q1=q2; end 
    
    angle=acos(p1*q1'/a^2);
    distance=norm3d(pnt1-p1)+norm3d(pnt2-q1)+angle*a;
else,
    distance=norm3d(pnt1-pnt2);
end
