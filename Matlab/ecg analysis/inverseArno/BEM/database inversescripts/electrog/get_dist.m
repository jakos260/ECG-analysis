% get_dist.m
% script of electrog
% 20051105
% computes distance between two ponts while traveling through an annulus
clear

a=0.7;% inner circle
b=1;% outer circle
% plot annulus
for i=1:101,
    phi=(i-1)/100*2*pi;
    ya(i)=a*sin(phi);
    yb(i)=b*sin(phi);
    xa(i)=a*cos(phi);
    xb(i)=b*cos(phi);
end
clf
plot(xa,ya)
hold on
plot(xb,yb)
axis square

pnt1=[0 1   0];

if norm(pnt1)<a | norm(pnt1)>b, 'invalid data', stop, end

pnt2=[0 -.9  0];
if norm(pnt2) < a | norm(pnt2) > b, 'invalid data', stop,end

plot(pnt1(1),pnt1(2),'r*');
plot(pnt2(1),pnt2(2),'r*');
line([pnt1(1) pnt2(1)],[pnt1(2) pnt2(2)],'col','r')

% test if line intersects inner circle
dist=dis2line([0 0 0],pnt1,pnt2);

if dist<a,
    'line intersects inner circle' 
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
   
    P=[r11;r12]
    p1=[P(1,1:2) 0];
    p2=[P(2,1:2) 0];
    c1=cross(pnt1-p1,pnt2-p1);
    c2=cross(pnt1-p2,pnt2-p2);
    if abs(c2(3))<abs(c1(3)), p1=p2; end 
        plot(p1(1),p1(2),'g*')
        plot([pnt1(1) p1(1)],[pnt1(2) p1(2)],'g')  
        
    % then from pnt2
    if abs(pnt2(2))>=abs(pnt2(1)),
       argroot=max(0,(2*a^2*pnt2(1))^2 - 4 * (pnt2*pnt2') * a^2*(a^2-pnt2(2)^2))
       r11(1)=(2*a^2*pnt2(1)+sqrt(argroot ))/(2*pnt2*pnt2');
       r11(2)=(a^2-r11(1)*pnt2(1))/pnt2(2);
       r12(1)=(2*a^2*pnt2(1)-sqrt((2*a^2*pnt2(1))^2 - 4 * (pnt2*pnt2') * a^2*(a^2-pnt2(2)^2)  ))/(2*pnt2*pnt2');
       r12(2)=(a^2-r12(1)*pnt2(1))/pnt2(2);
   else,
       argroot=max(0,(2*a^2*pnt2(2))^2 - 4 * (pnt2*pnt2') * a^2*(a^2-pnt2(1)^2))
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
        plot(q1(1),q1(2),'g*')
        plot([pnt2(1) q1(1)],[pnt2(2) q1(2)],'g')  
    angle=acos(p1*q1'/a^2)
    distance=norm3d(pnt1-p1)+norm3d(pnt2-q1)+angle*a
        dist=norm3d(pnt1-pnt2)
    else,
        distance=norm3d(pnt1-pnt2)
end
