function rtp=polar(v)
r=norm(v);
co=v(3)/r;
if co > 1
co=1;
elseif co < -1
co=1;
end
theta=acos(co);
if r*sin(theta) == 0
phi=0;
elseif r*sin(theta) ~= 0
co=v(1)/(r*sin(theta));
if co > 1
co=1;
elseif co < -1
co=-1;
end
phi=acos(co);
end
rtp=[r theta phi];

