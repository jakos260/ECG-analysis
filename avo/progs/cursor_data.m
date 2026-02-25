% file cursor_data.m
% date: 2014_04_12
% identifies data point i closest to cursor
% function [i, x y ]=cursor_data(fx,fy,fig)

function [i x y ]=cursor_data(fx,fy,fig)

figure(fig)
if size(fx,2)>size(fx,1), fx=fx'; end
if size(fy,2)>size(fy,1), fy=fy'; end
  

if size(fx,1)~=size(fy,1),
    'length fx and fy should be equal'
    pause
    return
end

ndat=size(fx,1);

DATA=[fx(:,1) fy(:,1)];
   [x y]=ginput(1);
   % identify node
DIST=DATA-ones(ndat,1)*[x y];

% identify node
DIST(:,1)=DIST(:,1)/(max(DIST(:,1))-min(DIST(:,1)));    
DIST(:,2)=DIST(:,2)/(max(DIST(:,2))-min(DIST(:,2)));  
d=DIST(:,1).^2+DIST(:,2).^2;
[~,i]=min(d);
disp(num2str([x y i]))
