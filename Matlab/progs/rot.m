% rot.m
% call callback routine of triplot/ecgsim
rota=get(ui1r1,'Value');
if rota ==1,
  set(ui2a,'visible','on');
  if viscene==3, axes(ax1d),
     rotate3d on;
     axes(ax1c)
     rotate3d on;
  else, 
    axes(ax1a)
    rotate3d on; 
  end
else
  [azim,elev]= view; 
  set(ui2a,'visible','off');
  if viscene==3, axes(ax1d),
     rotate3d off;
     axes(ax1c)
     rotate3d off;
  else, 
    axes(ax1a)
    rotate3d off; 
  end
end