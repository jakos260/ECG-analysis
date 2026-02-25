% gridsw.m
% call callback routine of triplot/ ecgsim
% grsw=get(ui1r2,'Value');
  grsw=get(ui1r2,'Value');
if grsw ==1,
  if exist('hs3') set(hs3,'EdgeColor',[0.5 0.5 0.5]); end
    if exist('hs3'),
      viscene==3, axes(ax1c),set(hs3,'EdgeColor',[0.5 0.5 0.5]);
  end
  set(hs,'EdgeColor',[0.5 0.5 0.5]);
  if exist('hbspm')  set(hbspm, 'EdgeColor',[0.5 0.5 0.5]); end
  if exist('hbspm3') set(hbspm3,'EdgeColor',[0.5 0.5 0.5]); end
else, 
  if exist('hs3') set(hs3,'EdgeColor','none'); end
  if exist('viscine'),
     if viscene==3, axes(ax1c),set(hs3,'EdgeColor','none'); end
  end
  set(hs,'EdgeColor','none');
  if exist('hbspm')   set(hbspm,'EdgeColor','none');end
  if exist('hbspm3') set(hbspm3,'EdgeColor','none');end
end