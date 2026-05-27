% cthrusw.m
%call callback routine of triplot
cthruval=get(ui1r6,'Value');
falpha=.6+(1-cthruval)*.4;
ealpha=falpha;
set(hs, 'FaceAlpha',falpha) % 0, <1;  1 for non transparent
set(hs, 'EdgeAlpha',ealpha) % 0, <1;  1 for non transparant
