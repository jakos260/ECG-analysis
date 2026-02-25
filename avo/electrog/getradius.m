% getradius
% script of electrog
% 20030130
% accept radius inner circle;
% and recompute geometry
rara=get(uie1,'string');
a=str2num(rara);
if r8==1, robs=a; end
setgeom
gettiming
getphi
if r6==1, plotphi,
else,
plotfront
end
setuis  