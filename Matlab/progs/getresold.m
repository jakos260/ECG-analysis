% getres.m; 
% function rd=getres(shift,scal,A,PHI,maprep,ttm,depin,repin,win,pol,nodes,t)
% used to optimize scaling and shifting of dep
function   rd=getres(shift,scal,A,PHI,maprep,ttm,depin,repin,win,pol,nodes,t)
phinorm=norm(PHI,'fro');
dep=shift+scal*depin;
rep=repin;
gets;
PHIA=A*S;
RES=PHI-PHIA;
rd=norm(RES,'fro')/phinorm;
