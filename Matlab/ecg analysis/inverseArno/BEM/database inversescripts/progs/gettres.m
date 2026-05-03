% gettres.m; script of invedl.m
% compute regtest, restest and trestest
% for times: testtims
% 2013-05-20

if strcmp(pol,'dep') dep=testtims; end
if strcmp(pol,'rep') rep=testtims; end
S=gets(T,dep,rep,p,mode);
S=S(:,1:usetime);
PHIA=A*S;
RES=PHI-PHIA;
rdtest=norm(RES,'fro')/norm(PHI(:,1:usetime),'fro')   % NOTE: unfiltered rd
%regtest=norm(REGOP*(rep-dep))/1000; % note: rep-dep=ari

regtest=norm(REGOP*(rep-dep))/1000; % note: rep-dep=ari

trestest=sqrt(rdtest^2+(regtest*mu)^2);
%rrms=std(PHIA);
%apex=max(rrms(100:usetime));
