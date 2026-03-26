% gettres.m; script of invedl.m
% computes regtest, restest and trestest
% based on  S(testtims)
% note: initial scaling by a factor 1000 is replaced by  n_S: the number of
%       nodes on Sh; scaling used to bring level up the apples and pears in trestest 
% 2015-02-09

if strcmp(pol,'dep') dep=testtims; end
if strcmp(pol,'rep') rep=testtims; end

S=gets(T,dep,rep,p,mode);
S=S(:,1:usetime);
PHIA=A*S;
RES=PHI-PHIA;
rdtest=norm(RES,'fro')/norm(PHI(:,1:usetime),'fro');   % NOTE: unfiltered rd
n_S=size(S,1);

if strcmp(pol,'rep'),
    if strcmp(rep_regtype(1),'a'),
        regtest=norm(REGOP*(rep-dep));
        %regtest=norm(REGOP*(rep-dep)/n_S); % note: rep-dep=ari
        %regtest=std(REGOP*(rep-dep));      % note: rep-dep=ari
        trestest=sqrt(rdtest^2+(regtest*muari)^2);
    else,
        %regtest=norm(REGOP*rep/n_S);
        regtest=norm(REGOP*rep);
        %regtest=std(REGOP*rep);
        trestest=sqrt(rdtest^2+(regtest*murep)^2);
    end   
else,
    %regtest=norm(REGOP*rep/n_S);
    regtest=std(REGOP*dep);
    trestest=sqrt(rdtest^2+(regtest*mudep)^2);
end


