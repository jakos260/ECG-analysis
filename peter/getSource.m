function S=getSource(dep,rep)

% gets.m
% version 04/10/123
% source strength 
% based on monophasic action potential 
% dep: a fifth order polynomial with width: win(1) and 
% rep: maprep
nn=length(dep); 
[ni nj]=size(dep);
if ni<nj, dep=dep'; rep=rep'; end
t=1:max(rep)+200;nt=length(t);
if min(dep)< 5, dep=dep+5;end
dSl=-2;
rSl=0.05;
Sdep=1./(1+exp(dSl*(ones(length(dep),1)*t-dep([1:length(dep)])*ones(1,nt))));
Srep= 1./(1+exp(rSl*(ones(length(rep),1)*t-rep([1:length(rep)])*ones(1,nt))));
S=(Sdep+Srep)-1;

% figure(5);
% clf; plot(S(1,:)-1,'r') 
% hold on; plot(Sdep(1,:),'g');plot(Srep(1,:),'b');  hold off

