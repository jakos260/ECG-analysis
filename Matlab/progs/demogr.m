% program demogr

% step=10 followed by 
% steptime
% advances time 10 years
% monitor shows time course of total population
clear
echo off
clf

filestm='ned81.stm'
filestf='ned81.stv'
filepmi='ned81.pmi'
filepfi='ned81.pvi'
filevrt='ned81.vrt'
fracmale=0.513;

pmi=loadasci(filepmi);
pfi=loadasci(filepfi);
vrt=loadasci(filevrt);
vrt(1:101,2)=0.001*vrt(1:101,2);
stm=loadasci(filestm);
stm(1:101,2)=0.001*stm(1:101,2);
stf=loadasci(filestf);
stf(1:101,2)=0.001*stf(1:101,2);
[nrow,ncol]=size(pmi);
year=pmi(nrow,2);
itim=1;
survey(itim,1)=year;
t=0:nrow-2;
ptot=pfi+pmi;
born=vrt(1:101,2).*pfi(1:101,2);
aanwas=sum(born)
mdeath=stm(1:101,2).*pmi(1:101,2);
fdeath=stf(1:101,2).*pfi(1:101,2);
exitmale=sum(mdeath)
exitfemale=sum(fdeath)
survey(itim,2)=sum(ptot(1:101,2));
survey(itim,3)=aanwas-exitmale-exitfemale;
survey(itim,4)=sum(ptot(21:66,2)/survey(itim,2));
aanwas-exitmale-exitfemale;

fourpl
step=10
saveasci('survey.lst',survey);
'specify step; followed by command: steptime'

