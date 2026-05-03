clearvars

subj='Pig10';

[FileName,PathName,FilterIndex] = uigetfile('*.mad','Select compound results file'...
    ,['/Users/peteroosterhoff/Documents/Werk/STW/Data/results/' subj '/cluster1'],'MultiSelect','off');
[VERv,ITRIv]=loadtri(['/Users/peteroosterhoff/Documents/Werk/STW/Data/geometries/' subj '/' subj '_ventricles.tri']);

load(fullfile(PathName,FileName),'-mat','C');

cath(:,1)=cell2mat({C.Cathposx});
cath(:,2)=cell2mat({C.Cathposy});
cath(:,3)=cell2mat({C.Cathposz});

focusinit(:,1)=cell2mat({C.Initx});
focusinit(:,2)=cell2mat({C.Inity});
focusinit(:,3)=cell2mat({C.Initz});

focusinv(:,1)=cell2mat({C.Invx});
focusinv(:,2)=cell2mat({C.Invy});
focusinv(:,3)=cell2mat({C.Invz});

qtriplot({
    {'reset'}
    {VERv,ITRIv,'ventricles'}
    {'transparent ventricles=0.4'}
    {cath,[],'cath'}
    {'color cath=blue'}
    {focusinit,[],'focusinit'}
    {'color focusinit=green'}
    {focusinv,[],'focusinv'}
    {'color focusinv=red'}
    })

if 0
cent=mean(VERv,1);
qtriplot({
    {cent,[],'O'}
    {'color O=black'}
    })

for i=1:length(C)
    qtriplot({
        {[cath(i,:);focusinit(i,:)],[]}
        {'color green'}
        {'edge yes'}
        {[focusinit(i,:);focusinv(i,:)],[]}
        {'color red'}
        {'edge yes'}
        {[cent;cent-cath(i,:)+focusinit(i,:)],[]}
        {'color green'}
        {'edge on'}
        {[cent;cent-cath(i,:)+focusinv(i,:)],[]}
        {'color red'}
        {'edge on'}
        })
    
end
end
% errinit=mean(cath-focusinit,1);
% errinv=mean(cath-focusinv,1);
% qtriplot({
%     {cent,[],'O'}
%     {'color O=black'}
%     {[cent;cent+errinit],[],'errinit'}
%     {'color errinit=black'}
%     {'edge errinit=yes'}
%     });



