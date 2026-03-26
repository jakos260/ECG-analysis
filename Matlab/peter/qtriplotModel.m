function qtriplotModel(dirIn, subject,dirout,varargin)


if ~exist(dirout)
    mkdir(dirout)
end
leadLoc = [];
if nargin >= 5
    leadLoc = varargin{5};
end;

outname = fullfile(dirout,[subject 'openthis.trp']);

vname = [subject '_ventricles.tri'];
vtname = [subject '_ventricles.typ'];
aname = [subject '_atria.tri'];
atname = [subject '_atria.typ'];
leadLocName= [subject '_leadloc.tri'];
ladname = [subject '_lad.tri'];
rcaname = [subject '_rca.tri'];
torname = [subject '_thorax.tri'];
v12name = [subject '_v12.tri'];

[GEOM.VER,GEOM.ITRI]=loadtri(fullfile(dirIn,subject, vname) );
[GEOM.VERA,GEOM.ITRIA]=loadtri(fullfile(dirIn,subject, aname) );
[GEOM.LADVER,GEOM.LADITRI]=loadtri(fullfile(dirIn,subject,ladname));
[GEOM.RCAVER,GEOM.RCAITRI]=loadtri(fullfile(dirIn,subject,rcaname));
[GEOM.TVER,GEOM.TITRI]=loadtri(fullfile(dirIn,subject,torname));
GEOM.typ = loadmat(fullfile(dirIn,subject,vtname) );
GEOM.typa = loadmat(fullfile(dirIn,subject,atname) );
GEOM.ADJ = graphdist(GEOM.ITRI,GEOM.VER,1);


savetri(fullfile(dirout,vname)  ,GEOM.VER,GEOM.ITRI);
savetri(fullfile(dirout,ladname),GEOM.LADVER,GEOM.LADITRI);
savetri(fullfile(dirout,rcaname),GEOM.RCAVER,GEOM.RCAITRI);

savetri(fullfile(dirout,v12name),GEOM.TVER(4:9,:),[]);

% vdepname = [subject '.dep'];
% saveasci(fullfile(dirout,vdepname),meas.depfinal- min(meas.depfinal));

aortaname = [subject '_aorta_mitral.tri'];
aortaAname = [subject '_aortaA.tri'];
tricusname = [subject '_tricus.tri'];
pulmname = [subject '_pulm.tri'];

a = find(GEOM.typ==4 | GEOM.typ== 7);
edges=a(1);
a(1)=[];
while ~isempty(a)
    ai= find(GEOM.ADJ(edges(end),a) > 0);
    ai=ai(GEOM.ADJ(edges(end),a(ai)) == min(GEOM.ADJ(edges(end),a(ai))));
    edges=[edges a(ai)];
    a(ai)=[];
end
edges=[edges edges(1)];
savetri(fullfile(dirout,aortaname),GEOM.VER(edges,:),[]);







a = find(GEOM.typ==5 );
edges=a(1);
a(1)=[];
while ~isempty(a)
    ai= find(GEOM.ADJ(edges(end),a) > 0);
    ai=ai(GEOM.ADJ(edges(end),a(ai)) == min(GEOM.ADJ(edges(end),a(ai))));
    edges=[edges a(ai)];
    a(ai)=[];
end
edges=[edges edges(1)];
savetri(fullfile(dirout,tricusname),GEOM.VER(edges,:),[]);

a = find(GEOM.typ==6 );
edges=a(1);
a(1)=[];
while ~isempty(a)
    ai= find(GEOM.ADJ(edges(end),a) > 0);
    ai=ai(GEOM.ADJ(edges(end),a(ai)) == min(GEOM.ADJ(edges(end),a(ai))));
    edges=[edges a(ai)];
    a(ai)=[];
end
edges=[edges edges(1)];
savetri(fullfile(dirout,pulmname),GEOM.VER(edges,:),[]); 

a = find(GEOM.typa==7 | GEOM.typa==4 );
ADJ=graphdist(GEOM.ITRIA,GEOM.VERA,1);
edges=a(1);
a(1)=[];
while ~isempty(a)
    ai= find(ADJ(edges(end),a) > 0);
    ai=ai(ADJ(edges(end),a(ai)) == min(ADJ(edges(end),a(ai))));
    edges=[edges a(ai)];
    a(ai)=[];
end
edges=[edges edges(1)];
savetri(fullfile(dirout,aortaAname),GEOM.VERA(edges,:),[]);

fid=fopen(outname,'wt');
fprintf(fid,'%s\n','bgdcolor white');
fprintf(fid,'%s\n','bgdtrans no');
% fprintf(fid,'%s\n',['file ' vlendoname]);
fprintf(fid,'%s\n',['file ' vname]);
fprintf(fid,'%s\n','cross y');
% fprintf(fid,'%s\n','back both ');
fprintf(fid,'%s\n','color heart');
% fprintf(fid,'%s\n',['funfile ' vdepname]);
% fprintf(fid,'%s\n','funcolor tims ');

fprintf(fid,'%s\n','step 10');
fprintf(fid,'%s\n',['file ' aortaname]);
fprintf(fid,'%s\n','edge y');
% fprintf(fid,'%s\n','color darkgrey');
fprintf(fid,'%s\n','bold 6');
fprintf(fid,'%s\n',['file ' aortaAname]);
fprintf(fid,'%s\n','edge y');
fprintf(fid,'%s\n','bold 6');
fprintf(fid,'%s\n',['file ' tricusname]);
fprintf(fid,'%s\n','edge y');
fprintf(fid,'%s\n','bold 6');
fprintf(fid,'%s\n',['file ' pulmname]);
fprintf(fid,'%s\n','edge y');
fprintf(fid,'%s\n','bold 6');
fprintf(fid,'%s\n',['file ' ladname]);
fprintf(fid,'%s\n','color red');
fprintf(fid,'%s\n',['file ' rcaname]);
fprintf(fid,'%s\n','color red');
% fprintf(fid,'%s\n',['file ' v12name]);
% fprintf(fid,'%s\n','color black');
fprintf(fid,'%s\n','mouse fun');
if ~isempty(leadLoc)
    fprintf(fid,'%s\n',['file ' leadLocName]);
    fprintf(fid,'%s\n','color red');

end
fclose(fid);

% qtriplot('delete *');
% qtriplot(outname);