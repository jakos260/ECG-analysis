function selectedLeads=multifociscanLeadselect(GEOM,refFoci,refDep,usecor)
global lpass

qrsduration=GEOM.specs(3)-GEOM.specs(2)+1;
cors=zeros(size(1:65,1));
rds=zeros(size(1:65,1));
t=0:qrsduration-1;
T=ones(length(GEOM.VER),1)*t;
leads=1:65;
for j=1:length(leads)
	tleads=leads;
	tleads(j)=[];
	PSIREF=baselinecor(zeromean(GEOM.BSM(tleads,GEOM.specs(2):GEOM.specs(5))));
	PSIREF=PSIREF(:,1:GEOM.specs(3)-GEOM.specs(2)+1);
	AMA=zeromean(GEOM.AMA(tleads,:));
	PSIA=lowpassma(AMA*getSmode(T,refDep,[],GEOM.pS,[],1),lpass);
	rds(j)=norm(PSIA-PSIREF,'fro')/norm(PSIREF,'fro');
	cor=corrcoef(PSIA,PSIREF);
	cors(j)=cor(2,1);
end
if ~usecor
	A=[(1:65)' rds'];A=sortrows(A,2);			
	selectedLeads=A(end-19:end,1);
else
	A=[(1:65)' cors'];A=sortrows(A,2);			
	selectedLeads=A(1:20,1);
end
