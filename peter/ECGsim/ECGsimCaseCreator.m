function ECGsimCaseCreator

global outdir
outdir='C:\PEACS\project\ECGsim\code\ECGsim2.0\casedata\';

dirout='C:\ECG_simulation\inverse\matlab\results\';
GEOM=[];
beat='01';
subj='wk0';GEOM.type='atria';

subj='fi';beat='04';GEOM.type=''; bsmfile='fi04.mes';

heartpath=['C:\ECG_simulation\geometries\' subj  '\' subj ];
alpha=2;
leadse='';


GEOM.heartpath=heartpath;
GEOM.subject=subj;
[GEOM.RVER,GEOM.RITRI]=loadtri([heartpath 'rho.tri']);
[GEOM.LVER,GEOM.LITRI]=loadtri([heartpath 'lho.tri']);
[GEOM.TVER,GEOM.TITRI]=loadtri([heartpath,'tor.tri']);
[GEOM.RLVER,GEOM.RLITRI]=loadtri([heartpath,'rlo.tri']);
[GEOM.LLVER,GEOM.LLITRI]=loadtri([heartpath,'llo.tri']);





load([dirout GEOM.subject '_' beat '_alpha' num2str(alpha) leadse GEOM.type '_deprep.mat'])


if strcmp(GEOM.type,'atria')
	saveasci([outdir GEOM.subject '.adep',meas.depfinal]);
	saveasci([outdir GEOM.subject '.arep',meas.repfinal]);
	[GEOM.AVER,GEOM.AITRI]=loadtri([heartpath '_atria' '.tri']);	
	savetri([outdir '_atria' '.tri'],GEOM.AVER,GEOM.AITRI);		
	A=loadmat([heartpath GEOM.subject  '_atria_edl.mat']);
	A=40*zeromean(A);
	lA=length(A);
	TORSO.AMA_A=A(1:length(GEOM.TVER),:);
	savemat([outdir 'atria2tor.edl'],TORSO.AMA_A);

	tl=length(GEOM.TVER);
	ATRIA.AMA_LC=A(tl+1:tl+length(GEOM.LVER),:);
	savemat([outdir 'atria2lcav.edl'],ATRIA.AMA_LC);
	tl=tl+length(GEOM.LVER);
	
	ATRIA.AMA_RC=A(tl+1:tl+length(GEOM.RVER),:);
	savemat([outdir 'atria2rcav.edl'],ATRIA.AMA_RC);
	
	tl=tl+length(GEOM.RVER);
	TORSO.AMA_A_RL=A(tl+1:tl+length(GEOM.RLVER),:);
	savemat([outdir 'atria2rlung.edl'],TORSO.AMA_A_RL);

	tl=tl+length(GEOM.RLVER);
	TORSO.AMA_A_LL=A(tl+1:tl+length(GEOM.LLVER),:);
	savemat([outdir 'atria2llung.edl'],TORSO.AMA_A_LL);


	tl=tl+length(GEOM.LLVER);
	VENTR.AMA_A=A(tl+1:tl+length(GEOM.AVER),:);
	ATRIA.AMA_A=A(lA-length(GEOM.AVER)+1:end,:);


	savemat([outdir 'atria2atria.edl'],ATRIA.AMA_A);
	savemat([outdir 'ventr2atria.edl'],VENTR.AMA_A);
	
else
	saveasci([outdir GEOM.subject '.vdep'],meas.depfinal);
	saveasci([outdir GEOM.subject '.vrep'],meas.repfinal);
	[GEOM.VVER,GEOM.VITRI]=loadtri([heartpath '_ventricle' '.tri']);
	savetri([outdir GEOM.subject '_ventricle' '.tri'],GEOM.VVER,GEOM.VITRI);
	V=loadmat([heartpath '_ventricle_edl.mat']);
	V=40*zeromean(V);
	lV=length(V);	
	TORSO.AMA_V=V(1:length(GEOM.TVER),:);
	savemat([outdir 'ventr2tor.edl'],TORSO.AMA_V);

	tl=length(GEOM.TVER);
	VENTR.AMA_LC=V(tl+1:tl+length(GEOM.LVER),:);
	savemat([outdir 'ventr2lcav.edl'],VENTR.AMA_LC);

	tl=tl+length(GEOM.LVER);
	VENTR.AMA_RC=V(tl+1:tl+length(GEOM.RVER),:);
	savemat([outdir 'ventr2rcav.edl'],VENTR.AMA_RC);

	tl=tl+length(GEOM.RVER);
	TORSO.AMA_V_RL=V(tl+1:tl+length(GEOM.RLVER),:);	
	savemat([outdir 'ventr2rlung.edl'],TORSO.AMA_V_RL);

	tl=tl+length(GEOM.RLVER);
	TORSO.AMA_V_LL=V(tl+1:tl+length(GEOM.LLVER),:);
	savemat([outdir 'ventr2llung.edl'],TORSO.AMA_V_LL);

	tl=tl+length(GEOM.LLVER);
% 	ATRIA.AMA_V=V(tl+1:tl+length(GEOM.VVER),:);
	VENTR.AMA_V=V(lV-length(GEOM.VVER)+1:end,:);
	
% 	//savemat([outdir 'atria2ventr.edl'],ATRIA.AMA_V);
	savemat([outdir 'ventr2ventr.edl'],VENTR.AMA_V);


end

%
CreateLeaysystem(GEOM,'ams65.mla');
savetri([outdir GEOM.subject 'rho.tri'],GEOM.RVER,GEOM.RITRI);
savetri([outdir GEOM.subject 'lho.tri'],GEOM.LVER,GEOM.LITRI);
savetri([outdir GEOM.subject 'tor.tri'],GEOM.TVER,GEOM.TITRI);
savetri([outdir GEOM.subject 'rlo.tri'],GEOM.RLVER,GEOM.RLITRI);
savetri([outdir GEOM.subject 'llo.tri'],GEOM.LLVER,GEOM.LLITRI);


GEOM.BSM=loadmat(['C:\ECG_simulation\geometries\' GEOM.subject '\' bsmfile]);
savemat([outdir bsmfile],GEOM.BSM);
fn=[outdir GEOM.subject '.elecs'];
VER=GEOM.TVER(1:65,:);
f = fopen(fn, 'w');
[nver dim]=size(VER);
fprintf(f,'%d\n ',nver);
for i=1:nver;
      fprintf(f,'%5d %8.6f %8.6f %8.6f\n',i,VER(i,1:3));
end
fclose(f);
disp('complete')


% savemat([GEOM.subject '.vrep'],GEOM.);

function CreateLeaysystem(GEOM,layfile)
global outdir

if strcmp(layfile,'nim64.mla')
    GEOM.leadname='nijmegen';
    GEOM.v12=[19 26 34 41 48 54 1 2];
	GEOM.wct=[1 2 65];	
	GEOM.barrB24=[7 10 12 4 18 24 25 26 27 31 33 34 36 34 47 46  2 45 54 60 58 57 52 61];
	GEOM.luxA=[61  6  5 13 15 17 19 20 22 25 26 27 28 29 30 32 34 35 36 38 41 43 44 53 54 55  2 59 58 66 63 ];% elec 137 ommitted
	GEOM.v12back=[19 26 34 41 48 54 1 2 27];  
	GEOM.finlay=[19	 35	  9	54 26 61 32  51  29   6 40  21 27  2  63  65 62  55 15 17 46  30 25   5  16 58 41 34  59 20  37 53  26 24  23] ;
	if strfind(GEOM.subject,'a'),	GEOM.v12=[GEOM.v12 65];end
	GEOM.franklabs=[53    33      26        64        61       7        62];
elseif strcmp(layfile,'ams65.mla')|| strcmp(layfile,'markams65.mla')
    GEOM.leadname='amsterdam';
    GEOM.v12=[12 18 25 31 40 45 63 64 65 ];
    GEOM.v12back=[12 18 25 31 40 45 63 64 65 19];
% 	finlayorder=[56 106 183 93 73 18 42 155 121 146 75 104 89 12 128 191 96 124 86 24 44 185 41 100 135 47 91 90 127 72 154 60 105  9 152 57    74    92   107   141];
	GEOM.finlay=[12	 27   9 46 18 57 24  48	 22  62 30  15 20 42  59  65 58  47  7  5 38  35 17  61   8 52 32 26  55 13  34 44  21 16  23];
	GEOM.luxA=[57 62 61  9  7 10 12 13 15 17 19 20 21 22 23 24 25 26 27 28 32 33 35 44 45 47 42 49 52 53 59 ];	
	GEOM.barrB24=[2  5  7 1 11 16 17 18 19 28 25 26 27 35 39 38 42 43 45 50 52 51 44 57];
	GEOM.wct=[64 63 65];
	GEOM.golem=1:65;GEOM.golem([12 13 14 18 19 46])=[]; 	
	GEOM.franklabs=[ 45        31      19        53        10      61        58];
end
if strcmp(GEOM.subject,'fi')
	GEOM.v12=[12 18 25 32 40 45 64 63 65 ];
end

VFrank.name='Frankleads';
VFrank.elecnames={ 'A' 'C' 'E'  'F'  'H'   'I' 'M'};
VFrank.elec=GEOM.TVER(GEOM.franklabs,:);
VFrank.lead=VFrank.elec;
VFrank.leadNames={ 'A' 'C' 'E'  'F'  'H'   'I' 'M' 'sig 1' 'sig 2' 'sig 3'  };
VFrank.ref=[];
VFrank.Ref=[];
VFrank.Transform=[];

VFrank.pos=[ [1 1]; [1 2]; [1 3]; [2 1]; [2 2]; [2 3]];
VFrank.showLeads={'sig 1' 'sig 2' 'sig 3' 'd1', 'd2', 'd3'};
VFrank.XVectorLead={'none' 'none' 'none' 'sig2' 'sig3' 'sig1'};

V12.name='standard_12';
V12.leadNames={'I' 'II' 'III' 'V1' 'V2' 'V3' 'V4' 'V5' 'V6' 'aVr' 'aVl' 'aVf' };
V12.elecnames={'V1' 'V2' 'V3' 'V4' 'V5' 'V6' 'Vr' 'Vl' 'Vf' };
V12.ref=[2 2 3 1 1 1 1 1 1 1 1 1];
V12.lead=GEOM.TVER([GEOM.wct(2) GEOM.wct(3) GEOM.wct(3) GEOM.v12],:);
V12.elec=GEOM.TVER(GEOM.v12,:);
V12.showLeads=V12.leadNames;
V12.XVectorLead=[];
V12.pos=[[ 1 1]; [1 2]; [1 3]; [3 1];[3 2];[3 3];[4 1];[4 2];[4 3];[2 1];[2 2];[2 3]; ];
V12.Ref(1).name='extremities';
V12.Ref(1).ver=GEOM.TVER(GEOM.wct,:);
V12.Ref(2).name='vr';
V12.Ref(2).ver=GEOM.TVER(GEOM.wct(1),:);
V12.Ref(3).name='vl';
V12.Ref(3).ver=GEOM.TVER(GEOM.wct(2),:);
V12.Transform=[];
LAY=loadmat(layfile);


VBSM.name=GEOM.leadname;
VBSM.leadNames=cell(length(LAY)-1,1);
for i=2:length(LAY)
	VBSM.leadNames(i-1)={num2str(LAY(i,3))};
end
VBSM.Ref(1).name='extremities';
VBSM.Ref(1).ver=GEOM.TVER(GEOM.wct,:);
VBSM.elec=[GEOM.TVER(1:65,:) ];
VBSM.lead=[GEOM.TVER(LAY(2:end,3),:) ];
A=LAY(2:end,3);
for i=1:9
	a=cell2mat(V12.leadNames(i+3));	
	if i>6
		a=a(2:end);
	end
	VBSM.leadNames(A==GEOM.v12(i))={a};
end
k=4;
VBSM.elecnames=VBSM.leadNames;
VBSM.showLeads=VBSM.leadNames;
VBSM.ref=ones(65,1);
VBSM.pos=LAY(2:end,1:2);
VBSM.Transform=[];
VBSM.XVectorLead=[];

saveLay([outdir GEOM.subject 'lay.v12'],V12);
saveLay([outdir GEOM.subject 'lay.frank'],VFrank);
saveLay([outdir GEOM.subject 'lay.bsm'],VBSM);


function saveLay(fn,V)

f=fopen(fn,'w');
if (f>0)

	fprintf(f,'%s\n',V.name);
	fprintf(f,'%d\n',size(V.elec,1));
	for i=1:size(V.elec,1);
		fprintf(f,'%d \t %8.3f %8.3f %8.3f \n',i,V.elec(i,:));
	end
	fprintf(f,'References %d\n',length(V.Ref));
	for i=1:length(V.Ref)
		fprintf(f,'%s\n',V.Ref(i).name);
		fprintf(f,'%d\n ',size(V.Ref(i).ver,1));
		for j=1:size(V.Ref(i).ver,1)
			fprintf(f,'%8.3f %8.3f %8.3f\n',V.Ref(i).ver(j,:));
		end
	end
	fprintf(f,'Transform %d %d\n',size(V.Transform,1),size(V.Transform,2));
	for i=1:size(V.Transform,1)
		for j=1:size(V.Transform,2)
			fprintf(f,'%8.3f /t',V.Transfrom);
		end
		fprintf(f,'\n');
	end
	fprintf(f,'Leads %d\n',length(V.lead));
	for i=1:length(V.lead);
		a=cell2mat(V.leadNames(i));
		if isempty(V.Ref)
			fprintf(f,'%d \t %s \t %8.3f %8.3f %8.3f \t %s \n',i,a,V.lead(i,:),'none');
		else
			fprintf(f,'%d \t %s \t %8.3f %8.3f %8.3f \t %s \n',i,a,V.lead(i,:),V.Ref(V.ref(i)).name);
		end
	end
	if isempty(strfind(V.name,'Frank'))

		ibeg=i;
		for i=1:size(V.Transform,1)
			a=cell2mat(V.leadNames(i+ibeg));
			fprintf(f,'%d \t %s \t 0 0 0 \t %s \n',i+ibeg,a,V.lead(i,:),'none');
		end

		fprintf(f,'LeadsPos %d\n',length(V.showLeads));
		for i=1:length(V.showLeads);
			a=cell2mat(V.showLeads(i));
			if isempty(V.XVectorLead)
				b='none';
			else
				b=cell2mat(V.XVectorLead(i));
			end
			fprintf(f,'%d \t %s \t %s\t %8.3f %8.3f  \n',i,a,b,V.pos(i,:));
		end
	end
end
fclose(f);