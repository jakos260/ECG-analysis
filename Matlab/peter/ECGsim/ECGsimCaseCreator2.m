function ECGsimCaseCreator2(GEOM,meas)

if ~exist(['.\' GEOM.subject]) ~= 7
	mkdir(['.\' GEOM.subject]);
end
outdir=[ '.\' GEOM.subject '\'];


wct=GEOM.wct;	
if ~isempty(strfind(GEOM.type,'atria'))
	saveasci([outdir GEOM.subject GEOM.beat '.adep'],meas.depfinal);
	saveasci([outdir GEOM.subject GEOM.beat '.arep'],meas.repfinal);
    if isfield(meas,'amplfinal')
        saveasci([outdir GEOM.subject GEOM.beat '.aampl'],meas.amplfinal);
    end  


	if exist( [GEOM.heartpath '_ventricles' '.tri'])
		[GEOM.VVER,GEOM.VITRI]=loadtri([GEOM.heartpath '_ventricles' '.tri']);	
	end
	[GEOM.AVER,GEOM.AITRI]=loadtri([GEOM.heartpath '_atria' '.tri']);
    
	A = GEOM.AMAORG;
	A=A-ones(length(A),1)*mean(A(wct,:)); % reference all to wct
	
	lA=length(A);
	TORSO.AMA_A=(A(1:length(GEOM.TVER),:));
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
	VENTR.AMA_A=A(tl+1:tl+length(GEOM.VVER),:);
	ATRIA.AMA_A=A(lA-length(GEOM.AVER)+1:end,:);


	savemat([outdir 'atria2atria.edl'],ATRIA.AMA_A);
	savemat([outdir 'atria2ventr.edl'],VENTR.AMA_A);
	
else
	saveasci([outdir GEOM.subject GEOM.beat '.vdep'],meas.depfinal);
	saveasci([outdir GEOM.subject GEOM.beat '.vrep'],meas.repfinal);
    if isfield(meas,'amplfinal')
        saveasci([outdir GEOM.subject beat '.vampl'],meas.amplfinal-0.85);
    end
    
	if exist( [GEOM.heartpath '_atria' '.tri'])
		[GEOM.AVER,GEOM.AITRI]=loadtri([GEOM.heartpath '_atria' '.tri']);	
	end
	[GEOM.VVER,GEOM.VITRI]=loadtri([GEOM.heartpath '_ventricles' '.tri']);

% 	savetri([outdir GEOM.subject '_ventricle' '.tri'],GEOM.VVER,GEOM.VITRI);
    V = GEOM.AMAORG;
	V=V-ones(length(V),1)*mean(V(wct,:));% reference all to wct
	lV=length(V);	
	TORSO.AMA_V=(V(1:length(GEOM.TVER),:));
	savemat([outdir 'ventr2tor.edl'],TORSO.AMA_V);

    if  size(V,1) > length(GEOM.TVER)
        tl=length(GEOM.TVER);
        VENTR.AMA_LC=V(tl+1:tl+length(GEOM.LVER),:);
        savemat([outdir 'ventr2lcav.edl'],VENTR.AMA_LC);
% 
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
        if exist( [GEOM.heartpath '_atria' '.tri'])
            ATRIA.AMA_V=V(tl+1:tl+length(GEOM.AVER),:);
            savemat([outdir 'ventr2atria.edl'],ATRIA.AMA_V);
        end
        VENTR.AMA_V=V(lV-length(GEOM.VVER)+1:end,:);
        savemat([outdir 'ventr2ventr.edl'],VENTR.AMA_V);
    end
end

%

% savetri([outdir GEOM.subject 'rho.tri'],GEOM.RVER,GEOM.RITRI);
% savetri([outdir GEOM.subject 'lho.tri'],GEOM.LVER,GEOM.LITRI);
% savetri([outdir GEOM.subject 'tor.tri'],GEOM.TVER,GEOM.TITRI);
% savetri([outdir GEOM.subject 'rlo.tri'],GEOM.RLVER,GEOM.RLITRI);
% savetri([outdir GEOM.subject 'llo.tri'],GEOM.LLVER,GEOM.LLITRI);
clear ATRIA
clear VENTR
clear TORSO
CreateLeadsystem(outdir, GEOM);
% T=intripol(GEOM.TVER,GEOM.TITRI,1:65);
GEOM.BSMALL=GEOM.BSMall;
savemat([outdir GEOM.beat '.mat'],GEOM.BSMALL);
fn=[outdir GEOM.subject '.elecs'];
VER=GEOM.TVER;%(1:65,:);
f = fopen(fn, 'w');
[nver dim]=size(VER);
fprintf(f,'%d\n ',nver);
for i=1:nver;
      fprintf(f,'%5d %8.6f %8.6f %8.6f\n',i,VER(i,1:3));
end
fclose(f);
disp('complete')


% savemat([GEOM.subject '.vrep'],GEOM.);

function CreateLeadsystem(outdir,GEOM)

% if strcmp(layfile,'nim64.mla')
%     GEOM.leadname='BSM_(nijmegen_64)';
%     GEOM.v12=[19 26 34 41 48 54 1 2 65];
% 	GEOM.wct=[1 2 65];	
% 	GEOM.barrB24=[7 10 12 4 18 24 25 26 27 31 33 34 36 34 47 46  2 45 54 60 58 57 52 61];
% 	GEOM.luxA=[61  6  5 13 15 17 19 20 22 25 26 27 28 29 30 32 34 35 36 38 41 43 44 53 54 55  2 59 58 66 63 ];% elec 137 ommitted
% 	GEOM.v12back=[19 26 34 41 48 54 1 2 27];  
% 	GEOM.finlay=[19	 35	  9	54 26 61 32  51  29   6 40  21 27  2  63  65 62  55 15 17 46  30 25   5  16 58 41 34  59 20  37 53  26 24  23] ;
% 	if strfind(GEOM.subject,'a'),	GEOM.v12=[GEOM.v12 65];end
% % 	GEOM.franklabs=[53    40      26        64        61       7        62];
%     if strcmp(GEOM.subject,'gu')
%         GEOM.franklabs=[53    40      19        66        228       7        62];%29=163
%         GEOM.v12=[19 26 65 41 48 54 1 2 66];
%         GEOM.wct=[1 2 66];	       
%     elseif strcmp(GEOM.subject,'wk6')
%         GEOM.v12=[19 26 65 41 48 54 1 2 65];
%         GEOM.franklabs=[53    40      183    65  227   174     92];
%         GEOM.bipv230=[235 121 205 280 190 26];       
%     elseif strcmp(GEOM.subject,'am')
%         GEOM.v12=[19 26 65 41 48 54 1 2 56];
%         GEOM.franklabs=[53    40      270    30  269   7     86];            
%     elseif strfind(GEOM.subject,'Haga')
%         GEOM.v12=[283 115 84 457 348 117 384 467 92];
%         GEOM.franklabs=[53    40      270    30  269   7     86];    % is nonsens
%     end
%     
%     
%     
% elseif strcmp(layfile,'ams65.mla')|| strcmp(layfile,'markams65.mla')
%     GEOM.leadname='BSM_(amsterdam_64)';
%     GEOM.v12=[12 18 25 31 40 45 63 64 65 ];
%     GEOM.v12back=[12 18 25 31 40 45 63 64 65 19];
% % 	finlayorder=[56 106 183 93 73 18 42 155 121 146 75 104 89 12 128 191 96 124 86 24 44 185 41 100 135 47 91 90 127 72 154 60 105  9 152 57    74    92   107   141];
% 	GEOM.finlay=[12	 27   9 46 18 57 24  48	 22  62 30  15 20 42  59  65 58  47  7  5 38  35 17  61   8 52 32 26  55 13  34 44  21 16  23];
% 	GEOM.luxA=[57 62 61  9  7 10 12 13 15 17 19 20 21 22 23 24 25 26 27 28 32 33 35 44 45 47 42 49 52 53 59 ];	
% 	GEOM.barrB24=[2  5  7 1 11 16 17 18 19 28 25 26 27 35 39 38 42 43 45 50 52 51 44 57];
% 	GEOM.wct=[64 63 65];
% 	GEOM.golem=1:65;GEOM.golem([12 13 14 18 19 46])=[]; 	
% % 	GEOM.franklabs=[ 45        31      19        53        10      61        58];
%     if  strcmp(GEOM.subject,'gd')
%         GEOM.franklabs=[ 212        31      19        219        16      160        98];
%     elseif strcmp(GEOM.subject,'fi')
%         GEOM.franklabs=[ 397        31      254        162        82      61        407];
%         GEOM.v12=[214 439 85 31 40 44 63 64 65 ];
%     end
%             
% end
% if strcmp(GEOM.subject,'fi')
% 	GEOM.v12=[12 18 25 32 40 45 64 63 65 ];
% end

% VFrank.name='VCG_(Frank)';
% VFrank.elecnames={ 'A' 'C' 'E'  'F'  'H'   'I' 'M'};
% VFrank.elec=GEOM.TVER(GEOM.franklabs,:);
% VFrank.lead=VFrank.elec;
% VFrank.leadNames={ 'A' 'C' 'E'  'F'  'H'   'I' 'M' 'sig 1' 'sig 2' 'sig 3'  };
% VFrank.ref=[];
% VFrank.Ref=[];
% VFrank.Transform=[];
% 
% VFrank.pos=[ [1 1]; [1 2]; [1 3]; [2 1]; [2 2]; [2 3]];
% VFrank.showLeads={'sig 1' 'sig 2' 'sig 3' 'd1', 'd2', 'd3'};
% VFrank.XVectorLead={'none' 'none' 'none' 'sig2' 'sig3' 'sig1'};

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




V9.name='minimap_montage';
V9.leadNames={'V1' 'V2' 'V3' 'V4' 'V5' 'V6' 'aVr' 'aVl' 'aVf' 'Vr' 'Vl' 'Vf'};
V9.elecnames={'V1' 'V2' 'V3' 'V4' 'V5' 'V6' 'aVr' 'aVl' 'aVf' 'Vr' 'Vl' 'Vf'};
V9.ref=[1 1 1 1 1 1 1 1 1 1 1 1];
V9.lead=GEOM.TVER([GEOM.v12 GEOM.v12(7:9)],:);
V9.elec=GEOM.TVER(GEOM.v12,:);
V9.showLeads=V9.leadNames([1:6 10:end]);
V9.XVectorLead=[];
V9.pos=[[ 2.5 2]; [3.5 2]; [4.5 2.5]; [5.5 3];[6.5 3];[7.5 3]; [1 1];[7.5 1];[7.5 5] ];
V9.Ref(1).name='extremities';
V9.Ref(1).ver=GEOM.TVER(GEOM.wct,:);
V9.Transform=[...
     0     0     0     0     0     0     1/1.5   0     0;...
     0     0     0     0     0     0     0     1/1.5   0;...
     0     0     0     0     0     0     0     0     1/1.5];
 
LAY=GEOM.LAY;
VBSM.name=GEOM.leadname;
VBSM.leadNames=cell(length(LAY)-1,1);
for i=2:length(LAY)
	VBSM.leadNames(i-1)={num2str(LAY(i,3))};
end
VBSM.Ref(1).name='extremities';
VBSM.Ref(1).ver=GEOM.TVER(GEOM.wct,:);
VBSM.elec=[GEOM.TVER(LAY(2:end,3),:) ];
VBSM.lead=[GEOM.TVER(LAY(2:end,3),:) ];
A=LAY(2:end,3);
for i=1:9
	a=cell2mat(V12.leadNames(i+3));	
	if i>6
		a=a(2:end);
    end
    if ~isempty(find(A==GEOM.v12(i)))
        VBSM.leadNames(A==GEOM.v12(i))={a};
    else
%         if strfind(layfile,'nim') && i == 3
%             VBSM.leadNames(A==33)={'V3'};
%         end
    end
end
k=4;
VBSM.elecnames=VBSM.leadNames;
VBSM.showLeads=VBSM.leadNames;
VBSM.ref=ones(length(VBSM.leadNames),1);
VBSM.pos=LAY(2:end,1:2);
VBSM.Transform=[];
VBSM.XVectorLead=[];



if isfield(GEOM,'bipv230')
    Vbip.name='exampleBipolar30mm';
    Vbip.leadNames={'1' '3' '2' '4' '5' };
    Vbip.elecnames={'1' '3' '2' '4' '5' 'ref' };
    Vbip.ref=ones(length(GEOM.bipv230)-1,1);
    Vbip.lead=GEOM.TVER(GEOM.bipv230(1:end-1),:);
    Vbip.elec=GEOM.TVER(GEOM.bipv230,:);
    Vbip.showLeads=Vbip.leadNames;
    Vbip.XVectorLead=[];
    Vbip.pos=[[ 2 3]; [3 1]; [3 3]; [2 1];[1 1];];
    Vbip.Ref(1).name='v2ref';
    Vbip.Ref(1).ver=GEOM.TVER(GEOM.bipv230(end),:);
    Vbip.Transform=[];
    saveLay([outdir GEOM.subject 'lay.bip'],Vbip);    
end


saveLay([outdir GEOM.subject 'lay.v9'],V9);
saveLay([outdir GEOM.subject 'lay.v12'],V12);
% saveLay([outdir GEOM.subject 'lay.frank'],VFrank);
saveLay([outdir GEOM.subject 'lay.bsm'],VBSM);

% /************************************************************/
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
        fprintf(f,'%8.3f \t',V.Transform(i,:));
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

		ibeg=10;
		for i=ibeg:size(V.Transform,1)
			a=cell2mat(V.leadNames(ibeg));
			fprintf(f,'%d \t %s \t 0 0 0 \t %s \n',ibeg,a,V.lead(i,:),'none');
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


function GEOM = loadGeometries(GEOM)
heartpath=GEOM.heartpath;
[GEOM.RVER,GEOM.RITRI]=loadtri([heartpath 'rho.tri']);
[GEOM.LVER,GEOM.LITRI]=loadtri([heartpath 'lho.tri']);
[GEOM.TVER,GEOM.TITRI]=loadtri([heartpath,'tor.tri']);
[GEOM.RLVER,GEOM.RLITRI]=loadtri([heartpath,'rlo.tri']);
[GEOM.LLVER,GEOM.LLITRI]=loadtri([heartpath,'llo.tri']);

