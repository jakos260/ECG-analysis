clearvars
leadset='all';
group=[];
saveCase =1;
sinkScan= 0;
layfile='pig64.mla';
type='ventricles';
initialVelocity=0.4;
initialActTime=15; % 15 for endo. 25 for epi?

geomdir='/Users/peteroosterhoff/Documents/Werk/STW/Data/geometries/';
subject='Pig09';
bsmdir='/Users/peteroosterhoff/Documents/Werk/STW/Data/Pig09/Biosemi/export/AVG/beats/';
%     bsmdir='/Users/peteroosterhoff/Documents/Werk/STW/Data/Pig09/Biosemi/export/beats/';
dirout= ['/Users/peteroosterhoff/Documents/Werk/STW/Data/results/' subject '/'];
heartpath=fullfile(geomdir, group,  subject,  subject );
AMAPath=[heartpath '_' type '_edl.mat']; % AMA loaded by invinit is not used!



bsmfile='/Users/peteroosterhoff/Documents/Werk/STW/Data/Pig09/Biosemi/export/AVG/beats/pig09_007_356_385_LVLatEndoThrx2SyncVoff_20130211T141537_beat000.selecg';

sigwidth=200; % number of samples to simulate

GEOM.specs(6)=0.005513;
GEOM.specs(7)=0.055512;
GEOM.specs(8)=1.317503;
pS=GEOM.specs(6:8);

velos=0.2:0.2:2;
% velos=0.4;
% vts=[0.3:0.1:0.8];
aniss=0.2:.2:2.6;
for anisi = 1:length(aniss)
    anis=aniss(anisi);
    GEOM=invInit('subject',subject,'layfile',layfile,'bsm',bsmfile,'type',type,'anisotropyRatio',anis,'group',group,'basedir',geomdir);
    
    %     A=loadmat(AMAPath);
    %     GEOM.AMA=zeromean( 40*A(1:length(GEOM.TVER),:)); % override A matrix
    useLeads=1:64;
    GEOM = selectLeads(GEOM,useLeads,1);
    
    if ~exist('fwsig','var')
        fwsig=zeros(length(aniss),length(velos),length(useLeads),sigwidth);
    end
    for veloi = 1:length(velos)
        %     for vti = 1:length(vts)
        %         vt=vts(vti);
        %         anis = velo/vt;
        velo = velos(veloi);
        fprintf('anisotropy %0.1f, velocity %0.1fm/s. %d of %d\n',anis,velo,(anisi-1)*length(velos)+veloi,length(aniss)*length(velos));
        vt=velo/anis;
        inode=1791; % Pig09-refined 1791; LVLateralndo, 1734 LVLateralEpi
        dep=GEOM.DIST2W(:,inode)/initialVelocity;
        
%         dep=min(GEOM.DIST2W(:,1791),GEOM.DIST2W(:,1556))/initialVelocity; % Endo midden in de wand
        depAbove=find(dep>initialActTime);
        dep(depAbove) = initialActTime + (dep(depAbove) - initialActTime)*initialVelocity / velo;

        
        fwsig(anisi,veloi,:,:) =GEOM.AMA*getSmode(ones(length(GEOM.VER),1) * (1:sigwidth),dep,300*ones(size(GEOM.VER,1),1),pS,[],1);
        
        
    end
end

% save('/Users/peteroosterhoff/Documents/Werk/STW/Data/Pig09/FWSIGPig09LVatEndo20130606.mat','GEOM','aniss','velos','fwsig');
% save('/Users/peteroosterhoff/Documents/Werk/STW/Data/Pig09/FWSIGPig09LVatEpi20130606.mat','GEOM','aniss','velos','fwsig');
% save('/Users/peteroosterhoff/Documents/Werk/STW/Data/Pig09/FWSIGPig09LVatEndoInWall20130610.mat','GEOM','aniss','velos','fwsig');
% save('/Users/peteroosterhoff/Documents/Werk/STW/Data/Pig09/FWSIGPig09LVatEndoInitVel0.4_1.25_20130611.mat','GEOM','aniss','velos','fwsig');