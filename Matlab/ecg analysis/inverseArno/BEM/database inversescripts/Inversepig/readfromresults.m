clearvars
Rinit={};
Rinv={};
fileinit={};
fileinv={};

% basepath='/Users/peteroosterhoff/Documents/MATLAB/InversePig/results/Pig07/v0/single/CORxcorrrmsshiftNewBSMNewGeom';
% basepath='/Users/peteroosterhoff/Documents/MATLAB/InversePig/results/Pig07/v0/single/CORNoShiftNewBSMNewGeom';
% basepath='/Users/peteroosterhoff/Documents/MATLAB/InversePig/results/Pig07/v0/single';% controle PvD Peters Oude Beats
% basepath='/Users/peteroosterhoff/Documents/MATLAB/InversePig/results/Pig07/v0/single/NewPeter20121022'; % opnieuw controle PvD. Eigen code useleads gefixed.
% basepath='/Users/peteroosterhoff/Documents/MATLAB/InversePig/results/Pig07/v0/single/PO_OptimRatio1.0_20121022';
% basepath='/Users/peteroosterhoff/Documents/MATLAB/InversePig/results/Pig07/v0/single/PvDHomogeenVarkenPig07_1.0_Vl1.5';
% basepath='/Users/peteroosterhoff/Documents/MATLAB/InversePig/results/Pig07/v0/single/homogeneous';
% basepath='/Users/peteroosterhoff/Documents/Werk/STW/Data/results/Pig07/v0/single/NotInWallmudep1e-4murep1e-4';
% basepath='/Users/peteroosterhoff/Documents/Werk/STW/Data/results/Pig07/v0/single/PO_Ratio1.0QRSScaleShiftInit10';
% basepath='/Users/peteroosterhoff/Documents/Werk/STW/Data/results/Pig08/single/20130214';
% basepath='/Users/peteroosterhoff/Documents/Werk/STW/Data/results/Pig09/single/20130214';
% basepath='/Users/peteroosterhoff/Documents/Werk/STW/Data/results/Pig07/v0/single/Ingetest';
% basepath='/Users/peteroosterhoff/Documents/Werk/STW/Data/results/Pig09/single/Single_20130409_mudep1.5e-4_murep1.5e-9_FixedRepol_NotInWall';


% [VERv ITRIv]=loadtri('/Users/peteroosterhoff/Documents/Werk/STW/Data/geometries/Pig07/Ingetest/Pig07/Pig07_ventricles.tri');
% [VERv ITRIv]=loadtri('/Users/peteroosterhoff/Documents/Werk/STW/Data/geometries/Pig08/Pig08_ventricles.tri');
% [VERv ITRIv]=loadtri('/Users/peteroosterhoff/Documents/Werk/STW/Data/geometries/Pig09/Pig09_ventricles.tri');

% surftype=loadmat('/Users/peteroosterhoff/Documents/Werk/STW/Data/geometries/Pig07/Ingetest/Pig07/Pig07_ventricles.typ');
% surftype=loadmat('/Users/peteroosterhoff/Documents/Werk/STW/Data/geometries/Pig08/Pig08_ventricles.typ');
% surftype=loadmat('/Users/peteroosterhoff/Documents/Werk/STW/Data/geometries/Pig09/Pig09_ventricles.typ');

% [VERv ITRIv]=loadtri('/Users/peteroosterhoff/Documents/Werk/STW/Data/geometries/Pig09_refined/Pig09_refined_ventricles.tri');
% basepath='/Users/peteroosterhoff/Documents/Werk/STW/Data/results/ESC/Pig09_refined/ESC/single';
% surftype=loadmat('/Users/peteroosterhoff/Documents/Werk/STW/Data/geometries/Pig09_refined/Pig09_refined_ventricles.typ');

[VERv ITRIv]=loadtri('/Users/peteroosterhoff/Documents/Werk/STW/Data/geometries/Pig09/Pig09_ventricles.tri');
basepath='/Users/peteroosterhoff/Documents/Werk/STW/Data/results/Pig09/single';
surftype=loadmat('/Users/peteroosterhoff/Documents/Werk/STW/Data/geometries/Pig09/Pig09_ventricles.typ');


D=dir(fullfile(basepath,'*.mat'));
focusinv=[];
for i=1:length(D)
%     t=load(fullfile(basepath,D(i).name));
    if D(i).name(end-7:end)=='init.mat'
        t=load(fullfile(basepath,D(i).name),'measinit');
        Rinit{end+1}=t.measinit;
        fileinit{end+1}=D(i).name;
    else
        t=load(fullfile(basepath,D(i).name),'meas');
        Rinv{end+1}=t.meas;
        fileinv{end+1}=D(i).name;
    end
    
end
fileinv=fileinv';
fileinit=fileinit';

for i=1:length(Rinv)
    [~, focusinv(i,1)]=min(Rinv{i}.depfinal);
    focusinvxyz(i,:)=VERv(focusinv(i),:);
    rdinv(i,1)=Rinv{i}.rdfinal; 
    corinv(i,1)=Rinv{i}.corfinal;
end

for i=1:length(Rinit)
     focusinit(i,1)=Rinit{i}.foci;
     focusinitxyz(i,:)=VERv(focusinit(i),:);
     rdinit(i,1)=Rinit{i}.rd;
     corinit(i,1)=Rinit{i}.cor;
end




% 
% M{1,1}=D(1).name;
% M{2,1}=D(2).name;
% M{:,2}={focusinit};
% M{:,3}=focusinitxyz(:,1);
% M{:,4}=focusinitxyz(:,2);
% M{:,5}=focusinitxyz(:,3);
% M{:,6}=focusinv;
% M{:,7}=focusinvxyz(:,1);
% M{:,8}=focusinvxyz(:,2);
% M{:,9}=focusinvxyz(:,3);
% dlmwrite('test.xls',M);

% Surftype 1-7: Epi, LVEndo,RVEndo,Mitral Valve, Tricuspid Valve,RVOT,Aorta
surfname={'Epi','LVEndo','RVEndo','LVEndo','RVEndo','RVEndo','LVEndo'}; 


focusinitsurf=surfname(surftype(focusinit))';
focusinvsurf=surfname(surftype(focusinv))';
