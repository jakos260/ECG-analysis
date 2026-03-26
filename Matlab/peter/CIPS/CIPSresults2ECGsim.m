function CIPSresults2ECGsim(fn, ecgfile, dirname)

[DATA,MODEL] = readCIPSresults(fn);
foundA=0;

for i=1:length(DATA)
    casename = DATA{i}.ecg.filename;
    ECG = loadmat(fullfile(ecgfile,[casename '.cipsecg']));
    for ib=1:length(DATA{i}.beats)
        if ~exist(fullfile(dirname,casename))
            mkdir(fullfile(dirname,casename));
        end
        saveasci(fullfile(dirname,casename,['beat' num2str(ib) '.dep']),DATA{i}.beats{ib}.dep);
        saveasci(fullfile(dirname,casename,['beat' num2str(ib) '.rep']),DATA{i}.beats{ib}.rep);
        if DATA{i}.beats{ib}.SPECS.onsetqrs > 0
            ecg = ECG(:,max(1,DATA{i}.beats{ib}.SPECS.onsetqrs-200):DATA{i}.beats{ib}.SPECS.endtwave);
            savemat(fullfile(dirname,casename,['beat' num2str(ib) '.ecg']),ecg(4:12,:));
        elseif DATA{i}.beats{ib}.SPECS.onsetP > 0
            ecg = ECG(:,DATA{i}.beats{ib}.SPECS.onsetP):min(size(ECG,2),DATA{i}.beats{ib}.SPECS.onsetP+800);
            savemat(fullfile(dirname,casename,['beat' num2str(ib) '.ecg']),ecg(4:12,:));
            foundA=1;
        end           
        SPECS = [DATA{i}.beats{ib}.SPECS.plateauslope DATA{i}.beats{ib}.SPECS.repslope]';
        saveasci(fullfile(dirname,casename,['beat' num2str(ib) '.specs']),SPECS);        
    end    
end

if ~foundA
    adep = MODEL.ATRIA.DIST3D(:,1);
    adep=adep*110/max(adep);
    
    saveasci(fullfile(dirname,casename,['beat' num2str(ib) '.adep']),adep);
    saveasci(fullfile(dirname,casename,['beat' num2str(ib) '.arep']),200 - adep);
end