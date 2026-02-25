function DATACASES = readCIPSresults(fn)


[basedir,name,~] = fileparts(fn);

XML=xml2struct(fn);
if length(XML.CIPS{2}.CIPSDATA.ECG_cases.CIPS_ECG_CASE) == 1
    CIPSCASE = XML.CIPS{2}.CIPSDATA.ECG_cases.CIPS_ECG_CASE;
    DATACASES.patientid=name;
    DATACASES.DATA = readCaseData(basedir,CIPSCASE);
    
    
    DATACASES.name = CIPSCASE.ECGcasename.Text;
    dd=dir( fullfile(basedir,XML.CIPS{2}.CIPSDATA.patientId.Text,'models') );
    for i=1:length(dd)
        if isempty(strfind(dd(i).name,'.'))
            break
        end
    end

    dleadNames = fullfile(basedir,'ECGcases',CIPSCASE.ECGcasename.Text,'Leadsys');
    if isfield(CIPSCASE,'modelname')
        modelN = CIPSCASE.modelname.Text;
    else
        modelN = XML.CIPS{2}.CIPSDATA.modelname.Text;
    end
    DATACASES.MODEL = readGeomPeacsModel(fullfile(basedir,'models',modelN),...
                                           modelN,...
                                           dleadNames);
else
    iCase=1;
    for k=1:length(XML.CIPS{2}.CIPSDATA.ECG_cases.CIPS_ECG_CASE)
        CIPSCASE = XML.CIPS{2}.CIPSDATA.ECG_cases.CIPS_ECG_CASE{k};
        DATA = readCaseData(fullfile(basedir,XML.CIPS{2}.CIPSDATA.patientId.Text),CIPSCASE);
        if ~isempty(DATA)
            DATACASES.DATA = DATA;

            DATACASES.patientid=name;

            DATACASES.name = CIPSCASE.ECGcasename.Text;
            
            dd=dir( fullfile(basedir,XML.CIPS{2}.CIPSDATA.patientId.Text,'models') );
            for i=1:length(dd)
                if isempty(strfind(dd(i).name,'.'))
                    break
                end
            end

            dleadNames = fullfile(basedir,XML.CIPS{2}.CIPSDATA.patientId.Text,'ECGcases',CIPSCASE.ECGcasename.Text,'Leadsys');

            DATACASES.MODEL = readGeomPeacsModel(fullfile(basedir,XML.CIPS{2}.CIPSDATA.patientId.Text,'models',CIPSCASE.modelname.Text),...
                            CIPSCASE.modelname.Text,...
                            dleadNames);
        end
        
    end
end




%%
function val= removeText(str)

val = str2num(struct2array(str));

%%
function DATA = readCaseData(basedir,CIPSCASE)

DATA=[];
if isfield(CIPSCASE,'ECGfiles')
    for iEcg=1:length(CIPSCASE.ECGfiles.ECGFiledata)
        if length(CIPSCASE.ECGfiles.ECGFiledata) == 1
            ecg=CIPSCASE.ECGfiles.ECGFiledata;
        else
            ecg=CIPSCASE.ECGfiles.ECGFiledata{iEcg};
        end
        ECG = loadmat(fullfile(basedir,'ECGcases',CIPSCASE.ECGcasename.Text,'ECG',[ecg.ECGfilename.Text '.cipsecg']));
        DATA{iEcg}.ecg.ECG = ECG;
        DATA{iEcg}.ecg.filename = ecg.ECGfilename.Text;
        DATA{iEcg}.ecg.onsetPwaves = removeText(ecg.onsetPwaves);
        DATA{iEcg}.ecg.endPwaves = removeText(ecg.endPwaves);
        DATA{iEcg}.ecg.onsetQRSs = removeText(ecg.onsetQRSs);
        DATA{iEcg}.ecg.endQRSs = removeText(ecg.endQRSs);
        DATA{iEcg}.ecg.endTwaves = removeText(ecg.endTwaves);
        DATA{iEcg}.ecg.peakTwaves = removeText(ecg.peakTwaves);
        
        leadsysnames= dir(fullfile(basedir,'ECGcases',CIPSCASE.ECGcasename.Text,'Leadsys','*.lead'));
        if length(leadsysnames)==1
            DATA{iEcg}.leadsysElec = loadmat(fullfile(basedir,'ECGcases',CIPSCASE.ECGcasename.Text,'Leadsys',leadsysnames.name));
            DATA{iEcg}.leadsysElec = DATA{iEcg}.leadsysElec(:,2:4);
        else            
            DATA{iEcg}.leadsysElec = loadmat(fullfile(basedir,'ECGcases',CIPSCASE.ECGcasename.Text,'Leadsys',leadsysnames(1).name));
            DATA{iEcg}.leadsysElec = DATA{iEcg}.leadsysElec(:,2:4);

            warning('more than one lead file')
        end
        if isfield(ecg,'selectedVentricularBeats')
            beats = ecg.selectedVentricularBeats.CIPS_OPTIMIZERESULT_CASE;
            ecg =rmfield(ecg,[{'useECGSignal'} {'selectedVentricularBeats'} {'Attributes'}]);
            
            
            i=0;
            if length(beats)==1
                if isfield(beats,'initdep') 
                    i=i+1;
                    selbeat = beats.SelectedECGbeat;
                    
                    SPECS.onsetP        = removeText(selbeat.Pwave_onset);
                    SPECS.onsetqrs      = removeText(selbeat.onsetQRS);
                    SPECS.endtwave      = removeText(selbeat.Twave_end);
                    SPECS.time_Jpoint   = removeText(selbeat.EndQRS);
                    SPECS.time_apexT    = removeText(selbeat.PeakTwave);
                    SPECS.time_apexU    = -1;
                    SPECS.depSlope      = 2;
                    SPECS.initialSlope  = removeText(selbeat.initialrepslope);
                    
                    SPECS.plateauslope  = removeText(selbeat.platslope);
                    if SPECS.plateauslope == 0
                        SPECS.plateauslope = 0.014;
                    end
                    SPECS.repslope      = removeText(selbeat.repslope);
                    if SPECS.repslope == 0
                        SPECS.repslope = 0.045;
                    end
                    SPECS.repCorrection = 0;
                    SPECS.useCumsum     = 0;
                    SPECS.qrsduration   = SPECS.time_Jpoint - SPECS.onsetqrs;
                    SPECS.qrstduration  = SPECS.endtwave - SPECS.onsetqrs;
                    
                    DATA{iEcg}.beats{i}.initdep = removeText(beats.initdep)';
                    if isfield(beats,'finaldep')
                        DATA{iEcg}.beats{i}.dep = removeText(beats.finaldep)';
                    end
                    if isfield(beats,'finalrep')
                        DATA{iEcg}.beats{i}.rep = removeText(beats.finalrep)';
                    end
                    if ( isfield(beats,'truthorigin') )
                        DATA{iEcg}.beats{i}.truthorigin = removeText(beats.truthorigin)';
                    end
                    DATA{iEcg}.beats{i}.SPECS        = SPECS;
                    DATA{iEcg}.beats{i}.ECG         = ECG(:,max(1,SPECS.onsetqrs):min(SPECS.endtwave,size(ECG,2)));
                end
            else
                for k=1:length(beats)
                    if isfield(beats{k},'initdep') %&& isfield(beats{k},'finalrep')
                        i=i+1;
                        selbeat = beats{k}.SelectedECGbeat;
                        
                        SPECS.onsetP        = removeText(selbeat.Pwave_onset);
                        SPECS.onsetqrs      = removeText(selbeat.onsetQRS);
                        SPECS.endtwave      = removeText(selbeat.Twave_end);
                        SPECS.time_Jpoint   = removeText(selbeat.EndQRS);
                        SPECS.time_apexT    = removeText(selbeat.PeakTwave);
                        SPECS.time_apexU    = -1;
                        SPECS.depSlope      = 2;
                        SPECS.initialSlope  = removeText(selbeat.initialrepslope);
                        
                        SPECS.plateauslope  = removeText(selbeat.platslope);
                        if SPECS.plateauslope == 0
                            SPECS.plateauslope = 0.014;
                        end
                        SPECS.repslope      = removeText(selbeat.repslope);
                        if SPECS.repslope == 0
                            SPECS.repslope = 0.045;
                        end
                        SPECS.repCorrection = 0;
                        SPECS.useCumsum     = 0;
                        SPECS.qrsduration   = SPECS.time_Jpoint - SPECS.onsetqrs;
                        SPECS.qrstduration  = SPECS.endtwave - SPECS.onsetqrs;
                        DATA{iEcg}.beats{i}.dep = removeText(beats{k}.initdep)';
                        if isfield(beats{k},'meanTSI')
                            DATA{iEcg}.beats{i}.meanTSI = removeText(beats{k}.meanTSI)';
                        end
                        if isfield(beats{k},'finaldep')
                            DATA{iEcg}.beats{i}.findep = removeText(beats{k}.finaldep)';
                        end
                        if isfield(beats{k},'finalrep')
                            DATA{iEcg}.beats{i}.finrep = removeText(beats{k}.finalrep)';
                        end
                        if ( isfield(beats{k},'truthorigin') )
                            DATA{iEcg}.beats{i}.truthorigin = removeText(beats{k}.truthorigin)';
                        end
                        DATA{iEcg}.beats{i}.SPECS        = SPECS;
                        DATA{iEcg}.beats{i}.ECG         = ECG(:,max(1,SPECS.onsetqrs):min(SPECS.endtwave,size(ECG,2)));
                    end
                end
            end
        end
        
        if isfield(ecg,'selectedAtrialBeats')
            beats = ecg.selectedAtrialBeats.CIPS_OPTIMIZERESULT_CASE;
            %         ecg =rmfield(ecg,[{'useECGSignal'} {'selectedAtrialBeats'} {'Attributes'}]);
            %             DATA{iEcg}.ecg.filename = ecg.ECGfilename.Text;
            %             DATA{iEcg}.ecg.onsetPwaves = removeText(ecg.onsetPwaves);
            %             DATA{iEcg}.ecg.endPwaves = removeText(ecg.endPwaves);
            %             DATA{iEcg}.ecg.onsetQRSs = removeText(ecg.onsetQRSs);
            %             DATA{iEcg}.ecg.endQRSs = removeText(ecg.endQRSs);
            %             DATA{iEcg}.ecg.endTwaves = removeText(ecg.endTwaves);
            %             DATA{iEcg}.ecg.peakTwaves = removeText(ecg.peakTwaves);
            
            
            if length(beats)>1
                i=0;
                for k=1:length(beats)
                    if isfield(beats{k},'finaldep') && isfield(beats{k},'finalrep')
                        i=i+1;
                        selbeat = beats{k}.SelectedECGbeat;
                        
                        SPECS.onsetP        = removeText(selbeat.Pwave_onset);
                        SPECS.onsetqrs      = removeText(selbeat.onsetQRS);
                        SPECS.endtwave      = removeText(selbeat.Twave_end);
                        SPECS.time_Jpoint   = removeText(selbeat.EndQRS);
                        SPECS.time_apexT    = removeText(selbeat.PeakTwave);
                        SPECS.time_apexU    = -1;
                        SPECS.depSlope      = 2;
                        SPECS.initialSlope  = removeText(selbeat.initialrepslope);
                        
                        SPECS.plateauslope  = removeText(selbeat.platslope);
                        if SPECS.plateauslope == 0
                            SPECS.plateauslope = 0.014;
                        end
                        SPECS.repslope      = removeText(selbeat.repslope);
                        if SPECS.repslope == 0
                            SPECS.repslope = 0.045;
                        end
                        SPECS.repCorrection = 0;
                        SPECS.useCumsum     = 0;
                        SPECS.qrsduration   = SPECS.time_Jpoint - SPECS.onsetqrs;
                        SPECS.qrstduration  = SPECS.endtwave - SPECS.onsetqrs;
                        
                        DATA{iEcg}.beats{i}.dep = removeText(beats{k}.finaldep)';
                        DATA{iEcg}.beats{i}.rep = removeText(beats{k}.finalrep)';
                        DATA{iEcg}.beats{i}.SPECS        = SPECS;
                        DATA{iEcg}.beats{i}.ECG         = ECG(:,max(1,SPECS.onsetP):min(SPECS.endtwave,size(ECG,2)));
                    end
                end
            else
                if isfield(beats,'finaldep') && isfield(beats,'finalrep')
                    selbeat = beats.SelectedECGbeat;
                    
                    SPECS.onsetP        = removeText(selbeat.Pwave_onset);
                    SPECS.onsetqrs      = removeText(selbeat.onsetQRS);
                    SPECS.endtwave      = removeText(selbeat.Twave_end);
                    SPECS.time_Jpoint   = removeText(selbeat.EndQRS);
                    SPECS.time_apexT    = removeText(selbeat.PeakTwave);
                    SPECS.time_apexU    = -1;
                    SPECS.depSlope      = 2;
                    SPECS.initialSlope  = removeText(selbeat.initialrepslope);
                    
                    SPECS.plateauslope  = removeText(selbeat.platslope);
                    if SPECS.plateauslope == 0
                        SPECS.plateauslope = 0.014;
                    end
                    SPECS.repslope      = removeText(selbeat.repslope);
                    if SPECS.repslope == 0
                        SPECS.repslope = 0.045;
                    end
                    SPECS.repCorrection = 0;
                    SPECS.useCumsum     = 0;
                    SPECS.qrsduration   = SPECS.time_Jpoint - SPECS.onsetqrs;
                    SPECS.qrstduration  = SPECS.endtwave - SPECS.onsetqrs;
                    
                    DATA{iEcg}.beats.dep = removeText(beats.finaldep)';
                    DATA{iEcg}.beats.rep = removeText(beats.finalrep)';
                    DATA{iEcg}.beats.SPECS        = SPECS;
                    DATA{iEcg}.beats{i}.ECG         = ECG(:,max(1,SPECS.onsetP):min(SPECS.endtwave,size(ECG,2)));
                end
            end
        end
    end
else
    disp('no analysed data found')
end