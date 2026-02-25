function DATACASES = readPollyresults(fn)

[basedir,name,~] = fileparts(fn);
if isempty(basedir)
    basedir = pwd;
end

XML=xml2struct(fn);
if isfield(XML,'POLLY')
    POLLY = XML.POLLY{2};
else
    POLLY = XML.POLLYMAPPER{2};
end

if length(POLLY.POLLYDATA.ECG_cases.ECG_CASE) == 1
    POLLYDATA = POLLY.POLLYDATA;
    POLLYCASE = POLLYDATA.ECG_cases.ECG_CASE;
    DATACASES.patientid=name;
    if isfield(POLLYDATA,'modelname')
        modelName = POLLYDATA.modelname.Text;
    end
    DATACASES.DATA = readCaseData(basedir,POLLYCASE,modelName);
    
    
    DATACASES.name = POLLYCASE.ECGcasename.Text;
%     dd=dir( fullfile(basedir,POLLY.POLLYDATA.patientId.Text,'models') );
    dd=dir( fullfile(basedir,'models') );
    for i=1:length(dd)
        if isempty(strfind(dd(i).name,'.'))
            break
        end
    end

    dleadNames = fullfile(basedir,'ECGcases',POLLYCASE.ECGcasename.Text,'Leadsys');
    
    DATACASES.MODEL = readGeomPeacsModel(fullfile(basedir,'models',modelName),...
                                           modelName,...
                                           dleadNames);
else
    iCase=1;
    for k=1:length(POLLY.POLLYDATA.ECG_cases.ECG_CASE)
        POLLYCASE = POLLY.POLLYDATA.ECG_cases.ECG_CASE{k};
        DATA = readCaseData(fullfile(basedir,POLLY.POLLYDATA.patientId.Text),POLLYCASE);
        if ~isempty(DATA)
            DATACASES.DATA = DATA;

            DATACASES.patientid=name;

            DATACASES.name = POLLYCASE.ECGcasename.Text;
            
            dd=dir( fullfile(basedir,POLLY.POLLYDATA.patientId.Text,'models') );
            for i=1:length(dd)
                if isempty(strfind(dd(i).name,'.'))
                    break
                end
            end

            dleadNames = fullfile(basedir,POLLY.POLLYDATA.patientId.Text,'ECGcases',POLLYCASE.ECGcasename.Text,'Leadsys');

            DATACASES.MODEL = readGeomPeacsModel(fullfile(basedir,POLLY.POLLYDATA.patientId.Text,'models',POLLYCASE.modelname.Text),...
                            POLLYCASE.modelname.Text,...
                            dleadNames);
        end
        
    end
end




%%
function val= removeText(str)

if isfield( str,'Text')
    val = str2double(strsplit(str.Text,','));
else
    val = str2double(struct2array(str));
end


%%
function DATA = readCaseData(basedir,POLLYCASE,modelname)

DATA=[];
if isfield(POLLYCASE,'ECGfiles')
    for iEcg=1:length(POLLYCASE.ECGfiles.ECGFiledata)
        if length(POLLYCASE.ECGfiles.ECGFiledata) == 1
            ecg=POLLYCASE.ECGfiles.ECGFiledata;
        else
            ecg=POLLYCASE.ECGfiles.ECGFiledata{iEcg};
        end
        ECG = loadmat(fullfile(basedir,'ECGcases',POLLYCASE.ECGcasename.Text,'ECG',[ecg.ECGfilename.Text '.ecg']));
        if size(ECG,1) ==1 % was empty
            ECG = loadmat(fullfile(basedir,'ECGcases',POLLYCASE.ECGcasename.Text,'ECG',[ecg.ECGfilename.Text '.cipsecg']));
        end
        
        DATA{iEcg}.ecg.ECG = ECG;
        DATA{iEcg}.ecg.filename = ecg.ECGfilename.Text;
        
        leadsysnames= dir(fullfile(basedir,'ECGcases',POLLYCASE.ECGcasename.Text,'Leadsys','*.lead'));
%         for i=1:length(leadsysnames)
%             ii=strfind(leadsysnames(i).name,'_');            
%             if strcmp(modelname,leadsysnames(i).name(1:ii(end)-1))
%                 DATA{iEcg}.leadsysElec = loadmat(fullfile(basedir,'ECGcases',POLLYCASE.ECGcasename.Text,'Leadsys',leadsysnames(i).name));
%                 DATA{iEcg}.leadsysElec = DATA{iEcg}.leadsysElec(:,2:4);    
%             end
%         end
%         if ~isfield(DATA{iEcg},'leadsysElec')
%             DATA{iEcg}.leadsysElec = loadmat(fullfile(basedir,'ECGcases',POLLYCASE.ECGcasename.Text,'Leadsys',leadsysnames(end).name));
%             DATA{iEcg}.leadsysElec = DATA{iEcg}.leadsysElec(:,2:4);                
%             warning('did not find a matching leadsystem file with the model');
%         end
        
        if isfield(ecg,'selectedVentricularBeats')
            beats = ecg.selectedVentricularBeats.CyncRESULT;
            ecg =rmfield(ecg,[{'useECGSignal'} {'selectedVentricularBeats'} {'Attributes'}]);
            
            
            i=0;
            for k=1:length(beats)
                if length(beats)==1
                    selbeat = beats;
                else
                    selbeat = beats{k};                
                end
                if isfield(selbeat,'initdep') %&& isfield(beats{k},'finalrep')
                    i=i+1;
                    SPECS.onsetqrs      = removeText(selbeat.onsetQRS);
                    SPECS.endtwave      = removeText(selbeat.Twave_end);
                    SPECS.time_Jpoint   = removeText(selbeat.EndQRS);
                    SPECS.time_apexT    = removeText(selbeat.PeakTwave);
                    SPECS.time_apexU    = -1;
                    SPECS.depSlope      = 2;
                    SPECS.plateauslope  = 0.014;
                    SPECS.repslope      = 0.045;
                    if isfield(selbeat,'initialSlope')
                        SPECS.initialSlope  = removeText(selbeat.initialrepslope);
                    end
                    if isfield(selbeat,'plateauslope')
                        SPECS.plateauslope  = removeText(selbeat.platslope);
                    end
                    SPECS.repCorrection = 0;
                    SPECS.useCumsum     = 0;
                    SPECS.qrsduration   = SPECS.time_Jpoint - SPECS.onsetqrs;
                    SPECS.qrstduration  = SPECS.endtwave - SPECS.onsetqrs;
                    
                    DATA{iEcg}.beats{i}.SPECS           = SPECS;                    
                    DATA{iEcg}.beats{i}.initdep         = removeText(selbeat.initdep)';
                    if isfield(selbeat,'finaldep')
                        DATA{iEcg}.beats{i}.finaldep    = removeText(selbeat.finaldep)';
                    end
                    if isfield(selbeat,'initialcorrelation')
                        DATA{iEcg}.beats{i}.initialcorrelation    = removeText(selbeat.initialcorrelation)';
                    end
                    if isfield(selbeat,'finalrep')
                        DATA{iEcg}.beats{i}.rep         = removeText(selbeat.finalrep)';
                    end
                    if ( isfield(selbeat,'truthorigin') )
                        DATA{iEcg}.beats{i}.truthorigin = removeText(selbeat.truthorigin)';
                    end                    
                    DATA{iEcg}.beats{i}.ECG             = ECG(:,max(1,SPECS.onsetqrs):min(SPECS.endtwave,size(ECG,2)));
                end            
            end
        end
        
        if isfield(ecg,'selectedAtrialBeats')
            beats = ecg.selectedAtrialBeats.CyncRESULT;
           
            i=0;
            for k=1:length(beats)
                if isfield(beats{k},'finaldep') && isfield(beats{k},'finalrep')
                    i=i+1;
                    selbeat = beats{k};

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
        end
    end
else
    disp('no analysed data found')
end