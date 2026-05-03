clearvars
% copy init code from create...

INITIALVELOCITY=[0,0.4,1.0];
INITIALVELOCITY=[0.4];

maxvelocity=[2.5];%1.6;
velo = 1.0;

ANIS = [0.01,1,2.5];
ANIS=[0.01 2.0 ]

INITANIS=[0.01,1,2.5];
INITANIS=[1.0];

MUDEP=[5e-5 5e-6 1e-6 5e-7]; % ignored when set to inf

% basepath='/Users/peteroosterhoff/Documents/Werk/STW/Data/results/Pig09/single20131113';
% basepath='/Users/peteroosterhoff/Documents/Werk/STW/Data/results/Pig09/single20131127Initmode1+dep0_-16plus';
% basepath='/Users/peteroosterhoff/Documents/Werk/STW/Data/results/Pig09/single20131127_Inverse+dep0-16plus';
% basepath='/Users/peteroosterhoff/Documents/Werk/STW/Data/results/Pig09/cluster16/Initmode0';
% basepath='/Users/peteroosterhoff/Documents/Werk/STW/Data/results/Pig09/cluster1/VeloScan';
% basepath='/Users/peteroosterhoff/Documents/Werk/STW/Data/results/Pig09/cluster1/NoVeloScan';
% basepath='/Users/peteroosterhoff/Documents/Werk/STW/Data/results/Pig07/cluster1/FullAVGScan';
% basepath='/Users/peteroosterhoff/Documents/Werk/STW/Data/results/Pig10/cluster1/FullAVGScan';
% basepath='/Users/peteroosterhoff/Documents/Werk/STW/Data/results/Pig07/cluster1/muScanAVG';
% basepath='/Users/peteroosterhoff/Documents/Werk/STW/Data/results/Pig09/cluster1/Paper';


% subj='Pig08'
for subj={'Pig07','Pig09','Pig08','Pig10'}
basepath=['/Users/peteroosterhoff/Documents/Werk/STW/Data/results_b/' subj{1} '/cluster1/minrd015_stuck_mode1_mudepscan'];
% basepath=['/Users/peteroosterhoff/Documents/Werk/STW/Data/results/' subj{1} '/cluster1/mudep5e-5_minrd0.15'];
DATE='';
% DATE='20140519';


initmode=0;
wrd=0;
SELSTRING={'','Endo','Epi'};%'Endo'; % If not empty, only select summary files containing this string (Endo, Epi, LV, RV)

header=[];

for selstringcell=SELSTRING
    Savg=[];
    selstring=cell2mat(selstringcell);
    for initialvelocity=INITIALVELOCITY;
        for initanis=INITANIS;
            for anis=ANIS
                for mudep=MUDEP
                    
                    if isempty(DATE)
                        if isinf(mudep)
                            dirpath= fullfile(basepath,sprintf('%s_im%d_wrd%d_iV%0.1f_iAnis%.01fAnis%0.2f','*',initmode,wrd,initialvelocity,initanis,anis)); %NOTE: DECIMAL PART OF ANIS LOST!!!
                        else
                            dirpath= fullfile(basepath,sprintf('%s_im%d_wrd%d_iV%0.1f_iAnis%.01fAnis%0.2fmudep%0.2e','*',initmode,wrd,initialvelocity,initanis,anis,mudep)); %NOTE: DECIMAL PART OF ANIS LOST!!!
                        end
                        D=dir(dirpath);
                        DATE=D(1).name(1:8);
                    end

                    
                    
                    if isinf(mudep)
                        summpath= fullfile(basepath,sprintf('%s_im%d_wrd%d_iV%0.1f_iAnis%.01fAnis%0.2f',DATE,initmode,wrd,initialvelocity,initanis,anis), ...
                            sprintf('%s_im%d_wrd%d_iV%0.1f_iAnis%.01fAnis%0.2f%s_summary.mad',DATE,initmode,wrd,initialvelocity,initanis,anis,selstring)); %NOTE: DECIMAL PART OF ANIS LOST!!!
                    else
                        summpath= fullfile(basepath,sprintf('%s_im%d_wrd%d_iV%0.1f_iAnis%.01fAnis%0.2fmudep%0.2e',DATE,initmode,wrd,initialvelocity,initanis,anis,mudep), ...
                            sprintf('%s_im%d_wrd%d_iV%0.1f_iAnis%.01fAnis%0.2fmudep%0.2e%s_summary.mad',DATE,initmode,wrd,initialvelocity,initanis,anis,mudep,selstring)); %NOTE: DECIMAL PART OF ANIS LOST!!!
                    end
                    %             summpath='/Users/peteroosterhoff/Documents/Werk/STW/Data/results/Pig09/single201231113/20131111_im1_wrd1_iV0.4_iAnis1.0Anis0.4/20131111_im1_wrd1_iV0.4_iAnis1.0Anis0.4_Summary.mad';
                    s=load(summpath,'-mat');
                    %             if isempty(header)
                    %                 header=fieldnames;
                    %             end
                    savg=s.avg(end-1);
                    n.initialvelocity=initialvelocity;
                    n.initanis=initanis;
                    n.anis=anis;
                    n.mudep=mudep;
                    %             pairs = [fieldnames(n), struct2cell(n); fieldnames(savg), struct2cell(savg)].';
                    sel = catstruct(n,savg);
                    if isempty(Savg)
                        Savg=sel;
                    else
                        Savg(end+1)=sel;
                    end
                end
            end
        end
    end
    fields=fieldnames(Savg);
    for i=1:length(fields)
        if isempty(getfield(Savg,fields{i}));
            Savg=rmfield(Savg,fields(i));
        end
    end
    WriteStructsToText(fullfile(basepath,['OverallSummary' selstring subj{1} '.txt']),Savg);
    save(fullfile(basepath,['OverallSummary' selstring subj{1} '.mad']),'Savg');
end
end