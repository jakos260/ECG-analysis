%
% ==========================================================
%               ||  GET MAPS  ||
% ==========================================================

subjects = [1,2];

for i = 1:size(subjects,2)
    subject_num = subjects(i);
    subject = sprintf('Pat%03d_12ECG.imap', subject_num);
    fprintf(newline+"_________ processing %s ________"+newline, subject);

    path = "C:\Users\Admin\Documents\Projects\ecg_project\Scripts\data\raw\Mapper\";
    MapResultPath = fullfile(path, subject);
    
    nRV = 4; % 3 for LV pacing or 4 for RV pacing
    nLV = 3;

    MAP = readPollyMapresults(MapResultPath);
end
