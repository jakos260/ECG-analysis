if exist('local_paths.m', 'file')
    run('local_paths.m');
    disp(['Working directory is: ', dir_name]);
else
    warning('File local_paths.m not found! Define the variable "dir_name".');
    dir_name = ''; % fallback
end


addpath (append(dir_name, 'sim'))
addpath (append(dir_name, 'pnet_dir'))
addpath (append(dir_name, 'progs'))
addpath (append(dir_name, 'peter'))
addpath (append(dir_name, 'peter/inverse/'))
addpath (append(dir_name, 'avo/progs'))
addpath (append(dir_name, 'avo/electrog'))
addpath (append(dir_name, 'ecg_analysis/helpers'))
addpath (append(dir_name, 'ecg_analysis/loaders'))
addpath (append(dir_name, 'ecg_analysis/loaders/OptimizerData'))
addpath (append(dir_name, 'ecg_analysis/qtriplot'))
addpath (append(dir_name, 'ecg_analysis/optimizer'))
addpath (append(dir_name, 'ecg_analysis/processing'))
addpath (append(dir_name, 'ecg_analysis/tmpGenerator'))
addpath (append(dir_name, 'ecg_analysis/tmpGenerator/ui_adjuster'))
addpath (append(dir_name, 'ecg_analysis/ap'))
addpath (append(dir_name, 'ecg_analysis/inverseMy'))
addpath (genpath (append (dir_name, 'ecg_analysis/inverseArno')))