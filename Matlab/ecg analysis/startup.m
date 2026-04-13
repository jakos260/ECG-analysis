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
addpath (append(dir_name, 'ecg analysis/helpers'))
addpath (append(dir_name, 'ecg analysis/loaders'))
addpath (append(dir_name, 'ecg analysis/loaders/OptimizerData'))
addpath (append(dir_name, 'ecg analysis/qtriplot'))
addpath (append(dir_name, 'ecg analysis/optimizer'))
addpath (append(dir_name, 'ecg analysis/processing'))
addpath (append(dir_name, 'ecg analysis/tmpGenerator'))
addpath (append(dir_name, 'ecg analysis/tmpGenerator/ui_adjuster'))
addpath (append(dir_name, 'ecg analysis/ap'))