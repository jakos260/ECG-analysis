global DATA;
num_subjects = 1; 

subject = sprintf('IKEM_Pat%03d', num_subjects);
path = 'C:\Users\Admin\Documents\Projects\ecg_project\Scripts\data\Dataset\';
bsm_name = 'x3_2019_11_15_12_36_17';
DATA = readGeomPeacsModelDataset(path, subject);
BSM = DATA.VENTR.SIGNALS.BSM.(bsm_name);
LAY  = loadmat('C:\Users\Admin\Documents\Projects\ecg_project\Scripts\Matlab\ecg_analysis\inverseArno\BEM\inverse\mla\prague99.mla');


[BSMOUT, badSigs, TIMES] = selectBSMBeats(BSM, LAY);

