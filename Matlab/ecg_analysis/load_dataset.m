num_subjects = 1; 

for i = 1:num_subjects
    subject = sprintf('IKEM_Pat%03d', i);
    printProgress(i, num_subjects, subject);

    path = 'C:\Users\Admin\Documents\Projects\ecg_project\Scripts\data\Dataset\';
    DATA = readGeomPeacsModelDataset(path, subject);
end

q = initQtripy();
q.reset();
q.surface(DATA.VENTR.geom.VER, DATA.VENTR.geom.ITRI, "ventricles");
q.transparency(0.3);
