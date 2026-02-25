clear all

q = initQtripy();

patient = 'normal_young_male'; offset = 169;
[ventri_ver, ventri_tri] = loadtri_ecgsim(append('data/ExportData/', patient, '/model/ventricle.tri'));
rep     = loadmat(append('data/ExportData/', patient, '/ventricular_beats/beat1/user.rep'));

Vpy = ventri_ver;
Tpy = ventri_tri;
% py.numpy.array(ventri_ver);
% py.numpy.array(ventri_tri);

q.cmd('reset')
q.surface(Vpy, Tpy, "ventricles");
q.transparency(0.3);
q.values(rep)
q.gradient_bins(20)