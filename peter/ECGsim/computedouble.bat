@ echo off
double -g %1_thorax.tri 1 -g %1_lcav.tri 3 -g %1_rcav.tri 3 -g %1_ventricles.tri 1 -s %1_ventricles.tri -e @ -t %1_ventricles_edl.mat

rem double -g %1_thorax.tri 1 -g %1_lcav.tri 3 -g %1_rcav.tri 3 -g %1_llung.tri 3 -g %1_rlung.tri 3 -g %1_ventricles.tri 1 -s %1_ventricles.tri -e @ -t %1_ventricles_edl.mat
getSpc  