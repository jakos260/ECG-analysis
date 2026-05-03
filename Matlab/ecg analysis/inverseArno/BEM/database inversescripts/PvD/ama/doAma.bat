@ echo off


ama -g %1_thorax65_Nijmegen.tri 1 -g %1_lcav.tri 3 -g %1_rcav.tri 3 -g %1_rlung.tri 0.2 -g %1_llung.tri 0.2 -g %1_ventricles.tri 1 -g %1_atria.tri     1 -s %1_atria.tri     -e @ -t %1_atria_edl.mat
ama -g %1_thorax65_Nijmegen.tri 1 -g %1_lcav.tri 3 -g %1_rcav.tri 3 -g %1_rlung.tri 0.2 -g %1_llung.tri 0.2 -g %1_atria.tri     1 -g %1_ventricles.tri 1 -s %1_ventricles.tri -e @ -t %1_ventricles_edl.mat


