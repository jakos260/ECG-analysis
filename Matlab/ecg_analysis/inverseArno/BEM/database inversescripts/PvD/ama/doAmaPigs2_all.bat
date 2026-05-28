@ echo off

ama -g %1_thoraxpigs_adam.tri 1 -g %1_lcav.tri 3 -g %1_rcav.tri 3 -g %1_rlung.tri 0.2 -g %1_llung.tri 0.2 -g %1_liver.tri 0.25 -g %1_ventricles.tri 1 -g %1_atria.tri     1 -s %1_atria.tri     -e @ -t %1_atria_edl_l.mat
ama -g %1_thoraxpigs_adam.tri 1 -g %1_lcav.tri 3 -g %1_rcav.tri 3 -g %1_rlung.tri 0.2 -g %1_llung.tri 0.2 -g %1_liver.tri 0.25 -g %1_atria.tri     1 -g %1_ventricles.tri 1 -s %1_ventricles.tri -e @ -t %1_ventricles_edl_l.mat

ama -g %1_thoraxpigs_adam.tri 1 -g %1_lcav.tri 3 -g %1_rcav.tri 3 -g %1_rlung.tri 0.2 -g %1_llung.tri 0.2 -g %1_ventricles.tri 1 -g %1_atria.tri     1 -s %1_atria.tri     -e @ -t %1_atria_edl_n.mat
ama -g %1_thoraxpigs_adam.tri 1 -g %1_lcav.tri 3 -g %1_rcav.tri 3 -g %1_rlung.tri 0.2 -g %1_llung.tri 0.2 -g %1_atria.tri     1 -g %1_ventricles.tri 1 -s %1_ventricles.tri -e @ -t %1_ventricles_edl_n.mat

ama -g %1_thoraxpigs_adam.tri 1 -g %1_rlung.tri 0.2 -g %1_llung.tri 0.2 -g %1_liver.tri 0.2 -g %1_ventricles.tri 1 -g %1_atria.tri     1 -s %1_atria.tri     -e @ -t %1_atria_edl_nb.mat
ama -g %1_thoraxpigs_adam.tri 1 -g %1_rlung.tri 0.2 -g %1_llung.tri 0.2 -g %1_liver.tri 0.2 -g %1_atria.tri     1 -g %1_ventricles.tri 1 -s %1_ventricles.tri -e @ -t %1_ventricles_edl_nb.mat

ama -g %1_thoraxpigs_adam.tri 1 -g %1_lcav.tri 3 -g %1_rcav.tri 3 -g %1_ventricles.tri 1 -g %1_atria.tri     1 -s %1_atria.tri     -e @ -t %1_atria_edl_nl.mat
ama -g %1_thoraxpigs_adam.tri 1 -g %1_lcav.tri 3 -g %1_rcav.tri 3 -g %1_atria.tri     1 -g %1_ventricles.tri 1 -s %1_ventricles.tri -e @ -t %1_ventricles_edl_nl.mat

ama -g %1_thoraxpigs_adam.tri 1 -g %1_ventricles.tri 1 -g %1_atria.tri     1 -s %1_atria.tri     -e @ -t %1_atria_edl_h.mat
ama -g %1_thoraxpigs_adam.tri 1 -g %1_atria.tri     1 -g %1_ventricles.tri 1 -s %1_ventricles.tri -e @ -t %1_ventricles_edl_h.mat
