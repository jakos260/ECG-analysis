@ echo off

 ama -g Bucket01_thorax.tri 6 -g Bucket01_ventricles.tri 1 -s Bucket01_ventricles.tri  -e @ -t Bucket01_ventricles_edl.mat

@rem ama -g Bucket01_thorax.tri 3 -g Bucket01_ventricles.tri 1 -s Bucket01_ventricles.tri  -e @ -t Bucket01_ventricles_edl.mat

 @rem ama -g %1_thoraxpigs_adam.tri 1 -g %1_lcav.tri 3 -g %1_rcav.tri 3 -g %1_rlung.tri 0.2 -g %1_llung.tri 0.2 -g %1_ventricles.tri 1 -s %1_ventricles.tri -e @ -t %1_ventricles_edl.mat