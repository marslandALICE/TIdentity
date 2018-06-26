# dataTree=/u/marsland/PHD/macros/marsland_EbyeRatios/schleching/data/trees/AnalysisResults_run1.root
# dataTree=/u/marsland/PHD/macros/marsland_EbyeRatios/schleching/data/trees/AnalysisResults_run2.root
# dataTree=/u/marsland/PHD/macros/marsland_EbyeRatios/schleching/data/trees/AnalysisResults_run2_new.root
# dataTree=/u/marsland/PHD/macros/marsland_EbyeRatios/schleching/data/trees/SubSample_cent0.00_ss1.root
# lineShapes=/u/marsland/PHD/macros/marsland_EbyeRatios/schleching/data/trees/LineShapes_ClonesArray.root
# lineShapes=/u/marsland/PHD/macros/marsland_EbyeRatios/schleching/data/trees/FitResults_piS0.7_kaS0.5_prS0.5_pikaprKauto.root


dataTree=/lustre/nyx/alice/users/marsland/pFluct/files/analysis/Data/PbPb/Real/RUN1/PHD/Systematics_cRows_80_16EtaBin_mombin20MeV/SubSamples/cent_0.0/SubSample_cent0.00_ss2.root
lineShapes=/lustre/nyx/alice/users/marsland/pFluct/files/analysis/Data/PbPb/Real/RUN1/PHD/Systematics_cRows_80_16EtaBin_mombin20MeV/ParamTrees/FitResults_piS0.7_kaS0.5_prS0.5_pikaprKauto.root

# ./testIden $dataTree $lineShapes $fSubsample $fSign $fCent $fpDown  $fpUp $fEtaDown $fEtaUp $fSystematic
./testIden   $dataTree $lineShapes      1        0       0     0.4     0.8    -0.5      0.5         0


# $fSystematic
# 0 -->  Reference
# 1 -->  CRows60
# 2 -->  CRows100
# 3 -->  Chi2TPC3
# 4 -->  Chi2TPC5
# 5 -->  DCAXYSmall
# 6 -->  DCAXYLarge
# 7 -->  VZSmall
# 8 -->  VZLarge
# 9 -->  EventVertexZSmall
# 10 --> EventVertexZLarge
# 11 --> ClusterRequirementITS
# 12 --> NewITSCut
