dataTree=/u/marsland/PHD/macros/marsland_EbyeRatios/schleching/data/trees/AnalysisResults.root
# dataTree=/u/marsland/PHD/macros/marsland_EbyeRatios/schleching/data/trees/SubSample_cent0.00_ss1.root
# lineShapes=/u/marsland/PHD/macros/marsland_EbyeRatios/schleching/data/trees/LineShapes_ClonesArray.root
lineShapes=/u/marsland/PHD/macros/marsland_EbyeRatios/schleching/data/trees/FitResults_piS0.7_kaS0.5_prS0.5_pikaprKauto.root

# ./testIden $dataTree $lineShapes $fSubsample $fSign $fCent $fpDown $fpUp $fEtaDown $fEtaUp $fSystematic
./testIden $dataTree $lineShapes 1 1     0    0.2 1.   -0.8  0.8     2


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
