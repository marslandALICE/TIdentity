#./testIdenTest PbPb/Real/PbPb_cRows_80_16EtaBin_mombin20MeV_LooseCuts_PxPyPz/SubSamples/cent_10.0/SubSample_cent10.00_ss1.root PbPb/Real/PbPb_cRows_80_16EtaBin_mombin20MeV_LooseCuts_PxPyPz/ParamTrees/FitResults_piS1_kaS1_prS1_piK1_pikaprKauto.root 10 1 6 -0.8 0.8 -1

# ./testIden PbPb/Real/RUN1/PHD/Systematics_cRows_80_16EtaBin_mombin20MeV/SubSamples/cent_0.0/SubSample_cent0.00_ss1.root PbPb/Real/RUN1/PHD/Systematics_cRows_80_16EtaBin_mombin20MeV/ParamTrees/FitResults_piS0.7_kaS0.5_prS0.5_pikaprKauto.root 0 1 6 0 0 1 0  0 0

#  ./testIden PbPb/Real/RUN1/PHD/Systematics_cRows_80_16EtaBin_mombin20MeV/SubSamples/cent_0.0/SubSample_cent0.00_ss1.list PbPb/Real/RUN1/PHD/Systematics_cRows_80_16EtaBin_mombin20MeV/ParamTrees/FitResults_piS0.7_kaS0.5_prS0.5_pikaprKauto.root 0 1 6 0 0 1 0  0 0

# ./testIden PbPb/Real/RUN1/PHD/Systematics_cRows_80_16EtaBin_mombin20MeV/SubSamples/cent_0.0/SubSample_cent0.00_ss1.root PbPb/Real/RUN1/PHD/Systematics_cRows_80_16EtaBin_mombin20MeV/ParamTrees/LineShapes_6.root 0 1 6 0 0 1 0  0 0

# ./testIden trees/SubSample_cent0.00_ss1.root trees/LineShapes_ClonesArray.root 0 1 6 1     0 1 0  0 0
 ./testIden trees/SubSample_cent0.00_ss1.root trees/FitResults_piS0.7_kaS0.5_prS0.5_pikaprKauto.root 0 1 6 1     0 1 0  0 0


#./testIdenTest PbPb/Real/PbPb_cRows_80_16EtaBin_mombin20MeV_TightCuts_PxPyPz/SubSamples/cent_40.0/SubSample_cent40.00_ss0.root PbPb/Real/PbPb_cRows_80_16EtaBin_mombin20MeV_TightCuts/ParamTrees/ParamTree_piS0.7_kaS0.5_prS0.5_Kurtosis2_fullStat.root 40 1 6 -0.8 0.8 -1









#        sprintf(fileName,"%s",argv[1]);
#        sprintf(fitFileName,"%s",argv[2]);
#        centInput    = atoi(argv[3]);
#        subsample    = atoi(argv[4]);
#        analyseIter  = atoi(argv[5]);
#        charge       = atoi(argv[6]);
#        systematics  = atoi(argv[7]);
#        isIntegrated = atoi(argv[8]);
#        rapInterval  = atoi(argv[9]);
#        momInterval  = atoi(argv[10]);
#        isSim        = atoi(argv[11]);

