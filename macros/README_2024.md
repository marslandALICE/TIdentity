# Datasets

| Dataset | System | List                        |
| ------- | ------ | --------------------------- |
| LHC18q  | Pb–Pb  | runs-2018-LHC18q-pass3.list |
| LHC18r  | Pb–Pb  | runs-2018-LHC18r-pass3.list |

# Productions

| Production  | System  | Period  | JIRA         | Details                                | Lists                                   |
| ----------- | ------- | ------- | ------------ | -------------------------------------- | --------------------------------------- |
| LHC20e3abc  | Pb–Pb   | LHC18qr | ALIROOT-8462 | General purpose HIJING min bias        | runsMC-2020-LHC20e3[a,b,c]-pass3.list   |
| LHC22b5     | Pb–Pb   | LHC18qr | ALIROOT-8783 | General purpose HIJING 50-90%          | runsMC-2022-LHC22b5-pass3.list          |
| LHC20k6abcd | Pb–Pb   | LHC18qr | ALIROOT-8580 | Cascade injected HIJING                | runsMC-2020-LHC20k6[a,b,c,d]-pass3.list |
| LHC22b7abcd | Pb–Pb   | LHC18qr | ALIROOT-8785 | V0 and cascade injected HIJING, GEANT4 | runsMC-2022-LHC22b7[a,b,c,d]-pass3.list |
| LHC20g4     | pp 5TeV | LHC18qr | ALIROOT-8505 | jet-jet, no ESDs                       |                                         |

# How to run runGrid.c
cd $RUN_ON_GRID_DIR/Ebye/test_fastGen_netParticles/
cp $RUN_ON_GRID_DIR/Ebye/code/*.* $RUN_ON_GRID_DIR/Ebye/code/READ*  .; rm *.so *.d

# runGrid.C macro parameters
void runGrid(Bool_t fRunLocalFiles = kTRUE,
             TString mode="test",
             Int_t localOrGrid=0,
             Int_t setType=3,
             TString list = "",
             Int_t isMC=0,
             Int_t lhcYear=2015,
             TString periodName="15o",
             Int_t passIndex=2,
             Int_t valgrindOption = 0
             )
# Run real data
aliroot -b -q 'runGrid.C(0,"test",0,  0, "$RUN_ON_GRID_DIR/Ebye/lists/runsIlya1run-2018-LHC18r-pass3.list",0,2018,"18r",3,0)'


# run full MC
aliroot -b -q 'runGrid.C(0,"test",0,  50, "$RUN_ON_GRID_DIR/Ebye/lists/runsTest-2020-LHC20e3a-pass3.list",1,2018,"18q",3,0)'


# run fastGen
aliroot -b -q 'runGrid.C(0,"test",0,  200, "$RUN_ON_GRID_DIR/Ebye/lists/runsGen1run-2022-LHC22d1c2-pass3.list",2,2022,"22d1c2",2,0)'

# copy data
alien_cp -T 6 -parent 99 -glob AnalysisResults.root /Alice/cern.ch/user/m/marsland/PWGPP695_MC_remapping/LHC20e3a_pass3_20240625_228/2020 file:




