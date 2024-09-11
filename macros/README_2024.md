# Datasets

| Dataset | System    | List                        |
| ------- | --------  | --------------------------- |
| LHC18q  | Pb–Pb     | runs-2018-LHC18q-pass3.list |
| LHC18r  | Pb–Pb     | runs-2018-LHC18r-pass3.list |
| LHC18b  | pp 13 TeV | runs-2018-LHC18b-pass2.list |

# Productions

| Production  | System   | Period  | JIRA         | Details                                  | Lists                                   | Events                |
| ----------- | -------- | ------- | ------------ | ---------------------------------------- | --------------------------------------- | --------------------- |
| LHC20e3abc  | Pb–Pb    | LHC18qr | ALIROOT-8462 | General purpose HIJING min bias          | runsMC-2020-LHC20e3[a,b,c]-pass3.list   | 3.6M, 0.3M, 3M        |
| LHC22b5     | Pb–Pb    | LHC18qr | ALIROOT-8783 | General purpose HIJING 50-90%            | runsMC-2022-LHC22b5-pass3.list          | 10M                   |
| LHC20d2abc  | Pb–Pb    | LHC18qr | ALIROOT-8423 | HIJING plus injected nuclei, hypernuclei | runsMC-2020-LHC20d2[a,b,c]-pass3.list   | 120k, 400k, 400k      |
| LHC20g4     | pp 5TeV  |         | ALIROOT-8505 | jet-jet, no ESDs                         |                                         |                       |
| LHC18g4     | pp 13TeV | LHC18b  | ALIROOT-7903 | General purpose PYTHIA8                  | runsMC-2018-LHC18g4-pass1.list          |                       |

# data location
| Production  | alien location                                                | #files | 
| ----------- | ------------------------------------------------------------- | ------ |
| LHC20e3a    | alien.py find alien:///alice/sim/2020/LHC20e3a   AliESDs.root | 398531 | 
| LHC20e3b    | alien.py find alien:///alice/sim/2020/LHC20e3b   AliESDs.root | 106955 | 
| LHC20e3c    | alien.py find alien:///alice/sim/2020/LHC20e3c   AliESDs.root | 551584 | 
| LHC22b5     | alien.py find alien:///alice/sim/2022/LHC22b5    AliESDs.root | 719926 | 
| LHC20d2a    | alien.py find alien:///alice/sim/2020/LHC20d2a   AliESDs.root |  63615 |
| LHC20d2b    | alien.py find alien:///alice/sim/2020/LHC20d2b   AliESDs.root |  59542 | 
| LHC20d2c    | alien.py find alien:///alice/sim/2020/LHC20d2c   AliESDs.root |   4583 |
| LHC20g4     | alien.py find alien:///alice/sim/2020/LHC20g4    AliESDs.root |      0 | 
| LHC18g4     | alien.py find alien:///alice/sim/2018/LHC18g4    AliESDs.root | 132421 |
| LHC22d1c2   | alien.py find alien:///alice/sim/2022/LHC22d1c2  galice.root  | 116637 |

# How to run runGrid.c
cd $RUN_ON_GRID_DIR/Ebye/test_fastGen_netParticles/
cp $RUN_ON_GRID_DIR/Ebye/code/*.* $RUN_ON_GRID_DIR/Ebye/code/READ*  .; rm *.so *.d

# test lxplus
cp $RUN_ON_GRID_DIR/Ebye/code/AliAnalysisTaskTIdentityPID.*   /afs/cern.ch/work/m/marsland/workdir/test_pp_MC
cp $RUN_ON_GRID_DIR/Ebye/code/AliAnalysisTaskTIdentityPID.*   /afs/cern.ch/work/m/marsland/workdir/test_pp_MC
cp $RUN_ON_GRID_DIR/Ebye/code/AddTask_marsland_TIdentityPID.C /afs/cern.ch/work/m/marsland/workdir/test_pp_MC
cp $RUN_ON_GRID_DIR/Ebye/code/Config_marsland_TIdentityPID.C  /afs/cern.ch/work/m/marsland/workdir/test_pp_MC
cp $RUN_ON_GRID_DIR/Ebye/code/README_2024.md                  /afs/cern.ch/work/m/marsland/workdir/test_pp_MC
cp $RUN_ON_GRID_DIR/Ebye/code/runGrid.C                       /afs/cern.ch/work/m/marsland/workdir/test_pp_MC
cp $RUN_ON_GRID_DIR/Ebye/code/AddTaskFilteredTreeLocal.C      /afs/cern.ch/work/m/marsland/workdir/test_pp_MC

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

# Fill Eff matrix
--> PbPb
aliroot -b -q 'runGrid.C(0,"test",0,  40, "$TIdentityDIRcommit/TIdentity/lists/runsMC-2020-LHC20e3a-pass3.list",1,2018,"18q",3,0)'
aliroot -b -q 'runGrid.C(0,"test",0,  40, "$TIdentityDIRcommit/TIdentity/lists/runsMC-2020-LHC20e3b-pass3.list",1,2018,"18q",3,0)'
aliroot -b -q 'runGrid.C(0,"test",0,  40, "$TIdentityDIRcommit/TIdentity/lists/runsMC-2020-LHC20e3c-pass3.list",1,2018,"18q",3,0)'
aliroot -b -q 'runGrid.C(0,"test",0,  40, "$TIdentityDIRcommit/TIdentity/lists/runsMC-2022-LHC22b5-pass3.list",1,2018,"18q",3,0)'

# Run real data
aliroot -b -q 'runGrid.C(0,"test",0,  0,  "$TIdentityDIRcommit/TIdentity/lists/runs-2018-LHC18q-pass3.list",0,2018,"18q",3,0)'
aliroot -b -q 'runGrid.C(0,"test",0,  0,  "$TIdentityDIRcommit/TIdentity/lists/runs-2018-LHC18r-pass3.list",0,2018,"18r",3,0)'
aliroot -b -q 'runGrid.C(0,"test",0,  11, "$TIdentityDIRcommit/TIdentity/lists/runs-2018-LHC18b-pass2.list",0,2018,"18b",3,0)'


# run full MC
--> pp
aliroot -b -q 'runGrid.C(0,"test",0,  51, "$TIdentityDIRcommit/TIdentity/lists/runsMC-2018-LHC18g4-pass1.list",1,2018,"18b",3,0)'

--> PbPb
aliroot -b -q 'runGrid.C(0,"test",0,  50, "$TIdentityDIRcommit/TIdentity/lists/runsMC-2020-LHC20e3a-pass3.list",1,2018,"18q",3,0)'
aliroot -b -q 'runGrid.C(0,"test",0,  50, "$TIdentityDIRcommit/TIdentity/lists/runsMC-2020-LHC20e3b-pass3.list",1,2018,"18q",3,0)'
aliroot -b -q 'runGrid.C(0,"test",0,  50, "$TIdentityDIRcommit/TIdentity/lists/runsMC-2020-LHC20e3c-pass3.list",1,2018,"18q",3,0)'
aliroot -b -q 'runGrid.C(0,"test",0,  50, "$TIdentityDIRcommit/TIdentity/lists/runsMC-2022-LHC22b5-pass3.list",1,2018,"18q",3,0)'

aliroot -b -q 'runGrid.C(0,"test",0,  50, "$TIdentityDIRcommit/TIdentity/lists/runsMC-2020-LHC20k6a-pass3.list",1,2018,"18q",3,0)'
aliroot -b -q 'runGrid.C(0,"test",0,  50, "$TIdentityDIRcommit/TIdentity/lists/runsMC-2020-LHC20k6b-pass3.list",1,2018,"18q",3,0)'
aliroot -b -q 'runGrid.C(0,"test",0,  50, "$TIdentityDIRcommit/TIdentity/lists/runsMC-2020-LHC20k6c-pass3.list",1,2018,"18q",3,0)'
aliroot -b -q 'runGrid.C(0,"test",0,  50, "$TIdentityDIRcommit/TIdentity/lists/runsMC-2020-LHC20k6d-pass3.list",1,2018,"18q",3,0)'



aliroot -b -q 'runGrid.C(0,"test",0,  51, "$TIdentityDIRcommit/TIdentity/lists/runsMC-2018-LHC18g4-pass1.list",1,2018,"18b",3,0)'

# run fastGen
aliroot -b -q 'runGrid.C(0,"test",0,  200, "$TIdentityDIRcommit/TIdentity/lists/runsGen-2022-LHC22d1c2-pass3.list",2,2022,"22d1c2",2,0)'
aliroot -b -q 'runGrid.C(0,"test",0,  200, "$TIdentityDIRcommit/TIdentity/lists/runsGen-2022-LHC22d1d2-pass3.list",2,2022,"22d1d2",2,0)'

# copy data
alien_cp -T 6 -parent 99 -glob AnalysisResults.root /Alice/cern.ch/user/m/marsland/PWGPP695_MC_remapping/LHC20e3a_pass3_20240625_228/2020 file:

# merge data e.g. only for ebye fluct. related objects
alihadd -i "cleanHists" -i "mcFull" -i "mcGen" -i "fTreeMC" -i "eventInfo" -i "eventInfoMC" -i "cutBased" -s 2000000000 AnalysisResults_mcTrees.root  @file.list

# kill the list of jobs 
for i in $(cat jobs.list); do alien.py kill $i; done

# resubmit all jobs of LHC20e3a in error --> see manual of "ps" in https://jalien.docs.cern.ch/jalien_commands/
for i in $(alien.py ps -E | grep LHC20e3a | awk '{print $2}'); do alien.py resubmit $i; done

# kill all masterjobs for a given period
for i in $(alien.py ps -M | grep LHC18b   | awk '{print $2}'); do alien.py kill $i; done
for i in $(alien.py ps -M | grep LHC18g4  | awk '{print $2}'); do alien.py kill $i; done
for i in $(alien.py ps -M | grep LHC20e3a | awk '{print $2}'); do alien.py kill $i; done






