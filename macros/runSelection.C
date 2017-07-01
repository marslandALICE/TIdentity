/*
 *  .L $NOTES/JIRA/PWGPP-340/runSelection.C
 */

void runSelection(const Int_t year, const TString period, const TString pass){
   
    /*
     
    cd /home/marsland/Desktop/RUN_ON_GRID/Ebye/lists
    aliroot -l 
    .L /home/marsland/Desktop/RUN_ON_GRID/Ebye/code/runSelection.C
    runSelection(2015,"LHC15o","pass2_lowIR")
    
    */
    
    TString runListSelected = Form("runs-%d-%s-%s.list",year,period.Data(),pass.Data());
    TString runListRejected = Form("runsRejected-%d-%s-%s.list",year,period.Data(),pass.Data());
   
    AliExternalInfo info;
    TTree * treeLogbook = info.GetTree("Logbook",period,pass,"QA.TPC;QA.TRD;QA.TOF;QA.ITS;QA.EVS;Logbook.detector:TPC:detector==\"TPC\";MonALISA.RCT");
    TTree * treeTPC = info.GetTree("QA.TPC",period,pass,"QA.TRD;QA.TPC;QA.ITS;QA.TOF;Logbook;QA.EVS;Logbook.detector:TPC:detector==\"TPC\";MonALISA.RCT");
    
    //
    // example PID, DCAR and DCA z within +-6 sigma
    //             PID_Status<100&&DCAz_Status<100&&DCAr_Status<100"
    
    treeTPC->Scan("run:PID_Status:DCAz_Status:DCAr_Status","PID_Status<100&&DCAz_Status<100&&DCAr_Status<100");
    // example PID, DCAR and DCA z within +-6 sigma + eff>0.5
    entries = treeTPC->Draw("EffTOTPt1","1");
    treeTPC->SetAlias("medEff",Form("(%f+0)",TMath::Median(entries,treeTPC->GetV1())));
    treeTPC->Scan("run:PID_Status:DCAz_Status:DCAr_Status:EffTOTPt1","PID_Status<100&&DCAz_Status<100&&DCAr_Status<100&&EffTOTPt1>medEff-0.1");
    //
    // Make run list as an csv file
    //
    //     AliTreePlayer::selectWhatWhereOrderBy(treeTPC, "run:PID_Status:DCAz_Status:DCAr_Status:EffTOTPt1:MonALISA.RCT.tpc_value","PID_Status<100&&DCAz_Status<100&&DCAr_Status<100&&EffTOTPt1>medEff-0.1","",0,10000,"csv",runListSelected);
    AliTreePlayer::selectWhatWhereOrderBy(treeTPC, "run:PID_Status:DCAz_Status:DCAr_Status:EffTOTPt1","PID_Status<100&&DCAz_Status<100&&DCAr_Status<100&&EffTOTPt1>medEff-0.1","",0,10000,"csv",runListSelected);
    //
    // Make rejected run list
    //
    //     AliTreePlayer::selectWhatWhereOrderBy(treeTPC, "run:PID_Status:DCAz_Status:DCAr_Status:EffTOTPt1:MonALISA.RCT.tpc_value","!(PID_Status<100&&DCAz_Status<100&&DCAr_Status<100&&EffTOTPt1>medEff-0.1)","",0,10000,"csv",runListRejected);
    AliTreePlayer::selectWhatWhereOrderBy(treeTPC, "run:PID_Status:DCAz_Status:DCAr_Status:EffTOTPt1","!(PID_Status<100&&DCAz_Status<100&&DCAr_Status<100&&EffTOTPt1>medEff-0.1)","",0,10000,"csv",runListRejected);
    //
    //
    //
    //     AliSysInfo::PrintJiraTable(treeTPC, "run:PID_Status:DCAz_Status:DCAr_Status:EffTOTPt1:MonALISA.RCT.tpc_value","!(PID_Status<100&&DCAz_Status<100&&DCAr_Status<100&&EffTOTPt1>medEff-0.1)","colsize=120","runRejected.jira");
    //     AliSysInfo::PrintJiraTable(treeTPC, "run:PID_Status:DCAz_Status:DCAr_Status:EffTOTPt1:totalEventsPhysics:Logbook.detector_TPC.eventCountPhysics:MonALISA.RCT.tpc_value","!(PID_Status<100&&DCAz_Status<100&&DCAr_Status<100&&EffTOTPt1>medEff-0.1)","colsize=30","runRejected.jira");
    
}
