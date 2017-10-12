
TF1* GetNBDLog(const char* name = "nbd")
{
        TF1* func = new TF1(name, "exp(log([0]) + TMath::LnGamma([2]+x) - TMath::LnGamma([2]) - TMath::LnGamma(x+1) + log([1] / ([1]+[2])) * x + log(1.0 + [1]/[2]) * -[2])");
        
        func->SetParNames("scaling", "averagen", "k");
        func->SetParLimits(0, 0.5, 2);
        func->SetParLimits(1, 1, 100);
        func->SetParLimits(2, 1, 20);
        func->SetParameters(1, 10, 2);        
        return func;
}

void predictionMultiplicity(){

  k = new TF1("k", "1.94", 1, 1e5); 

  TCanvas*mycanvas = new TCanvas("mycanvas", "mycanvas", 800, 800);
  gPad->SetRightMargin(0.05);
  gPad->SetTopMargin(0.10);
  
  dummy = new TH2F("dummy", ";N_{ch};P(N_{ch})", 100, 0.1, 22000, 100, 1e-7, 50.);
  dummy->GetXaxis()->SetNdivisions(8,true);
  dummy->GetYaxis()->SetNdivisions(8,true);
  dummy->SetStats(0);
  dummy->GetYaxis()->SetLabelSize(0.04);
  dummy->GetYaxis()->SetLabelSize(0.04);
  dummy->Draw("same");
  gPad->SetLogy();
  legend = new TLegend(0.58, 0.6, 0.91, 0.92);
  legend->SetFillColor(0);
  legend->SetLineColor(0);
  
  TLatex* texCms = new TLatex(0.18,0.93, "#scale[1.25]{CMS} Performance");
  texCms->SetNDC();
  texCms->SetTextAlign(12);
  texCms->SetTextSize(0.04);
  texCms->SetTextFont(42);
  texCms->Draw();


  TFile*fout=new TFile("outputALICE7TeV.root");
  TGraphErrors*graphALICE7TeV=(TGraphErrors*)fout->Get("graphALICE7TeV");

  TH1F*hpPb5TeV=(TH1F*)fillParticleDistribution("TGlauberMC-2.4/glau_ppb_smeared_ntuple.root");
  TH1F*hPbPb5TeV=(TH1F*)fillParticleDistribution("TGlauberMC-2.4/glau_pbpb_smeared_ntuple.root");
  TH1F*hXeXe5TeV=(TH1F*)fillParticleDistribution("TGlauberMC-2.4/glau_xexe_smeared_ntuple.root");
  graphALICE7TeV->SetMarkerColor(1);
  graphALICE7TeV->SetLineColor(1);
  graphALICE7TeV->SetMarkerStyle(28);
  graphALICE7TeV->Draw("psame");

  hpPb5TeV->SetMarkerColor(kRed+1);
  hpPb5TeV->SetLineColor(kRed+1);
  //hpPb5TeV->SetMarkerStyle(24);
  hpPb5TeV->SetLineWidth(2);
  hpPb5TeV->Draw("same");
  
  hXeXe5TeV->SetMarkerColor(kGreen+3);
  hXeXe5TeV->SetLineColor(kGreen+3);
  //hXeXe5TeV->SetMarkerStyle(26);
  hXeXe5TeV->SetLineWidth(2);
  hXeXe5TeV->Draw("same");

  hPbPb5TeV->SetMarkerColor(kBlue+2);
  hPbPb5TeV->SetLineColor(kBlue+2);
  //hPbPb5TeV->SetMarkerStyle(25);
  hPbPb5TeV->SetLineWidth(2);
  hPbPb5TeV->Draw("same");  
  
  TLegend *legendSigma=new TLegend(0.233871,0.6419492,0.5322581,0.815678,"");
  legendSigma->SetFillColor(0);
  legendSigma->SetLineColor(0);
  
  TLegendEntry *entpp=legendSigma->AddEntry(graphALICE7TeV,"ALICE pp 7 TeV","P");
  entpp->SetTextFont(42);
  entpp->SetTextColor(1);
  entpp->SetLineColor(1);
  entpp->SetMarkerColor(1);
  legendSigma->Draw();
  
  TLegendEntry *entpPb=legendSigma->AddEntry(hpPb5TeV,"NBD pPb 5 TeV","PL");
  entpPb->SetTextFont(42);
  entpPb->SetTextColor(kRed+1);
  entpPb->SetLineColor(kRed+1);
  entpPb->SetMarkerColor(kRed+1);
  legendSigma->Draw();

  TLegendEntry *entXeXe=legendSigma->AddEntry(hXeXe5TeV,"NBD XeXe 5 TeV","PL");
  entXeXe->SetTextFont(42);
  entXeXe->SetTextColor(kGreen+3);
  entXeXe->SetLineColor(kGreen+3);
  entXeXe->SetMarkerColor(kGreen+3);
  legendSigma->Draw();

  TLegendEntry *entPbPb=legendSigma->AddEntry(hPbPb5TeV,"NBD PbPb 5 TeV","PL");
  entPbPb->SetTextFont(42);
  entPbPb->SetTextColor(kBlue+2);
  entPbPb->SetLineColor(kBlue+2);
  entPbPb->SetMarkerColor(kBlue+2);
  legendSigma->Draw();

  TPad* pzoomhnTrkvshiBin = new TPad("pzoomhnTrkvshiBin","",0.60, 0.60, 0.95, 0.90);
  pzoomhnTrkvshiBin->SetLogy();
  pzoomhnTrkvshiBin->SetLogx();
  pzoomhnTrkvshiBin->SetFillColor(0);
  pzoomhnTrkvshiBin->SetBorderMode(0);
  pzoomhnTrkvshiBin->SetBorderSize(2);
  pzoomhnTrkvshiBin->SetRightMargin(0);
  pzoomhnTrkvshiBin->SetTopMargin(0);
  pzoomhnTrkvshiBin->Draw();
  pzoomhnTrkvshiBin->cd();

  TH2F* hemptyzoomhnTrkvshiBin = new TH2F("hemptyzoomhnTrkvshiBin","",1000, 0, 3500, 100, 1e-7, 50);
  hemptyzoomhnTrkvshiBin->GetXaxis()->SetNdivisions(4,true);
  hemptyzoomhnTrkvshiBin->GetYaxis()->SetNdivisions(5,true);
  hemptyzoomhnTrkvshiBin->GetXaxis()->SetLabelSize(0.12);
  hemptyzoomhnTrkvshiBin->GetYaxis()->SetLabelSize(0.12);
  hemptyzoomhnTrkvshiBin->Draw();
  hpPb5TeV->Draw("same");
  //hpPb5TeV->Draw("same");
  graphALICE7TeV->Draw("pesame");
  
}

TH1F* fillParticleDistribution(TString filename){

  func = GetNBDLog("nbd_5");
  func->SetParameters(1, 20.8, k->Eval(5000));
  func->SetRange(0, 170);

  TFile* fin = TFile::Open(filename.Data()); 
  TTree* mytree = (TTree*)fin->Get("nt");
  //TTree* mytree = (TTree*)fin->Get("HiTree");
  
  TH1F*histo=new TH1F("histo","histo",500,0,30000);
  double f_factor_5TeV=0.8;
  int nparticle;
  float Npart, Ncoll;

  mytree->SetBranchAddress("Npart", &Npart);
  mytree->SetBranchAddress("Ncoll", &Ncoll);
  
  int nancestors=0;
  Long64_t nentries = mytree->GetEntries();
  
  for(Long64_t i=0;i<nentries;i++){
    mytree->GetEntry(i);
    nancestors=(int) (Npart*f_factor_5TeV+Ncoll*(1-f_factor_5TeV));
    nparticle=0;
    for (int j=0;j<nancestors;j++){
      nparticle=nparticle+func->GetRandom();
    }
    histo->Fill(nparticle);
  }
  histo->Scale(1./histo->GetEntries());

  return histo;
}