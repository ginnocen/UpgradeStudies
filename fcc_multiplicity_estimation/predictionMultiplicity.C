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

  canvas = new TCanvas("c", "c", 800, 600);
  gPad->SetRightMargin(0.05);
  gPad->SetTopMargin(0.05);
  
  
  dummy = new TH2F("dummy", ";N_{ch};P(N_{ch})", 100, 0.1, 30000, 100, 1e-7, 50.);
  dummy->SetStats(0);
  dummy->GetYaxis()->SetTitleOffset(1.2);
  dummy->Draw();
  gPad->SetLogy();
  legend = new TLegend(0.58, 0.6, 0.91, 0.92);
  legend->SetFillColor(0);
  legend->SetLineColor(0);
  
  
  TH1F*hpPb5TeV=(TH1F*)fillParticleDistribution("TGlauberMC-2.4/glau_ppb_smeared_ntuple.root");
  TH1F*hPbPb5TeV=(TH1F*)fillParticleDistribution("TGlauberMC-2.4/glau_pbpb_smeared_ntuple.root");
  TH1F*hXeXe5TeV=(TH1F*)fillParticleDistribution("TGlauberMC-2.4/glau_xexe_smeared_ntuple.root");
  hpPb5TeV->Draw("same");
  hPbPb5TeV->Draw("same");
  hXeXe5TeV->Draw("same");
   

  TFile*fout=new TFile("outputALICE7TeV.root");
  TGraphErrors*graphALICE7TeV=(TGraphErrors*)fout->Get("graphALICE7TeV");
  graphALICE7TeV->Draw("same");
  
  func = GetNBDLog("nbd_7");
  func->SetParameters(1, 20.8, k->Eval(7000));
  func->SetRange(0.3 * 151, 170);
  func->SetLineColor(2);
// 	func->SetLineStyle(3);
  func->Draw("SAME");


  
}

TH1F* fillParticleDistribution(TString filename){

  func = GetNBDLog("nbd_5");
  func->SetParameters(1, 20.8, k->Eval(5000));
  func->SetRange(0, 170);

  TFile* fin = TFile::Open(filename.Data()); 
  TTree* mytree = (TTree*)fin->Get("nt");
  //TTree* mytree = (TTree*)fin->Get("HiTree");
  
  TH1F*histo=new TH1F("histo","histo",3000,0,30000);
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