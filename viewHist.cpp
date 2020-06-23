#include <iostream>
#include <string>
using namespace std;
#include <cmath>
#include <vector>

void viewHist(TString file)
{
  //View a stored histogram
  TFile *f = new TFile(file);

  TString newFile= file;
  TString extension = ".root";
  newFile = newFile.Remove(newFile.Length()-extension.Length(), extension.Length());
  TH1D* h1 = new TH1D(newFile, newFile, 400, 900, 500);
  h1 = (TH1D*)f->Get("hist");
  h1->GetXaxis()->SetTitle("PE");
  h1->GetYaxis()->SetTitle("Counts");
  
  TCanvas *c1 = new TCanvas("c1", "", 1000, 600);
  gStyle->SetOptStat(0);
  h1->Draw();
  c1->Update();
}
