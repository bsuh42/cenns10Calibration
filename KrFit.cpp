#include <iostream>
#include <string>
using namespace std;
#include <cmath>
#include <vector>

double bottomSPEFit(Double_t timestamp)
{
  //For fit SPE method
  const Double_t p0 = 178.43;
  const Double_t p1 = -0.72e-7;
  Double_t value = p0+p1*timestamp;
  return value;
}

double bottomSPESpline(Double_t timestamp, vector<Double_t> splineTimestamp, vector<Double_t> ch3Spline)
{
  Double_t readValue = splineTimestamp[0];
  Double_t value = ch3Spline[0];
  Int_t counter = 0;
  while (timestamp >= readValue)
  {
    counter += 1;
    readValue = splineTimestamp[counter];
    value = ch3Spline[counter];
  }
  return value;
}

double topSPESpline(Double_t timestamp, vector<Double_t> splineTimestamp, vector<Double_t> ch1Spline)
{
  Double_t readValue = splineTimestamp[0];
  Double_t value = ch1Spline[0];
  Int_t counter = 0;
  while (timestamp >= readValue)
  {
    counter += 1;
    readValue = splineTimestamp[counter];
    value = ch1Spline[counter];
  }
  return value;
}

void fitKr(TString KrFile, const Int_t injectionRun, const Double_t lowerBound, const Double_t upperBound)
{
  //injectionRun = runNumber when injection was started
  //Fit a krypton peak. Must define when injection run started to determine background
  TFile *FKr = new TFile(KrFile);

  TTree *TKr = (TTree *)FKr->Get("compressedTree");

  //Variables
  Double_t integralCh1 = 0;
  Double_t integralCh3 = 0;
  Double_t totalIntegral = 0;
  Int_t runNumber = 0;
  Double_t baselineCh1 = 0;
  Double_t baselineCh3 = 0;
  Double_t timestamp = 0;
  Double_t temp = 0;
  Int_t counter = 0;
  Int_t evtNumber = 0;
  Int_t index = 0;


  //Read the Spline file
  vector<Double_t> splineTimestamp;
  vector<Double_t> ch1Spline;
  vector<Double_t> ch3Spline;
  Double_t spetime;
  Double_t spe1;
  Double_t spe3;
  ifstream myFile;
  myFile.open("SPE_Calibration_TimeTable_Ch1Inferred_AvgRatioValue.txt");
  string line;
  while (myFile)
  {
    getline(myFile, line);
    stringstream s(line);
    s >> spetime >> spe1 >> spe3;
    splineTimestamp.push_back(spetime);
    ch1Spline.push_back(spe1);
    ch3Spline.push_back(spe3);
  }

  //Create histogram
  const Int_t histStart = 0;
  const Int_t histEnd = 1200;
  const Int_t numberBins = 200;
  TH1D* hist = new TH1D("Kr", "Kr", numberBins, histStart, histEnd);
  TH1D* unNormHist = new TH1D("", "", numberBins, histStart, histEnd);
  TH1D* bkgHist = new TH1D("bkg", "bkg", numberBins, histStart, histEnd);
  vector<Double_t> timeStart;
  vector<Double_t> timeEnd;
  vector<Double_t> numberOfEvents;
  vector<Int_t> KrRunNumbers;

  //Set branches to get
  TKr->SetBranchAddress("runNumber", &runNumber);
  TKr->SetBranchAddress("integralCh1", &integralCh1);
  TKr->SetBranchAddress("integralCh3", &integralCh3);
  TKr->SetBranchAddress("baselineCh1", &baselineCh1);
  TKr->SetBranchAddress("baselineCh3", &baselineCh3);
  TKr->SetBranchAddress("timestamp", &timestamp);

  Int_t numberEntries = TKr->GetEntries();
  for (Int_t eventNumber = 0; eventNumber < numberEntries; eventNumber++)
  {
    totalIntegral = 0;
    index = 0;
    TKr->GetEntry(eventNumber);
    if (eventNumber % 1000 == 0)
    {
      printf("Current on event %d of %d\n", eventNumber, numberEntries);
    }
    std::vector<int>::iterator it = std::find(KrRunNumbers.begin(), KrRunNumbers.end(), runNumber);
    if (integralCh1 != 0 && integralCh3 != 0 && baselineCh1 != 0 && baselineCh3 != 0)
    {
      if (it != KrRunNumbers.end())
      {
        index = std::distance(KrRunNumbers.begin(), it);
      }
      else
      {
        KrRunNumbers.push_back(runNumber);
	timeStart.push_back(timestamp);
	timeEnd.push_back(timestamp);
	numberOfEvents.push_back(0);
	index = 0;
      }
      //timing
      if (timestamp < timeStart[index])
      {
        timeStart[index] = timestamp;
      }
      if (timestamp > timeEnd[index])
      {
        timeEnd[index] = timestamp;
      }

      temp = topSPESpline(timestamp, splineTimestamp, ch1Spline);
      totalIntegral = integralCh1/temp;
      temp = bottomSPESpline(timestamp, splineTimestamp, ch3Spline);
      totalIntegral += integralCh3/temp;

      hist->Fill(totalIntegral);
      unNormHist->Fill(totalIntegral);
      hist->Fill(totalIntegral);
      unNormHist->Fill(totalIntegral);
      numberOfEvents[index] = numberOfEvents[index] + 1;
    }
  }

  //Normalize by time
  Double_t timeLength = 0;
  vector<Double_t> timeSince;
  Double_t totalTime = 0;
  for (Int_t runIndex = 0; runIndex < KrRunNumbers.size(); runIndex++)
  {
    timeLength = (timeEnd[runIndex]-timeStart[runIndex]);
    numberOfEvents[runIndex] = (numberOfEvents[runIndex]*1.0)/(timeLength*1.0);
    timeSince.push_back(timeStart[runIndex]-timeStart[0]);
    totalTime += timeLength;
  }
  hist->Scale(1.0/totalTime);

  Double_t center;
  Int_t binmax;
  //hist->GetXaxis()->SetRange(lowerBound, upperBound);
  binmax = hist->GetMaximumBin();
  center = hist->GetXaxis()->GetBinCenter(binmax);
  center = 200; //Define the center around where the peak is expected

  Double_t lowerLower = center - 110;
  Double_t upperLower = center - 10;
  Double_t lowerUpper = center + 10;
  Double_t upperUpper = center + 110;
  Double_t stepSize = 5;
  Int_t numberOfLower = ((upperLower-lowerLower)/stepSize)+1;
  Int_t numberOfUpper = ((upperUpper-lowerUpper)/stepSize)+1;
  Int_t totalNumber = numberOfLower*numberOfUpper;
  Double_t fitMean;
  Double_t fitError;
  Double_t minLower;
  Double_t minUpper;
  Double_t minMean;
  Double_t minError = 100;
  Double_t minChi2PerNDF;
  Double_t minChi2Desired = 100;
  Int_t minCounter;
  Double_t chi2 = 0;
  Double_t NDF = 0;
  Double_t chi2PerNDF = 0;
  Double_t chi2Desired = 0;
  Double_t tryLower = 0;
  Double_t tryUpper = 0;
  Double_t minChi2 = 0;
  Double_t minNDF;
  counter = 0;

  TF1** unNormFit = new TF1*[totalNumber];
  TF1** gausFit = new TF1*[totalNumber];

  for (Int_t lower = 0; lower < numberOfLower; lower++)
  {
    for (Int_t upper = 0; upper < numberOfUpper; upper++)
    {
      tryLower = lowerLower + stepSize*lower;
      tryUpper = lowerUpper + stepSize*upper;

      unNormFit[counter] = new TF1("", "gaus", tryLower, tryUpper);
      gausFit[counter] = new TF1("gausFit", "gaus", tryLower, tryUpper);
      unNormHist->Fit(unNormFit[counter], "QR");

      gausFit[counter]->SetParameter(1, unNormFit[counter]->GetParameter(1));
      gausFit[counter]->SetParameter(2, unNormFit[counter]->GetParameter(2));
      hist->Fit(gausFit[counter], "QR");

      fitMean = gausFit[counter]->GetParameter(1);
      fitError = gausFit[counter]->GetParError(1);
      chi2 = gausFit[counter]->GetChisquare();
      NDF = gausFit[counter]->GetNDF();
      chi2PerNDF = chi2/NDF;
      chi2Desired = abs(chi2PerNDF-1.0);
      if (chi2Desired <= minChi2Desired && fitMean >= tryLower && fitMean <= tryUpper)
      {
        minLower = tryLower;
	minUpper = tryUpper;
	minMean = fitMean;
	minError = fitError;
	minCounter = counter;
	minChi2Desired = chi2Desired;
	minChi2PerNDF = chi2PerNDF;
	minChi2 = chi2;
	minNDF = NDF;
      }
      counter += 1;
    }
  }
  printf("lower bound = %f, upper bound = %f, mean = %f, sigma = %f, chi2/DoF = %f\n", minLower, minUpper, minMean, minError, minChi2PerNDF);

  TCanvas *c1 = new TCanvas("c1", "", 1000, 600);
  TCanvas *c2 = new TCanvas("c2", "", 1000, 600);

  c1->cd();
  gStyle->SetOptStat(0);
  c1->SetLogy();
  hist->Draw("hist");
  gausFit[minCounter]->Draw("SAME");
  c1->Update();

  c2->cd();
  gStyle->SetOptStat(0);
  TGraph *gr = new TGraph(numberOfEvents.size(), &timeSince[0], &numberOfEvents[0]);
  gr->Draw("AP");
  c2->Update();
}
