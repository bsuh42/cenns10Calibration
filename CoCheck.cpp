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
    //Look for closest value to timestamp
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
    //Look for closest value to timestamp
    counter += 1;
    readValue = splineTimestamp[counter];
    value = ch1Spline[counter];
  }
  return value;
}

void fitCo(TString coFile, TString bkgFile, TString SPEMethod, const Int_t coDesiredRun, const Double_t lowerBound, const Double_t upperBound)
{
  //Fit a single co peak given a certain fit range.
  //SPEMethod can be fit, fixed, or spline
  //Should find bkg run automatically.
  //Does not find the best range. 
  //Saves output hist
  //For use with genTreeCal.cpp+

  //Load files
  TFile *FCo = new TFile(coFile);
  TFile *FBkg = new TFile(bkgFile);

  TTree *TCo = (TTree *)FCo->Get("compressedTree");
  TTree *TBkg = (TTree *)FBkg->Get("compressedTree");

  //Variables
  Double_t integralCh1 = 0; //output of top detector
  Double_t integralCh3 = 0; //output of bottom detecor
  Double_t totalIntegral = 0;
  Int_t runNumber = 0; 
  Double_t baselineCh1 = 0; //automatically found baseline, top detector
  Double_t baselineCh3 = 0; //automatically found baseline, bottom detector
  Double_t timestamp = 0;
  Double_t timeStartCo = 0;
  Double_t timeStartBkg = 0;
  Double_t timeEndCo = 0;
  Double_t timeEndBkg = 0; //want to store time length of run to normalize
  Int_t bkgDesiredRun = 0;
  
  Double_t temp = 0;
  Int_t counter = 0;
  Int_t index = 0;

  //SPE variables
  const Double_t speFactor = 1.18; //Use a fixed conversion from spe3 to spe1. Divide spe3 value by this factor. Not for Spline
  const Double_t speValue3 = 64.0; //fixed SPE
  //Spline
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

  //Set branches to get
  TCo->SetBranchAddress("runNumber", &runNumber);
  TCo->SetBranchAddress("integralCh1", &integralCh1);
  TCo->SetBranchAddress("integralCh3", &integralCh3);
  TCo->SetBranchAddress("baselineCh1", &baselineCh1);
  TCo->SetBranchAddress("baselineCh3", &baselineCh3);
  TCo->SetBranchAddress("timestamp", &timestamp);

  TBkg->SetBranchAddress("runNumber", &runNumber);
  TBkg->SetBranchAddress("integralCh1", &integralCh1);
  TBkg->SetBranchAddress("integralCh3", &integralCh3);
  TBkg->SetBranchAddress("baselineCh1", &baselineCh1);
  TBkg->SetBranchAddress("baselineCh3", &baselineCh3);
  TBkg->SetBranchAddress("timestamp", &timestamp);

  //Histogram details
  const Int_t histStart = 400;
  const Int_t histEnd = 900; //might need to change based on what region you're looking at
  const Int_t numberBins = 50;
  TH1D *CoHist = new TH1D(TString::Format("Run%d", coDesiredRun), TString::Format("Run%d", coDesiredRun), numberBins, histStart, histEnd);
  TH1D *BkgHist = new TH1D("", "", numberBins, histStart, histEnd);
  TH1D *UnNormHist = new TH1D("", "", numberBins, histStart, histEnd); //for getting fit parameters

  //Fill CoHist
  const Int_t CoNumberEntries = TCo->GetEntries();
  for (Int_t eventNumber = 0; eventNumber < CoNumberEntries; eventNumber++)
  {
    totalIntegral = 0;
    if (eventNumber % 1000 == 0)
    {
      printf("Co. Currently on event %d of %d\n", eventNumber, CoNumberEntries);
    }
    TCo->GetEntry(eventNumber);

    if (runNumber == coDesiredRun)
    {
      if (integralCh1 != 0 && integralCh3 != 0 && baselineCh1 != 0 && baselineCh3 != 0)
      {
        //Consider a good event where something is written and nothing funky happened to the baseline
	//Timing
	if (timeStartCo == 0)
	{
	  timeStartCo = timestamp;
	}
	if (timestamp <= timeStartCo)
	{
	  timeStartCo = timestamp;
	}
	if (timestamp >= timeEndCo)
	{
	  timeEndCo = timestamp;
	}

	//Calculate integral based on SPE method
	if (SPEMethod == "fit")
	{
	  temp = bottomSPEFit(timestamp);
	  totalIntegral = integralCh1/(temp/speFactor);
	  totalIntegral += integralCh3/temp;
	}
	if (SPEMethod == "fixed")
	{
	  totalIntegral = integralCh1/(speValue3/speFactor);
	  totalIntegral += integralCh3/speValue3;
	}
	if (SPEMethod == "spline")
	{
	  temp = bottomSPESpline(timestamp, splineTimestamp, ch3Spline);
	  totalIntegral = integralCh3/temp;
	  temp = topSPESpline(timestamp, splineTimestamp, ch1Spline);
	  totalIntegral += integralCh1/temp;
	}

	//Fill histograms
	CoHist->Fill(totalIntegral);
	UnNormHist->Fill(totalIntegral);
      }
    }
  }

  //Do the same for bkg
  const Int_t BkgNumberEntries = TBkg->GetEntries();
  for (Int_t eventNumber = 0; eventNumber < BkgNumberEntries; eventNumber++)
  {
    totalIntegral = 0;
    if (eventNumber % 1000 == 0)
    {
      printf("Bkg. Currently on event %d of %d\n", eventNumber, BkgNumberEntries);
    }
    TBkg->GetEntry(eventNumber);

    //search within 15 runs in either direction for corresponding bkg run
    if (runNumber >= coDesiredRun - 20 && runNumber <= coDesiredRun + 20 && bkgDesiredRun == 0)
    {
      bkgDesiredRun = runNumber;
    }

    //if the background run has been found, calculate total integral
    if (runNumber == bkgDesiredRun)
    {
      if (integralCh1 != 0 && integralCh3 != 0 && baselineCh1 != 0 && baselineCh3 != 0)
      {
        if (timeStartBkg == 0)
	{
	  timeStartBkg = timestamp;
	}
	if (timestamp <= timeStartBkg)
	{
	  timeStartBkg = timestamp;
	}
	if (timestamp >= timeEndBkg)
	{
	  timeEndBkg = timestamp;
	}

	//Calculate integral
	if (SPEMethod == "fit")
	{
	  temp = bottomSPEFit(timestamp);
	  totalIntegral = integralCh1/(temp/speFactor);
	  totalIntegral += integralCh3/temp;
	}
	if (SPEMethod == "fixed")
	{
	  totalIntegral = integralCh1/(speValue3/speFactor);
	  totalIntegral += integralCh3/speValue3;
	}
	if (SPEMethod == "spline")
	{
	  temp = bottomSPESpline(timestamp, splineTimestamp, ch3Spline);
	  totalIntegral = integralCh3/temp;
	  temp = topSPESpline(timestamp, splineTimestamp, ch1Spline);
	  totalIntegral += integralCh1/temp;
	}
	BkgHist->Fill(totalIntegral);
      }
    }
  }
  //Normalize by time and subtract bkg from Co
  Double_t timeCo = 0;
  Double_t timeBkg = 0;
  Double_t ratio = 0;

  timeCo = timeEndCo - timeStartCo;
  timeBkg = timeEndBkg - timeStartBkg;
  ratio = timeCo/timeBkg;
  CoHist->Add(BkgHist, -ratio);
  CoHist->GetXaxis()->SetTitle("PE");
  CoHist->GetYaxis()->SetTitle("Subtracted Counts");

  //Set fit range to loop over
  const Int_t stepSize = 5;
  Double_t lowerLower = lowerBound;
  Double_t upperLower = lowerBound;
  Double_t lowerUpper = upperBound;
  Double_t upperUpper = upperBound;
  const Int_t numberOfLower = ((upperLower-lowerLower)/stepSize)+1;
  const Int_t numberOfUpper = ((upperUpper-lowerUpper)/stepSize)+1;
  const Int_t totalNumber = numberOfLower*numberOfUpper;
  TF1** fullFit = new TF1*[totalNumber];
  TF1** expoFit = new TF1*[totalNumber];
  TF1** gausFit = new TF1*[totalNumber];

  //Store certain fit parameters
  Double_t expoPar[totalNumber][3];
  Double_t gausPar[totalNumber][3];
  Double_t fitMean;
  Double_t fitError;
  Double_t tryLower;
  Double_t tryUpper;
  Double_t minLower;
  Double_t minUpper;
  Double_t minMean;
  Double_t minError = 100; //Set an arbitrarily high value
  Double_t minChi2PerNDF;
  Double_t minChi2Desired = 100; //Set an arbitrarily high value to try to search under
  Int_t minCounter;
  counter = 0;

  Double_t chi2 = 0;
  Double_t NDF = 0;
  Double_t chi2PerNDF = 0;
  Double_t chi2Desired = 0;
  Double_t chi2Threshold = 100;
  Double_t minChi2;
  Double_t minNDF;

  //Loop through fit ranges
  for (Int_t i = 0; i < numberOfLower; i++)
  {
    for (Int_t j = 0; j < numberOfUpper; j++)
    {
      tryLower = lowerLower + stepSize*i;
      tryUpper = lowerUpper + stepSize*j;

      //Try fitting unnormalized hist with expo and gaus to get initial values
      expoFit[counter] = new TF1("expoFit", "[0]*exp([1]+[2]*x)", tryLower, tryUpper);
      gausFit[counter] = new TF1("gausFit", "gaus", tryLower, tryUpper);
      UnNormHist->Fit(expoFit[counter], "QR");
      UnNormHist->Fit(gausFit[counter], "QR");
      expoFit[counter]->GetParameters(&expoPar[counter][0]);
      gausFit[counter]->GetParameters(&gausPar[counter][0]);

      //Now fit with full function
      fullFit[counter] = new TF1("fullFit", "[0]*exp([1]+[2]*x)+[3]*exp(-0.5*((x-[4])/[5])**2)", tryLower, tryUpper);
      fullFit[counter]->SetParameter(1, expoPar[counter][1]);
      fullFit[counter]->SetParameter(2, expoPar[counter][2]);
      fullFit[counter]->SetParameter(4, gausPar[counter][1]);
      fullFit[counter]->SetParameter(5, gausPar[counter][2]);
      CoHist->Fit(fullFit[counter], "QR");

      //Extract necessary values
      fitMean = fullFit[counter]->GetParameter(4);
      fitError = fullFit[counter]->GetParError(4);
      chi2 = fullFit[counter]->GetChisquare();
      NDF = fullFit[counter]->GetNDF();
      chi2PerNDF = chi2/NDF;
      chi2Desired = abs(chi2PerNDF-1.0);
      if (chi2Desired <= minChi2Desired && fitMean >= tryLower && fitMean <= tryUpper)
      {
        //Try to get chi2/DoF as close to 1 as possible. Only accept mean values between fit range
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
  printf("For Co run %d and Bkg run %d, lower bound = %f, upper bound = %f, mean = %f, sigma = %f, chi2/DoF = %f\n", coDesiredRun, bkgDesiredRun, minLower, minUpper, minMean, minError, minChi2PerNDF);

  //Draw
  TCanvas *c1 = new TCanvas("c1", "", 1000, 600);
  gStyle->SetOptStat(0);
  CoHist->Draw("hist");
  fullFit[minCounter]->Draw("SAME");
  c1->Update();

  TString fileName = "/nfs_home/bvm41/cenns10Calibration/CoCheckPlots/"+TString::Format("Run%d",coDesiredRun)+".root"; 
  TFile *MyFile = new TFile(fileName, "RECREATE");
  CoHist->Write(TString::Format("CoHistRun%d", coDesiredRun));
  MyFile->Close();
}
