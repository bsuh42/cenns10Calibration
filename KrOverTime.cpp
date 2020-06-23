#include <iostream>
#include <string>
using namespace std;
#include <cmath>
#include <vector>

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

void fitKr(TString KrFile, const Int_t injectionRun)
{
  //Show Kr rate as a function of time. Before injection run is used for background
  TFile *FKr = new TFile(KrFile);
  TTree *TKr = (TTree *)FKr->Get("compressedTree");

  //Variables
  Double_t integralCh1 = 0;
  Double_t integralCh3 = 0;
  Double_t totalIntegral = 0;
  Int_t runNumber = 0;
  Double_t baselineCh1 = 0;
  Double_t baselineCh3 = 0;
  Double_t timestamp =0;
  Double_t temp = 0;
  Int_t counter = 0;
  Int_t evtNumber = 0;
  Int_t index = 0;

  //Read the spline file
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
  vector<TH1D*> histList;
  vector<TH1D*> rawHistList;
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
      printf("Currently on event %d of %d\n", eventNumber, numberEntries);
    }

    //Check to see if new run
    std::vector<int>::iterator it = std::find(KrRunNumbers.begin(), KrRunNumbers.end(), runNumber);
    //printf("ping\n");
    if (integralCh1 != 0 && integralCh3 != 0 && baselineCh1 != 0 && baselineCh3 != 0)
    {
      //printf("ping\n");
      if (it != KrRunNumbers.end())
      {
        index = std::distance(KrRunNumbers.begin(), it);
      }
      else
      {
        //printf("ping2\n");
        KrRunNumbers.push_back(runNumber);
	timeStart.push_back(timestamp);
	timeEnd.push_back(timestamp);
	//Create histograms
	TH1D *hist = new TH1D(TString::Format("Run%d", runNumber), TString::Format("Run%d", runNumber), numberBins, histStart, histEnd);
	TH1D *rawHist = new TH1D("", "", numberBins, histStart, histEnd);
	histList.push_back(hist);
	rawHistList.push_back(rawHist);
	index = 0;
	numberOfEvents.push_back(0);
	//printf("ping2\n");
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
      //printf("ping2\n");
      //Fill histogram
      histList.at(index)->Fill(totalIntegral);
      rawHistList.at(index)->Fill(totalIntegral);
      numberOfEvents[index] = numberOfEvents[index] + 1;
      //printf("ping\n");
    }
  }

  //Normalize by time
  Double_t timeLength = 0;
  Double_t ratio = 0;
  vector<Double_t> timeSince;
  Double_t totalTimeBkg = 0;
  Double_t totalTimeKr = 0;
  Double_t totalNumberBkg = 0;
  Double_t totalNumberKr = 0;
  Int_t krStartIndex = 0;
  totalNumberBkg = 0;

  //Combined histograms
  TH1D* combinedBkgHist = new TH1D("bkg", "bkg", numberBins, histStart, histEnd);
  TH1D* combinedKrHist = new TH1D("Kr", "Kr", numberBins, histStart, histEnd);
  TH1D* combinedRawHist = new TH1D("", "", numberBins, histStart, histEnd);
  for (Int_t runIndex = 0; runIndex < KrRunNumbers.size(); runIndex++)
  {
    timeLength = (timeEnd[runIndex] - timeStart[runIndex]);
    if (KrRunNumbers[runIndex] == injectionRun)
    {
      krStartIndex = runIndex;
    }
    if (KrRunNumbers[runIndex] < injectionRun)
    {
      totalNumberBkg += numberOfEvents[runIndex];
      totalTimeBkg += timeLength;
      combinedBkgHist->Add(histList.at(runIndex));
    }
    if (KrRunNumbers[runIndex] >= injectionRun)
    {
      totalNumberKr += numberOfEvents[runIndex];
      totalTimeKr += timeLength;
      combinedKrHist->Add(histList.at(runIndex));
      combinedRawHist->Add(rawHistList.at(runIndex));

      ratio = timeLength/(totalTimeBkg*1.0);
      histList.at(runIndex)->Add(combinedBkgHist, -ratio);
    }
    timeSince.push_back(timeStart[runIndex]-timeStart[0]);
    numberOfEvents[runIndex] = (numberOfEvents[runIndex]*1.0)/(timeLength*1.0);
  }

  ratio = totalTimeKr/(totalTimeBkg*1.0);
  combinedKrHist->Add(combinedBkgHist, -ratio);
  combinedBkgHist->Scale(1.0/totalTimeBkg);
  //combinedKrHist->Scale(1.0/totalTimeKr);
  totalNumberBkg = totalNumberBkg/(totalTimeBkg*1.0);

  for (Int_t runIndex = 0; runIndex < KrRunNumbers.size(); runIndex++)
  {
    numberOfEvents[runIndex] = numberOfEvents[runIndex]-totalNumberBkg;
  }

  Double_t center = 200;
  Double_t lowerLower = center - 110;
  Double_t upperLower = center - 10;
  Double_t lowerUpper = center + 10;
  Double_t upperUpper = center + 110;
  Double_t stepSize = 10;
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

  TF1** gausFit = new TF1*[totalNumber];
  TF1** unNormFit = new TF1*[totalNumber];
  for(Int_t lower = 0; lower < numberOfLower; lower++)
  {
    for (Int_t upper = 0; upper < numberOfUpper; upper++)
    {
      tryLower = lowerLower + stepSize*lower;
      tryUpper = lowerUpper + stepSize*upper;

      unNormFit[counter] = new TF1("", "gaus", tryLower, tryUpper);
      gausFit[counter] = new TF1("gausFit", "gaus", tryLower, tryUpper);
      combinedRawHist->Fit(unNormFit[counter], "QR");

      gausFit[counter]->SetParameter(1, unNormFit[counter]->GetParameter(1));
      gausFit[counter]->SetParameter(2, unNormFit[counter]->GetParameter(2));
      combinedKrHist->Fit(gausFit[counter], "QR");

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

  //Now fit each kr run to see if constant throughout run
  TF1** individualGausFit = new TF1*[KrRunNumbers.size()];
  TF1** individualUnNormFit = new TF1*[KrRunNumbers.size()];

  vector<Double_t> individualMean;
  vector<Double_t> individualError;
  vector<Double_t> timeSinceError;
  for (Int_t runIndex = 0; runIndex < KrRunNumbers.size(); runIndex++)
  {
    individualUnNormFit[runIndex] = new TF1("", "gaus", minLower, minUpper);
    individualGausFit[runIndex] = new TF1("", "gaus", minLower, minUpper);
    rawHistList.at(runIndex)->Fit(individualUnNormFit[runIndex], "QR");

    individualGausFit[runIndex]->SetParameter(1, individualUnNormFit[runIndex]->GetParameter(1));
    individualGausFit[runIndex]->SetParameter(2, individualUnNormFit[runIndex]->GetParameter(2));
    histList.at(runIndex)->Fit(individualGausFit[runIndex], "QR");

    individualMean.push_back(individualGausFit[runIndex]->GetParameter(1));
    individualError.push_back(individualGausFit[runIndex]->GetParError(1));
    timeSinceError.push_back(0);
    numberOfEvents[runIndex] = numberOfEvents[runIndex] - totalNumberBkg;
  }

  TCanvas *c1 = new TCanvas("c1", "", 1000, 600);
  TCanvas *c2 = new TCanvas("c2", "", 1000, 600);
  TCanvas *c3 = new TCanvas("c3", "", 1000, 600);

  c1->cd();
  gStyle->SetOptStat(0);
  c1->SetLogy();
  combinedKrHist->Draw("hist");
  gausFit[minCounter]->Draw("SAME");
  c1->Update();

  c2->cd();
  gStyle->SetOptStat(0);
  TGraph *gr = new TGraph(numberOfEvents.size(), &timeSince[0], &numberOfEvents[0]);

  TF1* expoFit = new TF1("", "[0]*exp([1]+[2]*x)", timeSince[krStartIndex+100], timeSince[timeSince.size()-1]);
  gr->Fit(expoFit, "R");
  gr->Draw("AP");
  expoFit->Draw("SAME");
  c2->Update();

  c3->cd();
  gStyle->SetOptStat(0);
  TGraphErrors *grErrors = new TGraphErrors(KrRunNumbers.size(), &timeSince[0], &individualMean[0], &timeSinceError[0], &individualError[0]);
  grErrors->Draw("AP");
  c2->Update();
}
