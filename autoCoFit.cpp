#include <iostream>
#include <fstream>
#include <string>
using namespace std;
#include <cmath>
#include <vector>

double bottomSPEFit(Double_t timestamp)
{
  //Get SPE value by fitting LED data
  Double_t a = 17.271;
  Double_t b = 1.895e-8;
  Double_t c = 17.271;
  Double_t d = 2.928e-4;
  Double_t value = 0;

  value = a + b*timestamp + c*exp(-d/timestamp);
  return value;
}

double bottomSPESpline(Double_t timestamp, vector<Double_t> splineTimestamp, vector<Double_t> ch3Spline)
{
  //Get SPE value from spline
  Double_t readValue = splineTimestamp[0];
  Double_t value = ch3Spline[0];
  Int_t counter = 0;
  while (timestamp >= readValue)
  {
    //Look for nearest timestamp value
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

void autoCoFit()
{
  //Using qsub, automatically find the Co peak and print output to text file
  //Recommend grouping files in a logical manner for ease of viewing

  //File list
  vector<TString> fileList;
  //fileList.push_back("/data6/coherent/data/LiqAr/benRedoneAnalysis2020/Run01Co08/Compressed_Run01_Co08.root");
  fileList.push_back("/data6/coherent/data/LiqAr/benRedoneAnalysis2020/Run01Co15/Compressed_Run01_Co15.root");
  //fileList.push_back("/data6/coherent/data/LiqAr/benRedoneAnalysis2020/Run01Co22/Compressed_Run01_Co22.root");
  //fileList.push_back("/data6/coherent/data/LiqAr/benRedoneAnalysis2020/Run02Co08/Compressed_Run02_Co08.root");
  fileList.push_back("/data6/coherent/data/LiqAr/benRedoneAnalysis2020/Run02Co15/Compressed_Run02_Co15.root");
  //fileList.push_back("/data6/coherent/data/LiqAr/benRedoneAnalysis2020/Run02Co22/Compressed_Run02_Co22.root");
  //fileList.push_back("/data6/coherent/data/LiqAr/benRedoneAnalysis2020/Run03Co08/CompressedCo.root");
  fileList.push_back("/data6/coherent/data/LiqAr/benRedoneAnalysis2020/Run03Co15/CompressedCo.root");
  //fileList.push_back("/data6/coherent/data/LiqAr/benRedoneAnalysis2020/Run03Co22/CompressedCo.root");
  //fileList.push_back("/data6/coherent/data/LiqAr/benRedoneAnalysis2020/Run04Co08/CompressedCo.root");
  fileList.push_back("/data6/coherent/data/LiqAr/benRedoneAnalysis2020/Run04Co15/CompressedCo.root");
  //fileList.push_back("/data6/coherent/data/LiqAr/benRedoneAnalysis2020/Run04Co22/CompressedCo.root");
  //fileList.push_back("/data6/coherent/data/LiqAr/benRedoneAnalysis2020/Run05Co08/CompressedCo.root");
  fileList.push_back("/data6/coherent/data/LiqAr/benRedoneAnalysis2020/Run05Co15/CompressedCo.root");
  //fileList.push_back("/data6/coherent/data/LiqAr/benRedoneAnalysis2020/Run05Co22/CompressedCo.root");
  
  //Background file list
  vector<TString> bkgFileList;
  bkgFileList.push_back("/data6/coherent/data/LiqAr/benRedoneAnalysis2020/Run01Co15/Compressed_Run01_Bkg.root");
  bkgFileList.push_back("/data6/coherent/data/LiqAr/benRedoneAnalysis2020/Run02Co15/Compressed_Run02_Bkg.root");
  bkgFileList.push_back("/data6/coherent/data/LiqAr/benRedoneAnalysis2020/Run03Co15/CompressedBkg.root");
  bkgFileList.push_back("/data6/coherent/data/LiqAr/benRedoneAnalysis2020/Run04Co15/CompressedBkg.root");
  bkgFileList.push_back("/data6/coherent/data/LiqAr/benRedoneAnalysis2020/Run05Co15/CompressedBkg.root");

  //Output file
  ofstream outputFile;
  outputFile.open("/nfs_home/bvm41/cenns10Calibration/avg15.txt");

  //Open input files and put into a TChain
  TChain chCo("compressedTree");
  TChain chBkg("compressedTree");
  for (Int_t i = 0; i < fileList.size(); i++)
  {
    chCo.Add(fileList[i]);
  }
  for (Int_t i = 0; i < bkgFileList.size(); i++)
  {
    chBkg.Add(bkgFileList[i]);
  }

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
  const Int_t histStart = 400;
  const Int_t histEnd = 900; //Adjust as necessary
  const Int_t numberBins = 100;
  vector<TH1D*> CoHistList;
  vector<TH1D*> BkgHistList;
  vector<TH1D*> RawHistList;
  vector<Double_t> timeStartCo;
  vector<Double_t> timeEndCo;
  vector<Double_t> timeStartBkg;
  vector<Double_t> timeEndBkg;
  vector<Int_t> CoRunNumbers;
  vector<Int_t> BkgRunNumbers;
  
  //Variables
  Double_t integralCh1 = 0;
  Double_t integralCh3 = 0;
  Double_t totalIntegral = 0;
  Int_t runNumber = 0;
  Double_t temp = 0;
  Int_t counter = 0;
  Int_t evtNumber = 0;
  Double_t baselineCh1 = 0;
  Double_t baselineCh3 = 0;
  Double_t timestamp = 0;
  Int_t index = 0;

  //For fixed value
  const Double_t SPE3 = 64.0;
  const Double_t SPE1Factor = 1.18; //For fixed and fit, multiply SPE by this factor to get ch1 SPE
  const TString SPEMethod = "spline"; //choose method for determining SPE. Can choose spline, fit, or fixed

  //Set branches to get
  chCo.SetBranchAddress("runNumber", &runNumber);
  chCo.SetBranchAddress("integralCh1", &integralCh1);
  chCo.SetBranchAddress("integralCh3", &integralCh3);
  chCo.SetBranchAddress("baselineCh1", &baselineCh1);
  chCo.SetBranchAddress("baselineCh3", &baselineCh3);
  chCo.SetBranchAddress("timestamp", &timestamp);
  
  chBkg.SetBranchAddress("runNumber", &runNumber);
  chBkg.SetBranchAddress("integralCh1", &integralCh1);
  chBkg.SetBranchAddress("integralCh3", &integralCh3);
  chBkg.SetBranchAddress("baselineCh1", &baselineCh1);
  chBkg.SetBranchAddress("baselineCh3", &baselineCh3);
  chBkg.SetBranchAddress("timestamp", &timestamp);

  Int_t CoNumberEntries = chCo.GetEntries();
  for (Int_t eventNumber = 0; eventNumber < CoNumberEntries; eventNumber++)
  {
    //Co run
    totalIntegral = 0;
    index = 0;
    chCo.GetEntry(eventNumber);
    std::vector<int>::iterator it = std::find(CoRunNumbers.begin(), CoRunNumbers.end(), runNumber);
    if (integralCh1 != 0 && integralCh3 != 0 && baselineCh1 != 0 && baselineCh3 != 0)
    {
      if (it != CoRunNumbers.end())
      {
        //check to see if this run is already in the vector
	index = std::distance(CoRunNumbers.begin(), it);
      }
      else
      {
        CoRunNumbers.push_back(runNumber);
	timeStartCo.push_back(timestamp);
	timeEndCo.push_back(timestamp);
	//Create histograms for subtractured histogram and unnormalized histogram
	TH1D *hist = new TH1D(TString::Format("Co Run%d", runNumber), TString::Format("Co Run%d", runNumber), numberBins, histStart, histEnd);
	TH1D *rawHist = new TH1D("", "", numberBins, histStart, histEnd);
	CoHistList.push_back(hist);
	RawHistList.push_back(rawHist);
	index = 0;
      }
      //timing
      if (timestamp < timeStartCo[index])
      {
        timeStartCo[index] = timestamp;
      }
      if (timestamp > timeEndCo[index])
      {
        timeEndCo[index] = timestamp;
      }

      if (SPEMethod == "spline")
      {
        temp = topSPESpline(timestamp, splineTimestamp, ch1Spline);
	totalIntegral = integralCh1/(temp*1.0);
	temp = bottomSPESpline(timestamp, splineTimestamp, ch3Spline);
	totalIntegral += integralCh3/(temp*1.0);
      }
      if (SPEMethod == "fit")
      {
        temp = bottomSPEFit(timestamp);
	totalIntegral = integralCh1/temp;
	totalIntegral += integralCh3/(temp*SPE1Factor);
      }
      if (SPEMethod == "fixed")
      {
        totalIntegral = integralCh1/(SPE3*SPE1Factor);
	totalIntegral += integralCh3/SPE3;
      }
      CoHistList.at(index)->Fill(totalIntegral);
      RawHistList.at(index)->Fill(totalIntegral);
    }
  }

  Int_t BkgNumberEntries = chBkg.GetEntries();
  for (Int_t eventNumber = 0; eventNumber < BkgNumberEntries; eventNumber++)
  {
    //Perform the same exercise for the background
    totalIntegral = 0;
    index = 0;
    chBkg.GetEntry(eventNumber);
    std::vector<int>::iterator it = std::find(BkgRunNumbers.begin(), BkgRunNumbers.end(), runNumber);
    if (integralCh1 != 0 && integralCh3 != 0 && baselineCh1 != 0 && baselineCh3 != 0)
    {
      if (it != BkgRunNumbers.end())
      {
        index = std::distance(BkgRunNumbers.begin(), it);
      }
      else
      {
        BkgRunNumbers.push_back(runNumber);
	timeStartBkg.push_back(timestamp);
	timeEndBkg.push_back(timestamp);
	TH1D *bkgHist = new TH1D(TString::Format("Bkg Run %d", runNumber), TString::Format("Bkg Run %d", runNumber), numberBins, histStart, histEnd);
	BkgHistList.push_back(bkgHist);
	index = 0;
      }
      if (timestamp < timeStartBkg[index])
      {
        timeStartBkg[index] = timestamp;
      }
      if (timestamp > timeEndBkg[index])
      {
        timeEndBkg[index] = timestamp;
      }

      if (SPEMethod == "spline")
      {
        temp = topSPESpline(timestamp, splineTimestamp, ch1Spline);
	totalIntegral = integralCh1/(temp*1.0);
	temp = bottomSPESpline(timestamp, splineTimestamp, ch3Spline);
	totalIntegral += integralCh3/(temp*1.0);
      }
      if (SPEMethod == "fit")
      {
        temp = bottomSPEFit(timestamp);
	totalIntegral = integralCh1/temp;
	totalIntegral += integralCh3/(temp*SPE1Factor);
      }
      if (SPEMethod == "fixed")
      {
        totalIntegral = integralCh1/(SPE3*SPE1Factor);
	totalIntegral += integralCh3/SPE3;
      }
      BkgHistList.at(index)->Fill(totalIntegral);
    }
  }

  //For each Co run, normalize, bkg subtract, then fit
  //Timing variables
  Double_t timeCo = 0;
  Double_t timeBkg = 0;
  Double_t ratio = 0;
  Int_t binmax;
  Double_t center;
  Bool_t foundBkg = 0;
  Int_t currentBkgIndex = 0;

  //Fit variables
  Double_t lowerLower;
  Double_t upperLower;
  Double_t lowerUpper;
  Double_t upperUpper;
  Int_t numberOfLower;
  Int_t numberOfUpper;
  Int_t totalNumber;
  Double_t stepSize = 5; //Adjust as needed
  
  Double_t expo1;
  Double_t expo2;
  Double_t gaus1;
  Double_t gaus2;
  Double_t fitMean;
  Double_t fitError;
  Double_t minLower;
  Double_t minUpper;
  Double_t minMean;
  Double_t minError;
  Double_t minChi2PerNDF;
  Double_t minChi2Desired;
  Int_t minCounter;
  Double_t chi2 = 0;
  Double_t NDF = 0;
  Double_t chi2PerNDF = 0;
  Double_t chi2Desired = 0;
  Double_t tryLower = 0;
  Double_t tryUpper = 0;
  Double_t minChi2;
  Double_t minNDF;
  Double_t minExpo1;
  Double_t minExpo2;
  Double_t minGaus1;
  Double_t minGaus2;

  //Save hist to an output file
  TFile *MyFile;
  TString fileName;

  //For each Co Run, find the nearest Bkg run
  for (Int_t CoIndex = 0; CoIndex < CoRunNumbers.size(); CoIndex++)
  {
    foundBkg = 0;
    timeCo = 0;
    timeBkg = 0;
    ratio = 0;
    currentBkgIndex = 0;
    for (Int_t BkgIndex = 0; BkgIndex < BkgRunNumbers.size(); BkgIndex++)
    {
      if (BkgRunNumbers[BkgIndex] >= CoRunNumbers[CoIndex]-20 && BkgRunNumbers[BkgIndex] <= CoRunNumbers[CoIndex]+20)
      {
        foundBkg = 1;
	currentBkgIndex = BkgIndex;
      }
    }

    if (foundBkg)
    {
      timeCo = timeEndCo[CoIndex] - timeStartCo[CoIndex];
      timeBkg = timeEndBkg[currentBkgIndex]-timeStartBkg[currentBkgIndex];
      ratio = timeCo/(timeBkg*1.0);
      CoHistList.at(CoIndex)->Add(BkgHistList.at(currentBkgIndex), -ratio);

      //Find max points
      binmax = CoHistList.at(CoIndex)->GetMaximumBin();
      center = CoHistList.at(CoIndex)->GetXaxis()->GetBinCenter(binmax);
      
      lowerLower = center - 170;
      upperLower = center - 10;
      lowerUpper = center + 10;
      upperUpper = center + 170;
      numberOfLower = ((upperLower-lowerLower)/stepSize)+1;
      numberOfUpper = ((upperUpper-lowerUpper)/stepSize)+1;

      //Fit with expo+gaus
      totalNumber = numberOfLower*numberOfUpper;
      TF1** fullFit = new TF1*[totalNumber];
      TF1** expoFit = new TF1*[totalNumber];
      TF1** gausFit = new TF1*[totalNumber];
      
      minError = 100;
      minChi2Desired = 100; //Set arbitarily high to normalize
      counter = 0;

      for (Int_t lower = 0; lower < numberOfLower; lower++)
      {
        for (Int_t upper = 0; upper < numberOfUpper; upper++)
	{
	  tryLower = lowerLower + stepSize*lower;
	  tryUpper = lowerUpper + stepSize*upper;
	  
	  expoFit[counter] = new TF1("expoFit", "[0]*exp([1]+[2]*x)", tryLower, tryUpper);
	  gausFit[counter] = new TF1("gausFit", "gaus", tryLower, tryUpper);
	  RawHistList.at(CoIndex)->Fit(expoFit[counter], "QR");
	  RawHistList.at(CoIndex)->Fit(gausFit[counter], "QR");
	  
	  expo1 = expoFit[counter]->GetParameter(1);
	  expo2 = expoFit[counter]->GetParameter(2);
	  gaus1 = gausFit[counter]->GetParameter(1);
	  gaus2 = gausFit[counter]->GetParameter(2); //fit with expo and gaus first to get initial parameters
	  
	  fullFit[counter] = new TF1("fullFit", "[0]*exp([1]+[2]*x)+[3]*exp(-0.5*((x-[4])/[5])**2)", tryLower, tryUpper);
	  fullFit[counter]->SetParameter(1, expo1);
	  fullFit[counter]->SetParameter(2, expo2);
	  fullFit[counter]->SetParameter(4, gaus1);
	  fullFit[counter]->SetParameter(5, gaus2);
	  CoHistList.at(CoIndex)->Fit(fullFit[counter], "QR");
	  
	  fitMean = fullFit[counter]->GetParameter(4);
	  fitError = fullFit[counter]->GetParError(4);
	  chi2 = fullFit[counter]->GetChisquare();
	  NDF = fullFit[counter]->GetNDF();
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
	    minExpo1 = expo1;
	    minExpo2 = expo2;
	    minGaus1 = gaus1;
	    minGaus2 = gaus2;
	  }
	  counter += 1;
	}
      }
      outputFile << BkgRunNumbers[currentBkgIndex] << " " << CoRunNumbers[CoIndex] << " " << minLower << " " << minUpper << " " << minMean << " " << minError << " " << minChi2PerNDF << "\n";
      //CoHistList.at(CoIndex)->Fit(fullFit[minCounter], "QR");
      //
      //Fit again to save hist
      TF1 *finalFit = new TF1("", "[0]*exp([1]+[2]*x)+[3]*exp(-0.5*((x-[4])/[5])**2)", minLower, minUpper);
      finalFit->SetParameter(1, minExpo1);
      finalFit->SetParameter(2, minExpo2);
      finalFit->SetParameter(4, minGaus1);
      finalFit->SetParameter(5, minGaus2);
      CoHistList.at(CoIndex)->Fit(finalFit, "QR");

      fileName = "/nfs_home/bvm41/cenns10Calibration/CoCheckPlots/"+TString::Format("Run%d",CoRunNumbers[CoIndex])+".root";
      MyFile = new TFile(fileName, "RECREATE");
      CoHistList.at(CoIndex)->Write("hist");
      //fullFit[minCounter]->Write("function");
      MyFile->Close();
    }
  }
}
