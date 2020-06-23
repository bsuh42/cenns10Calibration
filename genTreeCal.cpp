#include <iostream>
#include <string>
using namespace std;
#include <cmath>
#include <vector>
#include "TROOT.h"
#include "TFile.h"
#include "TTree.h"
#include "TH1.h"
#include "TH2.h"
#include "TF1.h"
#include "TLegend.h"
#include "TSystem.h"
#include "TCanvas.h"
#include "TChain.h"
#include <map>
#include <limits>
#include <cstdio>
#include <cstdlib>
#include <getopt.h>
#include <unistd.h>
#include <typeinfo>

#include "EventData.hh"
#include "ChannelData.hh"
#include "Pulse.hh"
#include "ScinEvt.hh"

void makeTree(char *infile, char *outputname)
{
  //Compressed for cobalt calibrations
  //Load using .L genTreeCal.cpp+
  //Grabs scintillation events

  TChain ch("Events");
  ch.Add(infile);
  EventData *event = 0;
  ch.SetBranchAddress("event", &event);

  Int_t numEntries = ch.GetEntries();
  printf("There are %d events in the file\n", numEntries);

  //Variables
  Double_t temp = 0;
  Int_t counter = 0;
  Int_t channelID = 0;
  Int_t runNumber = 0;
  Double_t timestamp = 0;
  Double_t eventTime = 0;
  Double_t baselineCh1 = 0;
  Double_t baselineCh3 = 0;
  Double_t integralCh1 = 0;
  Double_t integralCh3 = 0;
  Int_t evtNumber = 0;
  Bool_t writeEvent = 0;
  const Double_t polarity = -1.0;

  //Create target file and tree to store branches
  TFile *outfile = TFile::Open(outputname, "RECREATE");
  TTree *outtree = new TTree("compressedTree", "stores data");
  outtree->Branch("evtNumber", &evtNumber);
  outtree->Branch("timestamp", &timestamp);
  outtree->Branch("eventTime", &eventTime);
  outtree->Branch("runNumber", &runNumber);
  outtree->Branch("baselineCh1", &baselineCh1);
  outtree->Branch("baselineCh3", &baselineCh3);
  outtree->Branch("integralCh1", &integralCh1);
  outtree->Branch("integralCh3", &integralCh3);

  for (Int_t eventNumber = 0; eventNumber <numEntries; eventNumber++)
  {
    integralCh1 = 0;
    integralCh3 = 0;
    writeEvent = 0;
    if (eventNumber % 1000 == 0)
    {
      printf("Currently on event %d of %d\n", eventNumber, numEntries);
    }
    ch.GetEntry(eventNumber);
    evtNumber = eventNumber;
    eventTime = event->event_time;
    timestamp = event->timestamp;
    runNumber = event->run_id;

    //Loop through channels
    std::vector<ChannelData> chanVector;
    chanVector = event->channels;
    std::vector<ChannelData>::iterator chanBeg = chanVector.begin();
    std::vector<ChannelData>::iterator chanEnd = chanVector.end();
    std::vector<ChannelData>::iterator chanIt = chanVector.begin();

    for (chanIt = chanBeg; chanIt != chanEnd; ++chanIt)
    {
      channelID = chanIt->channel_id;
      //Loop through scintillation events
      std::vector<ScinEvt> ScinEvts = chanIt->vScinEvts;
      std::vector<ScinEvt>::iterator ScinEvtBeg = ScinEvts.begin();
      std::vector<ScinEvt>::iterator ScinEvtEnd = ScinEvts.end();
      std::vector<ScinEvt>::iterator ScinEvtIt = ScinEvts.begin();

      for (; ScinEvtIt != ScinEvtEnd; ScinEvtIt++)
      {
        if (channelID == 1)
	{
	  integralCh1 = ScinEvtIt->GetFullIntegral();
	  baselineCh1 = chanIt->baseline.mean;
	  writeEvent = 1;
	}
	if (channelID == 3)
	{
	  integralCh3 = ScinEvtIt->GetFullIntegral();
	  baselineCh3 = chanIt->baseline.mean;
	  writeEvent = 1;
	}
      }
    }
    if (writeEvent)
    {
      outtree->Fill();
    }
  }
  outfile->Write();
}
