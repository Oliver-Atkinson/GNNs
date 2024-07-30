#include <TApplication.h>
#include <TGClient.h>

#include <fstream>
#include "TLegend.h"
#include "TCanvas.h"
#include "TPad.h"
#include "TLatex.h"

#include "SimpleAnalysis.h"
#include "SimpleAnalysisHEPMC.h"
#include "data1.h"


int main(int argc, char **argv)
{ 
  
  double sevents2[10],sevents[10],seventscuts[10],seventscuts2[10];
  double eventscuts[10],events[10];
  TApplication theApp("App",&argc,argv);
  double max[10],min[10];
  
  for (unsigned ii=0; ii<10; ii++){
    max[ii]=0;
    min[ii]=0;
    sevents[ii]=0;
    sevents2[ii]=0;
    events[ii]=0;   
    
    seventscuts[ii]=0;
    seventscuts2[ii]=0;
    eventscuts[ii]=0;   
    
  }

  // Define Data Files  
  SimpleAnalysisHEPMC <data1>signal;     
  SimpleAnalysisHEPMC <data1>signal2;     
  SimpleAnalysisHEPMC <data1>signal3;     
  SimpleAnalysisHEPMC <data1>signal4;     
  SimpleAnalysisHEPMC <data1>signal5;     
  SimpleAnalysisHEPMC <data1>signal6;     
  SimpleAnalysisHEPMC <data1>signal7;     
  SimpleAnalysisHEPMC <data1>signal8;     
  SimpleAnalysisHEPMC <data1>signal9;     
  SimpleAnalysisHEPMC <data1>bkg;     
  
  // Define Colors
  //
  signal.SetColor(1);
  bkg.SetColor(2);
  signal2.SetColor(3);
  signal3.SetColor(4);
  signal4.SetColor(4);
  signal5.SetColor(4);
  signal6.SetColor(4);
  signal7.SetColor(4);
  signal8.SetColor(4);
  signal9.SetColor(4);
  
  // Analyze
  //
  cout << "#--------------------------------------------------------------------------" << endl;
  cout << "#---   SimpleAnalysis                                                   ---" << endl;
  cout << "#---                                                                    ---" << endl;
  cout << "#---                                                         enjoy....  ---" << endl;
  cout << "#--------------------------------------------------------------------------" << endl;   
  
  double weight=1.0;  


  signal.Analyze("topjets.hepmc",sevents[0],seventscuts[0],true,weight,2);
  //  signal.Analyze("Tops_100k_550GeV/tag_1_pythia8_events_1.hepmc",sevents[0],seventscuts[0],false,weight);

  bkg.Analyze("qcdjets.hepmc",events[0],eventscuts[0],true,weight,1);
  //  bkg.Analyze("QCD_200k_550GeV/tag_1_pythia8_events_1.hepmc",events[0],eventscuts[0],false,weight);
  //  bkg.Analyze("QCD_200k_550GeV/tag_1_pythia8_events_2.hepmc",events[0],eventscuts[0],false,weight);
  //  bkg.Analyze("QCD_200k_550GeV/tag_1_pythia8_events_3.hepmc",events[0],eventscuts[0],false,weight);


  signal2.Analyze("wjets.hepmc",sevents[1],seventscuts[1],true,weight,3);
   //  signal2.Analyze("WZ_100k_550GeV/tag_1_pythia8_events_1.hepmc",sevents[1],seventscuts[1],false,weight);

  // Loop over Observables
  //
  map <string, MyChart*> :: iterator it;
  gStyle-> SetPalette(1,0);
  TCanvas *can2 = new TCanvas("plot2", "");
  string filename;
  for (it = signal.charts.begin(); it != signal.charts.end(); ++it )
    {
    
      cout << "Observable: " << it->first << endl;
      can2->Update();
      
      
      // Draw Observables
      string str=it->first;
      cout << str << endl;
      if ( str.find("2d") !=string::npos)
	{
	  signal.Draw(it->first, "COLZ");  //Cp even
	  signal2.Draw(it->first, "same");  //Cp even
	  signal3.Draw(it->first, "same");  //Cp even
	  bkg.Draw(it->first, "same");  //Cp even


	}
      else
	{

	  //Get maxima for nice-looking plots
	  max[0]=signal.GetMax(it->first);    
	  min[0]=signal.GetMin(it->first);    
	  max[1]=signal2.GetMax(it->first);    
	  min[1]=signal2.GetMin(it->first);    
	  max[2]=bkg.GetMax(it->first);    
	  min[2]=bkg.GetMin(it->first);    
	  int cnt;
	  double histmax=0;
	  double histmin=0;
	  for (int k=0;k<=3; k++) 
	    {
	      if (max[k]>histmax) 
		{
		  histmax=max[k];
		  cnt=k;
		}
	      if (min[k]<histmin) 
		{
		  histmin=min[k];
		}
	    }


	  signal.Range(it->first,histmin*1.1, histmax*1.1);
	  signal.Draw(it->first, "");
	  signal2.Draw(it->first, "same"); 
	  bkg.Draw(it->first, "same"); 
	  
	  // Define Legend and output observable_leg.eps
	  //
	  TLegend *leg = new TLegend(0.7, 0.7, 0.98, 0.99);
	  leg->AddEntry(dynamic_cast<TH1F *>(bkg.charts[it->first]->histo), "QCD", "L");
	  leg->AddEntry(dynamic_cast<TH1F *>(signal.charts[it->first]->histo), "T", "same");
	  leg->AddEntry(dynamic_cast<TH1F *>(signal2.charts[it->first]->histo), "W", "same");
	  leg->SetFillColor(10);
	  leg->Draw();	  



	}



    
    
      //Save as pdf
      filename = "./"+it->first+".pdf";
      can2->SaveAs( filename.c_str() );
      
    }
    

  //  theApp.Run();  
  return 0;
}

