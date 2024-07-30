//ROOT STUFF
//
#include "BaseAnalysis.h"
#include "TRandom.h"
#include "TRandom2.h"
#include "math.h"
#include "lester_mt2_bisect.h"

//FASTJET STUFF
#include "fastjet/ClusterSequence.hh"

//HEPMC STUFF
#include "HepMC/IO_GenEvent.h"
#include "HepMC/GenEvent.h"

#include "Angles.h"


//Load namespaces
using namespace fastjet; 
using namespace std;
using namespace HepMC;
using namespace angles;

template<class T>
class SimpleAnalysisHEPMC : public BaseAnalysis<T>{

 public:

    SimpleAnalysisHEPMC() : BaseAnalysis<T>() {
    Init();
    }
    void Init();
    void Analyze(const string hepmcfile, double & eventnumber, double & eventpasscuts, const bool scale, double &  weight, const int procid);
    void DefinePreJets(HepMC::GenEvent* evt, vector<fastjet::PseudoJet> & PartOfPrejets, vector<fastjet::PseudoJet> & Leptons,vector<fastjet::PseudoJet> & Muons, 
		       vector<fastjet::PseudoJet> & Neutrinos);
    vector<fastjet::PseudoJet> DefineJets(const vector<fastjet::PseudoJet> vec1, const double ptmin, const double ymax);
    double Rsep(const fastjet::PseudoJet vec1, const fastjet::PseudoJet vec2);
    double Phisep(const fastjet::PseudoJet vec1, const fastjet::PseudoJet vec2);
};




template <class T>
void SimpleAnalysisHEPMC<T>::Init(){


  //1d histogram
  this->Add("Ptj",50, 0.0, 1000.0, "p_Tjet 1 [GeV]", "d#sigma/dp_T [fb GeV^{-1}]");
  this->Add("mj",50, 0.0, 1000.0, "p_Tjet 1 [GeV]", "d#sigma/dp_T [fb GeV^{-1}]");
  this->Add("mj2",50, 0.0, 1000.0, "p_Tjet 1 [GeV]", "d#sigma/dp_T [fb GeV^{-1}]");
  this->Add("Ptt",50, 0.0, 1000.0, "p_Tjet 1 [GeV]", "d#sigma/dp_T [fb GeV^{-1}]");
  this->Add("Rtj",20, 0.0, 5.0, "p_Tjet 1 [GeV]", "d#sigma/dp_T [fb GeV^{-1}]");
  this->Add("ntt",20, 0.0, 200.0, "p_Tjet 1 [GeV]", "d#sigma/dp_T [fb GeV^{-1}]");
  this->Add("nj",20, 0.0, 200.0, "p_Tjet 1 [GeV]", "d#sigma/dp_T [fb GeV^{-1}]");

}


template <class T>
void SimpleAnalysisHEPMC<T>::Analyze(const string hepmcfile, double & eventnumber, double & eventpasscuts, const bool scale, double &  weight, const int procid ) 
{

  double lweight=weight;
  string filename;
  filename="./"+hepmcfile;
  Long64_t icount=0;
  Long64_t icountcuts=0; 
  Long64_t checkcuts[10]; 
  HepMC::IO_GenEvent ascii_in(filename.c_str(),std::ios::in);
  HepMC::GenEvent* evt = ascii_in.read_next_event();

  for (unsigned kk=0; kk<10; kk++)
    {
      checkcuts[kk]=0;
    }

// Loop over Events in the HEPMC File
  while ( evt )// && icount < 10000 ) 
    {

      icount++;
     


      vector<fastjet::PseudoJet> PartOfPrejets,recjets,recjetsusort,leptons,muons,neutrinos,jets,bjets,mcbjets,nicejets;
      fastjet::PseudoJet missing;
      fastjet::Strategy strategy = Best;
      
      ///////////////////////////////////////////////
      // JETS & LEPTONS & PHOTONS
      // Fix Jetdefinition and Cone size
      ///////////////////////////////////////////////
      double Rparam =1.0;
      JetDefinition jet_def(cambridge_algorithm, Rparam, strategy);    
      ///////////////////////////////////////////////
      // Prejets, Clustering &  Reconstruction
      ///////////////////////////////////////////////
      
      DefinePreJets(evt, PartOfPrejets, leptons,muons,neutrinos);
      fastjet::ClusterSequence cs(PartOfPrejets, jet_def);

      for (unsigned kk=0; kk<muons.size(); kk++)
	{
	  leptons.push_back(muons[kk]);
	}
      muons.clear();
      
      recjetsusort = cs.inclusive_jets();
      recjets = fastjet::sorted_by_pt(recjetsusort);
      
      //select four vectors with pt>50 and |rapidity|<4.5
      jets=DefineJets(recjets,550,4.5);      
      leptons=DefineJets(leptons,12.,2.5);      


          

      ///////////////////////////////////////////////
      //CUTS
      ///////////////////////////////////////////////
      bool cuts=jets.size()>0 ; 

      
      if (cuts)	
	{
	  int doonce=1;
	  for ( HepMC::GenEvent::particle_iterator p = evt->particles_begin();
                p != evt->particles_end(); ++p )
	    {
	      
	      if ( (procid==2 && fabs((*p)->pdg_id())==6) || (procid==3 && fabs((*p)->pdg_id())==23) )
		{	

		  //		  if (procid==2) cout << "fount top" << endl;
				  //		  doonce==2;
		  double maxperp=140000;
		  for ( HepMC::GenVertex::particle_iterator des=(*p)->end_vertex()->particles_begin(HepMC::descendants);
			des != (*p)->end_vertex()->particles_end(HepMC::descendants); ++des )
		    {
		      //		      cout << (*des)->pdg_id() << " " << (*des)->status() << endl;
		      if ( (*des)->status() == 1 )
			{
			  fastjet::PseudoJet temp=PseudoJet( (*des)->momentum().px(),
							     (*des)->momentum().py(),
							     (*des)->momentum().pz(),
							     (*des)->momentum().e());
			  bool dopush=true;
			  for (unsigned qq=0; qq<nicejets.size(); qq++)
			    {
			      if (temp.e()==nicejets[qq].e())
				{
				  dopush=false;
				  break;
				}
			    }
			  if (dopush)
			    {
			      nicejets.push_back(temp);
			    }
			}
		    }
		}		
	      else if  (procid==1)
		{
		  nicejets=jets[0].constituents();
		}
	    }    

	  /*	  
	  for (unsigned kk=0; kk<jets.size(); kk++)
	    {
	      vector<fastjet::PseudoJet> constituents=jets[kk].constituents();
	      cout << "consitutents of jet " << kk << " of event " << icount << endl;
	      for (unsigned qq=0; qq< constituents.size(); qq++)
		{
		  cout << qq+1 <<" "<< constituents[qq].px() << " " << constituents[qq].py() << " " << constituents[qq].pz() << " " << constituents[qq].e() << " " << constituents[qq].user_index() << endl;
		}
	    }
	  */
	  this->Fill("Ptj",jets[0].perp(),lweight);
	  this->Fill("mj",jets[0].m(),lweight);
	  this->Fill("ntt",jets[0].constituents().size(),lweight);
	  this->Fill("nj",jets.size(),lweight);

	  //	  cout << nicejets.size() << endl;

	  if (nicejets.size()>0)
	    {
	      JetDefinition jet_def2(cambridge_algorithm, Rparam, strategy);    
	      fastjet::ClusterSequence cs2(nicejets, jet_def2);      
	      recjetsusort = cs2.inclusive_jets();
	      recjets = fastjet::sorted_by_pt(recjetsusort);
	      
	      this->Fill("mj2",recjets[0].m(),lweight);
	      this->Fill("Ptt",recjets[0].perp(),lweight);
	      this->Fill("Rtj",Rsep(recjets[0],jets[0]),lweight);
	    }
	  
	  icountcuts++;
	      
	      
	}
      
      PartOfPrejets.clear();
      jets.clear();
      leptons.clear();
      bjets.clear();

	
      if (icount==1){ 
	cout << "--------------------------------------------------------------------------------" << endl;
	cout << "**** Analyze: Loop over HEPMC File "<< filename  << endl;
	}

      if ( icount%1000==1 ) cout << "Processing Event Number " << icount << endl; 		  	    
      
      
      //
      // Delete event and get a new one
      //
      
      delete evt;
      ascii_in >> evt;
      
    }

  //////////////////////////////////////
  //Enter Scaling for normalized plots
  //////////////////////////////////////
  if (scale) 
    {
      this->Scale("Ptj", 1. / dynamic_cast<TH1F *>(this->charts["Ptj"]->histo)->Integral() );
      this->Scale("Ptt", 1. / dynamic_cast<TH1F *>(this->charts["Ptt"]->histo)->Integral() );
      this->Scale("Rtj", 1. / dynamic_cast<TH1F *>(this->charts["Rtj"]->histo)->Integral() );
      this->Scale("ntt", 1. / dynamic_cast<TH1F *>(this->charts["ntt"]->histo)->Integral() );
      this->Scale("nj", 1. / dynamic_cast<TH1F *>(this->charts["nj"]->histo)->Integral() );
      this->Scale("mj", 1. / dynamic_cast<TH1F *>(this->charts["mj"]->histo)->Integral() );
      this->Scale("mj2", 1. / dynamic_cast<TH1F *>(this->charts["mj2"]->histo)->Integral() );


    }
  
  eventnumber=eventnumber+double(icount);
  eventpasscuts=eventpasscuts+double(icountcuts);
  
  cout << "#--------------------------------------------------------------------------" << endl;
  cout << eventnumber << " Events in sample" << endl;
  cout << eventpasscuts << " Events pass cuts, that is "<< eventpasscuts/eventnumber  << endl;
  cout << "#--------------------------------------------------------------------------" << endl;
  
  
  
}


template <class T>
vector<fastjet::PseudoJet> SimpleAnalysisHEPMC<T>::DefineJets(const vector<fastjet::PseudoJet> vec1, const double ptmin, const double ymax)
{
  vector<fastjet::PseudoJet> particlesloc;
  particlesloc.clear();
  for (unsigned i=0 ; i<vec1.size(); i++){
    if (vec1[i].perp()>ptmin && fabs(vec1[i].rapidity()) <= ymax){
      particlesloc.push_back(vec1[i]);
    }
  }
  return particlesloc;
}


template <class T>
void SimpleAnalysisHEPMC<T>::DefinePreJets(HepMC::GenEvent* evt,vector<fastjet::PseudoJet> & PartOfPrejets2, vector<fastjet::PseudoJet> & Leptons, vector<fastjet::PseudoJet> & Muons,vector<fastjet::PseudoJet> & Neutrinos)
{
  vector<fastjet::PseudoJet> PartOfPrejets;
  vector<fastjet::PseudoJet> LeptonsLoc;
  vector<fastjet::PseudoJet> MuonsLoc;
  vector<fastjet::PseudoJet> NeutrinosLoc;
  fastjet::PseudoJet temp;

  for ( HepMC::GenEvent::particle_const_iterator p= evt->particles_begin();
	p != evt->particles_end(); ++p ){


    if ( ((*p)->status())==1) {


      if ( fabs((*p)->pdg_id()) < 111 && (*p)->pdg_id()!=22 && fabs((*p)->pdg_id()) != 11 && fabs((*p)->pdg_id()) != 13 
	   && fabs((*p)->pdg_id()) != 12 && fabs((*p)->pdg_id()) != 14 && fabs((*p)->pdg_id()) != 16 && fabs((*p)->pdg_id())!=25
	   && fabs((*p)->pdg_id())!=411 && fabs((*p)->pdg_id())!=3112 && fabs((*p)->pdg_id())!=3122 && fabs((*p)->pdg_id())!=3222 
	   && fabs((*p)->pdg_id())!=3312 &&  fabs((*p)->pdg_id())!=3322 && fabs((*p)->pdg_id())!=3334
	   && fabs((*p)->pdg_id())!=511 && fabs((*p)->pdg_id())!=531 && fabs((*p)->pdg_id())!=521 && fabs((*p)->pdg_id())!=5122
	   && fabs((*p)->pdg_id())!=5232 && fabs((*p)->pdg_id())!=5132 && fabs((*p)->pdg_id())!=5332 && fabs((*p)->pdg_id())!=5312
	   && fabs((*p)->pdg_id())!=533 && fabs((*p)->pdg_id())!=513 && fabs((*p)->pdg_id())!=523 && fabs((*p)->pdg_id())!=541
	   && fabs((*p)->pdg_id())!=543 && fabs((*p)->pdg_id())!=5114 && fabs((*p)->pdg_id())!=5212 && fabs((*p)->pdg_id())!=5224
	   && fabs((*p)->pdg_id())!=5222 && fabs((*p)->pdg_id())!=5214 && fabs((*p)->pdg_id())!=5112 && fabs((*p)->pdg_id())!=5324
	   && fabs((*p)->pdg_id())!=5314 && fabs((*p)->pdg_id())!=5334 && fabs((*p)->pdg_id())!=5322   )
	{
	  cout << fabs((*p)->pdg_id())<< " " << ((*p)->status()) << " error in DefinePreJets: not all particles clustered!" <<endl;
	}
      


      if ( fabs((*p)->pdg_id())==211 || fabs((*p)->pdg_id())==2212|| fabs((*p)->pdg_id())==2112  || fabs((*p)->pdg_id())==310
	   || fabs((*p)->pdg_id())==130 || fabs((*p)->pdg_id())==321 || fabs((*p)->pdg_id())==22 
	   || fabs((*p)->pdg_id())==411 || fabs((*p)->pdg_id())==3112 || fabs((*p)->pdg_id())==3122 || fabs((*p)->pdg_id())==3222 
	   || fabs((*p)->pdg_id())==3312 ||  fabs((*p)->pdg_id())==3322 || fabs((*p)->pdg_id())==3334
	   || fabs((*p)->pdg_id())==511 || fabs((*p)->pdg_id())==531 || fabs((*p)->pdg_id())==521 || fabs((*p)->pdg_id())==5122
	   || fabs((*p)->pdg_id())==5232 || fabs((*p)->pdg_id())==5132 || fabs((*p)->pdg_id())==5332 || fabs((*p)->pdg_id())==5312
	   || fabs((*p)->pdg_id())==533 || fabs((*p)->pdg_id())==513 || fabs((*p)->pdg_id())==523 || fabs((*p)->pdg_id())==541
	   || fabs((*p)->pdg_id())==543 || fabs((*p)->pdg_id())==5114 || fabs((*p)->pdg_id())==5212 || fabs((*p)->pdg_id())==5224
	   || fabs((*p)->pdg_id())==5222 || fabs((*p)->pdg_id())==5214 || fabs((*p)->pdg_id())==5112 || fabs((*p)->pdg_id())==5324
	   || fabs((*p)->pdg_id())==5314 || fabs((*p)->pdg_id())==5334 || fabs((*p)->pdg_id())==5322 )
	{
	
	  temp=PseudoJet( (*p)->momentum().px(),
			  (*p)->momentum().py(),
			  (*p)->momentum().pz(),
			  (*p)->momentum().e());
	  temp.set_user_index(((*p)->pdg_id()));
	  PartOfPrejets.push_back(  temp );
	  
	}
      
      
      if ( fabs((*p)->pdg_id()) == 11 ) 
	{
	
	  temp=PseudoJet( (*p)->momentum().px(),
			  (*p)->momentum().py(),
			  (*p)->momentum().pz(),
			  (*p)->momentum().e());
	  temp.set_user_index(((*p)->pdg_id()));
	  
	  LeptonsLoc.push_back(  temp );
	}
      
      if ( fabs((*p)->pdg_id()) == 13 )
	{
	
	  temp= PseudoJet( (*p)->momentum().px(),
			   (*p)->momentum().py(),
			   (*p)->momentum().pz(),
			   (*p)->momentum().e());
	  temp.set_user_index(((*p)->pdg_id()));
	  
	  MuonsLoc.push_back( temp);
	}
      
      
      if ( fabs((*p)->pdg_id()) == 12 || fabs((*p)->pdg_id()) == 14 || fabs((*p)->pdg_id()) == 16 
	   || fabs((*p)->pdg_id()) == 1000022 || fabs((*p)->pdg_id()) == 1000023 || fabs((*p)->pdg_id()) == 1000025 
	   || fabs((*p)->pdg_id()) == 1000035 ) 
	{
	  
	  temp=PseudoJet( (*p)->momentum().px(),
			  (*p)->momentum().py(),
			  (*p)->momentum().pz(),
			  (*p)->momentum().e());
	  
	  temp.set_user_index(((*p)->pdg_id()));
	  NeutrinosLoc.push_back( temp );
	  
	}
      
    }
  }
  

  //  cout << "finished " << Photons.size() << endl;
  // look for isolated electrons wrst hadrons incl pi0
  for (unsigned kk=0; kk<LeptonsLoc.size(); kk++)
    {
      double ethad=0;
      for (unsigned kkk=0; kkk<PartOfPrejets.size(); kkk++)
	{
	  if (Rsep(LeptonsLoc[kk],PartOfPrejets[kkk])<0.3)
	    {
	      ethad=ethad+PartOfPrejets[kkk].perp();
	    }
	}
      if (true) //ethad<0.1*LeptonsLoc[kk].perp())
	{
	  Leptons.push_back(LeptonsLoc[kk]);
	}
      else
	{
	  PartOfPrejets.push_back(LeptonsLoc[kk]);
	}
    }
  // finish with isolated muons
  for (unsigned kk=0; kk<MuonsLoc.size(); kk++)
    {
      double ethad=0;
      for (unsigned kkk=0; kkk<PartOfPrejets.size(); kkk++)
	{
	  if (Rsep(MuonsLoc[kk],PartOfPrejets[kkk])<0.3)
	    {
	      ethad=ethad+PartOfPrejets[kkk].perp();
	    }
	}
      if (true)//ethad<0.1*MuonsLoc[kk].perp())
	{
	  Muons.push_back(MuonsLoc[kk]);
	}
      else
	{
	  PartOfPrejets.push_back(MuonsLoc[kk]);
	}
    }
  
  
  PartOfPrejets2=DefineJets(PartOfPrejets,0.,5.);
  if (Leptons.size()>1) {
    Leptons=sorted_by_pt(Leptons);
  }
  
  if (Neutrinos.size()>1) {
    Neutrinos=sorted_by_pt(Neutrinos);
  }
  
  
  if (Muons.size()>1) {
    Muons=sorted_by_pt(Muons);
  }
  
  
  LeptonsLoc.clear();
  NeutrinosLoc.clear();
  PartOfPrejets.clear();
  
}




template <class T>
double SimpleAnalysisHEPMC<T>::Rsep(const fastjet::PseudoJet vec1, const fastjet::PseudoJet vec2)
{
  double phi1,phi2,y1,y2,result,dphi;
  double const pi(3.14159265358979323);

  phi1=vec1.phi_std();
  phi2=vec2.phi_std();
  y1=vec1.rap();
  y2=vec2.rap();
  dphi=(phi1-phi2);
  if (dphi<=-pi) dphi=dphi+2.*pi;

  result= sqrt( pow(dphi,2) + pow((y1-y2),2) );
  return result;
  
}

template <class T>
double SimpleAnalysisHEPMC<T>::Phisep(const fastjet::PseudoJet vec1, const fastjet::PseudoJet vec2)
{
  double phi1,phi2,result,dphi;
  double const pi(3.14159265358979323);

  phi1=vec1.phi_std();
  phi2=vec2.phi_std();
  dphi=(phi1-phi2);
  if (dphi<=-pi) dphi=dphi+2.*pi;

  result= dphi;
  return result;
  
}
