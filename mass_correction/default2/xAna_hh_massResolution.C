// example code to run Bulk Graviton->ZZ->ZlepZhad selections on electron-channel

#include <vector>
#include <iostream>
#include <fstream>
#include <algorithm>
#include <TString.h>
#include <map>
#include <TH1D.h>
#include <TFile.h>
#include <TGraph.h>
#include "../../untuplizer.h"
#include <TClonesArray.h>
#include <TLorentzVector.h>
#include <TF1.h>

using namespace std;


float getPUPPIweight(float puppipt, float puppieta ){

  double barrleMean[5]={1.1533,1.14386,1.13982,1.14885,1.16188};
  double endcapMean[5]={1.29727,1.25759,1.23614,1.23004,1.22874};
  double barrlex[5]={364.602,452.286,627.553,883.213,1310.44};
  double endcapx[5]={337.376,408.777,531.384,683.63,910.163};
 TGraph* tg1=new TGraph(5,barrlex,barrleMean);
 TGraph* tg2=new TGraph(5,endcapx,endcapMean);

  float totalWeight=1;
   if( fabs(puppieta)  <= 1.3 ){
	return tg1->Eval(puppipt);
   }
else{
	 return tg2->Eval(puppipt);
}
  // file->Close();
  //return totalWeight;
}

float getPUPPIweightBin(float puppipt, float puppieta ){

  float totalWeight=1;
   if( fabs(puppieta)  <= 1.3 ){
	   if(puppipt>200 && puppipt<=300)totalWeight=1.32257;
	   else if(puppipt>300 && puppipt<=400)totalWeight=1.20493;
	    else if(puppipt>400 && puppipt<=500)totalWeight=1.14726;
	    else if(puppipt>500 && puppipt<=600)totalWeight=1.14608;
	   else  if(puppipt>600 && puppipt<=700)totalWeight=1.14016;
	   else  if(puppipt>700 && puppipt<=800)totalWeight=1.13022;
	   else  if(puppipt>800 && puppipt<=900)totalWeight=1.12412;
	   else  if(puppipt>900 && puppipt<=1000)totalWeight=1.14825;
	  else   if(puppipt>1000 && puppipt<=1250)totalWeight=1.13804;
	   else  if(puppipt>1250 && puppipt<=1500)totalWeight=1.09761;
	   else  if(puppipt>1500 && puppipt<=1750)totalWeight=1.10551;
	   else  if(puppipt>1750)totalWeight=1.09126;
	   //else  if(puppipt>2000)totalWeight=1.10885;
   }
else{
	  if(puppipt>200 && puppipt<=300)totalWeight=1.39712;
	  else   if(puppipt>300 && puppipt<=400)totalWeight=1.27803;
	  else   if(puppipt>400 && puppipt<=500)totalWeight=1.22741;
	  else   if(puppipt>500 && puppipt<=600)totalWeight=1.18521;
	  else   if(puppipt>600 && puppipt<=700)totalWeight=1.17565;
	  else   if(puppipt>700 && puppipt<=800)totalWeight=1.14685;
	  else  if(puppipt>800 && puppipt<=900)totalWeight=1.12692;
	  else   if(puppipt>900 && puppipt<=1000)totalWeight=1.14497;
	  else   if(puppipt>1000 && puppipt<=1250)totalWeight=1.16743;
	  else   if(puppipt>1250)totalWeight=8.92857;
	  // else  if(puppipt>1500 && puppipt<=1750)totalWeight=1.1231;
	  // else  if(puppipt>1750 && puppipt<=2000)totalWeight=1.27456;
	  //    else  if(puppipt>2000)totalWeight=1.99435;
}
  // file->Close();
  return totalWeight;
}

float getPUPPIweight_o(float puppipt, float puppieta ){

   TF1* puppisd_corrGEN      = new TF1("puppisd_corrGEN","[0]+[1]*pow(x*[2],-[3])");
  puppisd_corrGEN->SetParameters(
   				 1.00626,
   				 -1.06161,
   				 0.07999,
   				 1.20454
   				 );
  TF1* puppisd_corrRECO_cen = new TF1("puppisd_corrRECO_cen","[0]+[1]*x+[2]*pow(x,2)+[3]*pow(x,3)+[4]*pow(x,4)+[5]*pow(x,5)");
  puppisd_corrRECO_cen->SetParameters(
   				      1.05807,
   				      -5.91971e-05,
   				      2.296e-07,
   				      -1.98795e-10,
   				      6.67382e-14,
   				      -7.80604e-18
   				      );

  TF1* puppisd_corrRECO_for = new TF1("puppisd_corrRECO_for","[0]+[1]*x+[2]*pow(x,2)+[3]*pow(x,3)+[4]*pow(x,4)+[5]*pow(x,5)");
  puppisd_corrRECO_for->SetParameters(
   				      1.26638,
   				      -0.000658496,
   				      9.73779e-07,
   				      -5.93843e-10,
   				      1.61619e-13,
   				      -1.6272e-17);

  float genCorr  = 1.;
  float recoCorr = 1.;
  float totalWeight = 1.;

  genCorr =  puppisd_corrGEN->Eval( puppipt );
  if( fabs(puppieta)  <= 1.3 ) recoCorr = puppisd_corrRECO_cen->Eval( puppipt );
  else
    if( fabs(puppieta) > 1.3 ) recoCorr = puppisd_corrRECO_for->Eval( puppipt );

  totalWeight = genCorr * recoCorr;
  // file->Close();
  return totalWeight;
}


void xAna_hh_massResolutionBase(std::string inputFile,TString outputFile,int pt_up=60,int pt_down=0, bool matchb=false, bool debug=false, bool cut=false){


  cout << "output file name = " << outputFile.Data() << endl;


  //get TTree from file ...
  TFile* f;
  f = TFile::Open(Form("%s",inputFile.data()));
 TDirectory * dir; 
 dir = (TDirectory*)f->Get(Form("%s:/tree",inputFile.data()));
 TTree* tree;
 dir->GetObject("treeMaker",tree);
		
		
  TreeReader data(tree);
  
  Long64_t nTotal=0;
  Long64_t nPass[20]={0};

  const int nHistos=3;

  TH1F* h_massDiff = new TH1F("h_massDiff","",100,-0.5,0.5);
  //TH1F* h_mass     = new TH1F("h_mass","",100,0,200);
  TH1F* h_mass     = new TH1F("h_mass","",100,62.5,187.5);

  TH1F* h_SD[nHistos];
  TH1F* h_SDCorr[nHistos];
  TH1F* h_SDCorrThea[nHistos];
  
  TH1F* h_AK8SD[nHistos];
  TH1F* h_AK8SDCorrThea[nHistos];
  TH1F* h_AK8SDHCorr[nHistos];
  TH1F* h_AK8SDThealikeHCorr[nHistos];
  
  TH1F* h_PR[nHistos];
  TH1F* h_PRCorr[nHistos];

  TH1F* h_diff_SD[nHistos];
  TH1F* h_diff_SDCorr[nHistos];
  TH1F* h_diff_SDCorrThea[nHistos];

  TH1F* h_diff_PR[nHistos];
  TH1F* h_diff_PRCorr[nHistos];
  
  TH1F* h_diff_AK8SD[nHistos];
  TH1F* h_diff_AK8SDCorrThea[nHistos];
  TH1F* h_diff_AK8SDHCorr[nHistos];
  TH1F* h_diff_AK8SDThealikeHCorr[nHistos];
  std::string prefix[]={"leading","subleading","both"};
  for(int i=0; i<nHistos; i++)
    {

      h_SD[i] = (TH1F*)h_mass->Clone(Form("h_SD_%s",prefix[i].data()));
      h_SD[i]->SetXTitle("Raw Puppi+Softdrop mass [GeV]");
	
	h_AK8SD[i] = (TH1F*)h_mass->Clone(Form("h_AK8SD_%s",prefix[i].data()));
      h_AK8SD[i]->SetXTitle("AK8 Raw Puppi+Softdrop mass [GeV]");

      h_SDCorr[i] = (TH1F*)h_mass->Clone(Form("h_SDCorr_%s",prefix[i].data()));
      h_SDCorr[i]->SetXTitle("L2L3-corrected Puppi+Softdrop mass [GeV]");

      h_SDCorrThea[i] = (TH1F*)h_mass->Clone(Form("h_SDCorrThea_%s",prefix[i].data()));
      h_SDCorrThea[i]->SetXTitle("Thea-corrected Puppi+Softdrop mass [GeV]");
	
	h_AK8SDCorrThea[i] = (TH1F*)h_mass->Clone(Form("h_AK8SDCorrThea_%s",prefix[i].data()));
      h_AK8SDCorrThea[i]->SetXTitle("AK8 Thea-corrected Puppi+Softdrop mass [GeV]");
	
	h_AK8SDHCorr[i] = (TH1F*)h_mass->Clone(Form("h_AK8SDHCorr_%s",prefix[i].data()));
      h_AK8SDHCorr[i]->SetXTitle("AK8 H-corrected Puppi+Softdrop mass [GeV]");
	
	h_AK8SDThealikeHCorr[i] = (TH1F*)h_mass->Clone(Form("h_AK8SDThealikeHCorr_%s",prefix[i].data()));
      h_AK8SDThealikeHCorr[i]->SetXTitle("AK8 Thealike H-corrected Puppi+Softdrop mass [GeV]");

      h_PR[i] = (TH1F*)h_mass->Clone(Form("h_PR_%s",prefix[i].data()));
      h_PR[i]->SetXTitle("Raw CHS+Pruned mass [GeV]");

      h_PRCorr[i] = (TH1F*)h_mass->Clone(Form("h_PRCorr_%s",prefix[i].data()));
      h_PRCorr[i]->SetXTitle("L2L3-corrected CHS+Pruned mass [GeV]");

      // study the difference with respect to 125 GeV

      h_diff_SD[i] = (TH1F*)h_massDiff->Clone(Form("h_diff_SD_%s",prefix[i].data()));
      h_diff_SD[i]->SetXTitle("Raw Puppi+Softdrop (m-125)/125");
	
	h_diff_AK8SD[i] = (TH1F*)h_massDiff->Clone(Form("h_diff_AK8SD_%s",prefix[i].data()));
      h_diff_AK8SD[i]->SetXTitle("AK8 Raw Puppi+Softdrop (m-125)/125");

      h_diff_SDCorr[i] = (TH1F*)h_massDiff->Clone(Form("h_diff_SDCorr_%s",prefix[i].data()));
      h_diff_SDCorr[i]->SetXTitle("L2L3-corrected Puppi+Softdrop (m-125)/125");

      h_diff_SDCorrThea[i] = (TH1F*)h_massDiff->Clone(Form("h_diff_SDCorrThea_%s",prefix[i].data()));
      h_diff_SDCorrThea[i]->SetXTitle("Thea-corrected Puppi+Softdrop (m-125)/125");
	
	h_diff_AK8SDCorrThea[i] = (TH1F*)h_massDiff->Clone(Form("h_diff_AK8SDCorrThea_%s",prefix[i].data()));
      h_diff_AK8SDCorrThea[i]->SetXTitle("AK8 Thea-corrected Puppi+Softdrop (m-125)/125");
	
	h_diff_AK8SDHCorr[i] = (TH1F*)h_massDiff->Clone(Form("h_diff_AK8SDHCorr_%s",prefix[i].data()));
      h_diff_AK8SDHCorr[i]->SetXTitle("AK8 H-corrected Puppi+Softdrop (m-125)/125");
	
	h_diff_AK8SDThealikeHCorr[i] = (TH1F*)h_massDiff->Clone(Form("h_diff_AK8SDThealikeHCorr_%s",prefix[i].data()));
      h_diff_AK8SDThealikeHCorr[i]->SetXTitle("AK8 Thealike H-corrected Puppi+Softdrop (m-125)/125");

      h_diff_PR[i] = (TH1F*)h_massDiff->Clone(Form("h_diff_PR_%s",prefix[i].data()));
      h_diff_PR[i]->SetXTitle("Raw CHS+Pruned (m-125)/125");

      h_diff_PRCorr[i] = (TH1F*)h_massDiff->Clone(Form("h_diff_PRCorr_%s",prefix[i].data()));
      h_diff_PRCorr[i]->SetXTitle("L2L3-corrected CHS+Pruned (m-125)/125");

    }


  for(Long64_t jEntry=0; jEntry<data.GetEntriesFast() ;jEntry++){

  if((jEntry+1)%2)continue;
    if (jEntry % 1000 == 0)
      fprintf(stderr, "Processing event %lli of %lli\n", jEntry + 1, data.GetEntriesFast());


    if(debug && jEntry>10)break;

    data.GetEntry(jEntry);
    nTotal++;

    //2. pass electron or muon trigger                                                                                                                                          
    std::string* trigName = data.GetPtrString("hlt_trigName");
    vector<bool> &trigResult = *((vector<bool>*) data.GetPtr("hlt_trigResult"));
    const Int_t nsize = data.GetPtrStringSize();

    bool passTrigger=false;
    for(int it=0; it< nsize; it++)
      {
	std::string thisTrig= trigName[it];
        bool results = trigResult[it];

        // std::cout << thisTrig << " : " << results << std::endl;                                                                                                              

        if( (thisTrig.find("HLT_PFHT900_v")!= std::string::npos && results==1)
            )
          {
            passTrigger=true;
            break;
          }


      }

    if(!passTrigger && cut)continue;
    nPass[4]++;



    Int_t nGenPar        = data.GetInt("nGenPar");
    Int_t* genParId      = data.GetPtrInt("genParId");
    Int_t* genParSt      = data.GetPtrInt("genParSt");
    Int_t* genMomParId   = data.GetPtrInt("genMomParId");
    Int_t* genMo1      = data.GetPtrInt("genMo1");
    Int_t* genDa1      = data.GetPtrInt("genDa1");
    Int_t* genDa2      = data.GetPtrInt("genDa2");
    TClonesArray* genParP4 = (TClonesArray*) data.GetPtrTObject("genParP4");

    int genHIndex[2]={-1,-1};
    int genbIndex[2][2]={{-1,-1},
			 {-1,-1}};		       

    std::vector<int> neutrinos;
    neutrinos.clear();
    std::vector<int> Higgses;
    Higgses.clear();
    std::vector<int> bhadrons;
    bhadrons.clear();
    for(int ig=0; ig < nGenPar; ig++){

      int pid   = abs(genParId[ig]);
      int mpid  = abs(genMomParId[ig]);

      if( (pid==12 || pid==14 || pid==16) && genParSt[ig]==1
	  )
	{
			      
	  // Now start checking the mother all the way and make sure we only 
	  // consider neutrinos from Higgs
	  bool hasHMother=false;
	  int mpid=-9999;
	  int dindex=ig;
	  int dpid=-9999;
	  int mindex=-1;
	  int bquarkIndex=-1;
	  int bhadronIndex=-1;
	  while(!hasHMother && dindex>=0 && abs(mpid)!=2212){
			      
	    mpid = genMomParId[dindex];
	    mindex= genMo1[dindex];
	    dpid = genParId[dindex];
	    //				cout << ig << "\t" << dindex << "\t" << mindex << "\t" << mpid << "\t" << dpid << endl;
	    if( mpid==25)
	      {
		hasHMother=true;
		bquarkIndex=dindex;
	      }
	    else
	      {
		dindex=mindex;
		hasHMother=false;
		if((abs(mpid)>=500 && abs(mpid)<600) ||
		   (abs(mpid)>=5000 && abs(mpid)<6000))
		  {
		    bhadronIndex=mindex;
		  }
	      }
	  }
			    
			      
	  if(!hasHMother)continue;
	  if(bhadronIndex<0)continue;
	  if(bquarkIndex<0)continue;

	  TLorentzVector* nu_l4 = (TLorentzVector*)genParP4->At(ig);
	  TLorentzVector* b_l4  = (TLorentzVector*)genParP4->At(bquarkIndex);
	  TLorentzVector* B_l4  = (TLorentzVector*)genParP4->At(bhadronIndex);

	  neutrinos.push_back(ig);
	  Higgses.push_back(mindex);
	  bhadrons.push_back(bhadronIndex);
	}

      if(genParId[ig]!=25)continue;

      if(genHIndex[0]<0)
	{
	  genHIndex[0]=ig;
	  genbIndex[0][0]=genDa1[ig];
	  genbIndex[0][1]=genDa2[ig];
	}

      else if(genHIndex[1]<0)
	{
	  genHIndex[1]=ig;
	  genbIndex[1][0]=genDa1[ig];
	  genbIndex[1][1]=genDa2[ig];
	}

    }    


    if(genHIndex[0]<0 || genHIndex[1]<0)continue;
    if(genbIndex[0][0]<0 || genbIndex[0][1]<0)continue;
    if(genbIndex[1][0]<0 || genbIndex[1][1]<0)continue;

    nPass[0]++;

    if(genHIndex[0]==genHIndex[1])continue;
    nPass[1]++;

    TLorentzVector genH_l4[2];
    TLorentzVector genb_l4[2][2];

    for(int ih=0; ih<2; ih++)
      {
	genH_l4[ih] = *((TLorentzVector*)genParP4->At(genHIndex[ih]));
	for(int ib=0; ib<2; ib++)
	  {
	    genb_l4[ih][ib] = *((TLorentzVector*)genParP4->At(genbIndex[ih][ib]));
	  }
      }


    if(debug){
      cout << genHIndex[0] << "\t" << genHIndex[1] << endl;
      genH_l4[0].Print();
      genH_l4[1].Print();

      cout << genbIndex[0][0] << "\t" << genbIndex[0][1] << "\t" 
	   << genbIndex[1][0] << "\t" << genbIndex[1][1] << endl;
      genH_l4[0].Print();
      genH_l4[1].Print();
      genb_l4[0][0].Print();
      genb_l4[0][1].Print();
      genb_l4[1][0].Print();
      genb_l4[1][1].Print();

    }
    int nFATJet         = data.GetInt("FATnJet");
    const int nJets=nFATJet;
    TClonesArray* fatjetP4   = (TClonesArray*) data.GetPtrTObject("FATjetP4");
    TClonesArray* puppijetP4 = (TClonesArray*) data.GetPtrTObject("FATjetPuppiP4");
    TClonesArray* AK8PuppijetP4 = (TClonesArray*) data.GetPtrTObject("AK8PuppijetP4");

    // check matching first

    bool findAMatch=false;
    const float dRMax=0.4;
    const float dRbMax=0.8;
    int matchedHJetIndex[2]={-1,-1};
    int matchedHiggsIndex[2]={-1,-1};

    for(int ij=0; ij<nJets; ij++)
      {
	TLorentzVector* thisJet = (TLorentzVector*)fatjetP4->At(ij);

	for(int jj=0; jj<nJets; jj++)
	  {

	    if(ij==jj)continue;
	    TLorentzVector* thatJet = (TLorentzVector*)fatjetP4->At(jj);
	    
	    if(thisJet->DeltaR(genH_l4[0])<dRMax && 
	       (!matchb || (matchb && 
			    thisJet->DeltaR(genb_l4[0][0])<dRbMax && 
			    thisJet->DeltaR(genb_l4[0][1])<dRbMax)) &&
	       thatJet->DeltaR(genH_l4[1])<dRMax &&
	       (!matchb || (matchb &&
			    thatJet->DeltaR(genb_l4[1][0])<dRbMax &&
			    thatJet->DeltaR(genb_l4[1][1])<dRbMax)))
	      {
		if(debug)
		  {
		    cout << "dRhb00= " <<  thisJet->DeltaR(genb_l4[0][0]) << endl;
		    cout << "dRhb01= " <<  thisJet->DeltaR(genb_l4[0][1]) << endl;
		    cout << "dRhb10= " <<  thatJet->DeltaR(genb_l4[1][0]) << endl;
		    cout << "dRhb11= " <<  thatJet->DeltaR(genb_l4[1][1]) << endl;
		  }
		if(ij<jj){
		  matchedHJetIndex[0]=ij;
		  matchedHJetIndex[1]=jj;
		  matchedHiggsIndex[0]=genHIndex[0];
		  matchedHiggsIndex[1]=genHIndex[1];
		}
		else
		  {
		    matchedHJetIndex[0]=jj;
		    matchedHJetIndex[1]=ij;
		    matchedHiggsIndex[0]=genHIndex[1];
		    matchedHiggsIndex[1]=genHIndex[0];
		  }
		findAMatch=true;
		break;
	      }

	    if(findAMatch)break;

	  }	

	if(findAMatch)break;

      }

    if(!findAMatch)continue;
    if(debug)
      cout << matchedHJetIndex[0] << "\t" << matchedHJetIndex[1] << endl;
    nPass[2]++;

    
    bool findAK8Match=false;
    //const float dRMax=0.4;
    //const float dRbMax=0.8;
    int matchedHAK8JetIndex[2]={-1,-1};
		      int AK8nJet=data.GetInt("AK8PuppinJet");
    for(int ij=0; ij<AK8nJet; ij++)
      {
	TLorentzVector* thisJet = (TLorentzVector*)AK8PuppijetP4->At(ij);

	for(int jj=0; jj<AK8nJet; jj++)
	  {

	    if(ij==jj)continue;
	    TLorentzVector* thatJet = (TLorentzVector*)AK8PuppijetP4->At(jj);
	    
	    if(thisJet->DeltaR(genH_l4[0])<dRMax && 
	       (!matchb || (matchb && 
			    thisJet->DeltaR(genb_l4[0][0])<dRbMax && 
			    thisJet->DeltaR(genb_l4[0][1])<dRbMax)) &&
	       thatJet->DeltaR(genH_l4[1])<dRMax &&
	       (!matchb || (matchb &&
			    thatJet->DeltaR(genb_l4[1][0])<dRbMax &&
			    thatJet->DeltaR(genb_l4[1][1])<dRbMax)))
	      {
		if(debug)
		  {
		    cout << "dRhb00= " <<  thisJet->DeltaR(genb_l4[0][0]) << endl;
		    cout << "dRhb01= " <<  thisJet->DeltaR(genb_l4[0][1]) << endl;
		    cout << "dRhb10= " <<  thatJet->DeltaR(genb_l4[1][0]) << endl;
		    cout << "dRhb11= " <<  thatJet->DeltaR(genb_l4[1][1]) << endl;
		  }
		if(ij<jj){
		  matchedHAK8JetIndex[0]=ij;
		  matchedHAK8JetIndex[1]=jj;
		}
		else
		  {
		    matchedHAK8JetIndex[0]=jj;
		    matchedHAK8JetIndex[1]=ij;
		  }
		findAK8Match=true;
		break;
	      }

	    if(findAK8Match)break;

	  }	

	if(findAK8Match)break;

      }

    if(!findAK8Match)continue;
   
    //0. has a good vertex
    Int_t nVtx        = data.GetInt("nVtx");
    if(nVtx<1 && cut)continue;
    nPass[3]++;



    Float_t*  fatjetTau1 = data.GetPtrFloat("FATjetTau1");
    Float_t*  fatjetTau2 = data.GetPtrFloat("FATjetTau2");
    Float_t*  fatjetCISVV2 = data.GetPtrFloat("FATjetCISVV2");
    Float_t*  fatjetPRmass = data.GetPtrFloat("FATjetPRmass");
    Float_t*  fatjetPRmassL2L3Corr = data.GetPtrFloat("FATjetPRmassL2L3Corr");
    Float_t*  fatjetSDmass = data.GetPtrFloat("FATjetPuppiSDmass");
    Float_t*  fatjetSDmassL2L3Corr = data.GetPtrFloat("FATjetPuppiSDmassL2L3Corr"); 
    Float_t*  AK8PuppijetSDmass = data.GetPtrFloat("AK8PuppijetSDmass");

    vector<bool>    &passFatJetLooseID = *((vector<bool>*) data.GetPtr("FATjetPassIDLoose"));
    
    TLorentzVector recoH_l4[2];
    int nGoodJets=0;
    
    
    
    for(int i=0; i<2; i++)
      {
    	
	int ij = matchedHJetIndex[i];
     	TLorentzVector* thisJet = (TLorentzVector*)fatjetP4->At(ij);
	recoH_l4[i]= (*thisJet);
    	if(thisJet->Pt()<200)continue;
	if(fabs(thisJet->Eta())>2.4)continue;
	nGoodJets++;
      }

    if(nGoodJets<2)continue;
    nPass[5]++;
    
    if(debug)
      {
	recoH_l4[0].Print();
	recoH_l4[1].Print();
      }

    float dEta=fabs(recoH_l4[0].Eta()-recoH_l4[1].Eta());
    if(dEta>1.3 && cut)continue;
    nPass[6]++;

    float M=(recoH_l4[0] + recoH_l4[1]).M();
    if(M<800 && cut)continue;
    nPass[7]++;


    int nHP=0;
    int nLP=0;
    for(int i=0; i<2; i++)
      {
    	
	int ij = matchedHJetIndex[i];

	float tau21_i = fatjetTau2[ij]/fatjetTau1[ij];
	bool isHP= (tau21_i < 0.6);
	if(isHP)nHP++;
      }

    if(nHP<2 && cut)continue;
    nPass[8]++;

	long eventId=data.GetLong64("eventId");
	
	if(AK8nJet<2)continue;
    // now plot mass
    for(int i=0; i<2;i++)
      {
		
		//cout<< eventId<<endl;
	int jet=matchedHJetIndex[i];
	int AK8jet=matchedHAK8JetIndex[i];
	TLorentzVector* thisJet = (TLorentzVector*)puppijetP4->At(jet);
	float thea_corr = getPUPPIweight(thisJet->Pt(),thisJet->Eta());
	float thea_mass = fatjetSDmass[jet]*thea_corr;
	
	

	//if(jet>AK8nJet-1)break;
	if(thisJet->Pt()>99998)break;
	//cout<< eventId<<endl;
	
	TLorentzVector* thisAK8Jet = (TLorentzVector*)AK8PuppijetP4->At(AK8jet);
	thea_corr = getPUPPIweight_o(thisAK8Jet->Pt(),thisAK8Jet->Eta());
	
	float H_corr= getPUPPIweight(thisAK8Jet->Pt(),thisAK8Jet->Eta());
	float H_corrBin= getPUPPIweightBin(thisAK8Jet->Pt(),thisAK8Jet->Eta());
	
	if(thisAK8Jet->Pt()<pt_down || thisAK8Jet->Pt()>pt_up )continue;

	if(debug)
	  cout << thisJet->Pt() << "\t" << thisJet->Eta() << "\t" << thea_corr << endl;
	
	h_SD[i]->Fill(fatjetSDmass[jet]);
	h_SDCorr[i]->Fill(fatjetSDmassL2L3Corr[jet]);
	h_SDCorrThea[i]->Fill(thea_mass);
	h_PR[i]->Fill(fatjetPRmass[jet]);
	h_PRCorr[i]->Fill(fatjetPRmassL2L3Corr[jet]);
	
	h_AK8SD[i]->Fill(AK8PuppijetSDmass[AK8jet]);
	h_AK8SDCorrThea[i]->Fill(AK8PuppijetSDmass[AK8jet]*thea_corr);
	
	h_AK8SDHCorr[i]->Fill(AK8PuppijetSDmass[AK8jet]*H_corr);
	h_AK8SDThealikeHCorr[i]->Fill(AK8PuppijetSDmass[AK8jet]*H_corr);

	h_SD[2]->Fill(fatjetSDmass[jet]);
	h_SDCorr[2]->Fill(fatjetSDmassL2L3Corr[jet]);
	h_SDCorrThea[2]->Fill(thea_mass);
	h_PR[2]->Fill(fatjetPRmass[jet]);
	h_PRCorr[2]->Fill(fatjetPRmassL2L3Corr[jet]);
	
	h_AK8SD[2]->Fill(AK8PuppijetSDmass[AK8jet]);
	h_AK8SDCorrThea[2]->Fill(AK8PuppijetSDmass[AK8jet]*thea_corr);
	h_AK8SDHCorr[2]->Fill(AK8PuppijetSDmass[AK8jet]*H_corr);
	h_AK8SDThealikeHCorr[2]->Fill(AK8PuppijetSDmass[AK8jet]*H_corr);

	h_diff_SD[i]->Fill((fatjetSDmass[jet]-125)/125);
	h_diff_SDCorr[i]->Fill((fatjetSDmassL2L3Corr[jet]-125)/125);
	h_diff_SDCorrThea[i]->Fill((thea_mass-125)/125);
	h_diff_PR[i]->Fill((fatjetPRmass[jet]-125)/125);
	h_diff_PRCorr[i]->Fill((fatjetPRmassL2L3Corr[jet]-125)/125);
	
	h_diff_AK8SD[i]->Fill((AK8PuppijetSDmass[AK8jet]-125)/125);
	h_diff_AK8SDCorrThea[i]->Fill((AK8PuppijetSDmass[AK8jet]*thea_corr-125)/125);
	
	h_diff_AK8SDHCorr[i]->Fill((AK8PuppijetSDmass[AK8jet]*H_corr-125)/125);
	h_diff_AK8SDThealikeHCorr[i]->Fill((AK8PuppijetSDmass[AK8jet]*H_corr-125)/125);

	h_diff_SD[2]->Fill((fatjetSDmass[jet]-125)/125);
	h_diff_SDCorr[2]->Fill((fatjetSDmassL2L3Corr[jet]-125)/125);
	h_diff_SDCorrThea[2]->Fill((thea_mass-125)/125);
	h_diff_PR[2]->Fill((fatjetPRmass[jet]-125)/125);
	h_diff_PRCorr[2]->Fill((fatjetPRmassL2L3Corr[jet]-125)/125);
	
	h_diff_AK8SD[2]->Fill((AK8PuppijetSDmass[AK8jet]-125)/125);
	h_diff_AK8SDCorrThea[2]->Fill((AK8PuppijetSDmass[AK8jet]*thea_corr-125)/125);
	
	h_diff_AK8SDHCorr[2]->Fill((AK8PuppijetSDmass[AK8jet]*H_corr-125)/125);
	h_diff_AK8SDThealikeHCorr[2]->Fill((AK8PuppijetSDmass[AK8jet]*H_corr-125)/125);

      }
    

  } // end of loop over entries

  std::cout << "nTotal    = " << nTotal << std::endl;
  for(int i=0;i<20;i++)
    if(nPass[i]>0)
      std::cout << "nPass[" << i << "]= " << nPass[i] << std::endl;

  TFile* outFile = new TFile(Form("%d/%s",pt_down,outputFile.Data()),"recreate");

  for(int i=0; i<nHistos; i++)
    {
      h_diff_SD[i]->Write();
      h_diff_SDCorr[i]->Write();
      h_diff_SDCorrThea[i]->Write();
      h_diff_PR[i]->Write();
      h_diff_PRCorr[i]->Write();

      h_SD[i]->Write();
      h_SDCorr[i]->Write();
      h_SDCorrThea[i]->Write();
      h_PR[i]->Write();
      h_PRCorr[i]->Write();
	
	h_AK8SD[i]->Write();
	h_AK8SDCorrThea[i]->Write();
	h_AK8SDHCorr[i]->Write();
	h_AK8SDThealikeHCorr[i]->Write();
	
	h_diff_AK8SD[i]->Write();
	h_diff_AK8SDCorrThea[i]->Write();
	
	h_diff_AK8SDHCorr[i]->Write();
	h_diff_AK8SDThealikeHCorr[i]->Write();
    }
    
   

  outFile->Close();

 for(int i=0; i<nHistos; i++)
    {
      delete h_diff_SD[i] ;
      delete h_diff_SDCorr[i] ;
      delete h_diff_SDCorrThea[i] ;
      delete h_diff_PR[i] ;
      delete h_diff_PRCorr[i] ;

      delete h_SD[i] ;
      delete h_SDCorr[i] ;
      delete h_SDCorrThea[i] ;
      delete  h_PR[i] ;
      delete h_PRCorr[i] ;
	
	delete h_AK8SD[i] ;
	delete h_AK8SDCorrThea[i] ;
	delete h_AK8SDHCorr[i] ;
	delete h_AK8SDThealikeHCorr[i] ;
	
	delete h_diff_AK8SD[i] ;
	delete h_diff_AK8SDCorrThea[i] ;
	
	delete h_diff_AK8SDHCorr[i] ;
	delete h_diff_AK8SDThealikeHCorr[i] ;
    }
    delete h_massDiff;
    delete h_mass;
    delete outFile;
}

void xAna_hh_massResolution(int pt_down=60,int pt_up=0){
  xAna_hh_massResolutionBase("/data7/syu/special_study/80X_miniAOD/80X_neutrino/BulkGravTohhTohbbhbb/B800.root","B800.root",pt_up,pt_down);
  xAna_hh_massResolutionBase("/data7/syu/special_study/80X_miniAOD/80X_neutrino/BulkGravTohhTohbbhbb/B1000.root","B1000.root",pt_up,pt_down);
  xAna_hh_massResolutionBase("/data7/syu/special_study/80X_miniAOD/80X_neutrino/BulkGravTohhTohbbhbb/B1400.root","B1400.root",pt_up,pt_down);
  xAna_hh_massResolutionBase("/data7/syu/special_study/80X_miniAOD/80X_neutrino/BulkGravTohhTohbbhbb/B2000.root","B2000.root",pt_up,pt_down);
  xAna_hh_massResolutionBase("/data7/syu/special_study/80X_miniAOD/80X_neutrino/BulkGravTohhTohbbhbb/B3000.root","B3000.root",pt_up,pt_down);



}
/*
void xAna_hh_massResolution(){
	//xAna_hh_massResolutionPt(200,300);
	xAna_hh_massResolutionPt(0,5000);
	xAna_hh_massResolutionPt(300,400);
	xAna_hh_massResolutionPt(400,500);
	xAna_hh_massResolutionPt(500,600);
	xAna_hh_massResolutionPt(600,700);
	xAna_hh_massResolutionPt(700,800);
	xAna_hh_massResolutionPt(800,900);
	xAna_hh_massResolutionPt(900,1000);
	xAna_hh_massResolutionPt(1250,1500);
	xAna_hh_massResolutionPt(1500,5000);
	
}
*/
