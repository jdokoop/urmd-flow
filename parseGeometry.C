//-----------------------------------------------
// Code to parse the OSCAR199A output file from
// Ultrarelativistic Quantum Molecular Dynamics
// and calculate anisotropic flow coefficients
//   (1) with respect to the event plane
//   (2) with four-particle cumulants
//
// Authors: J. Orjuela-Koop
//          P. Yin
//----------------------------------------------

#include <iostream>
#include <iomanip>
#include <fstream>
#include <string>
#include <sstream>
#include <vector>
#include <numeric>
#include <iterator>

#include "TLorentzVector.h"
#include "TFile.h"
#include "TStyle.h"
#include "TNtuple.h"
#include "TCanvas.h"
#include "TH1.h"
#include "TH2.h"
#include "TH3.h"
#include "TProfile.h"

using namespace std;

//-------------------------------------------
// Variables
//-------------------------------------------

//Data structure to contain particles
struct particle
{
	float px;
	float py;
	float pz;
	float x;
	float y;
	float z;
	float phi;
	float rsquare;
	float pT;
	float eta;
};

//Complex number 
struct cmpx
{
	float re;
	float im;
};

//Number of bins in histograms
const int NOB = 100;

//Number of nucleons in system
// --> p+Au = 198
// --> d+Au = 199
// --> d+Pb = 210
// --> p+Pb = 209
const int NUCL = 199;

//Event characterization variables
int npart          = 0;
int npartsum       = 0;
int nspectator     = 0;
int numevent       = 0;
int process_number = 0;

//Variables for reference four-particle cumulant analysis
float q2x  = 0;
float q2y  = 0;
float q4x  = 0;
float q4y  = 0;
int n      = 0;
float cn_2 = 0;
float cn_4 = 0;
float db_2 = 0;
float db_4 = 0;

//Variables for differential four-particle cumulant analysis
TH1F *hpT_Template;
TProfile *hdb_2prime;
TProfile *hdb_4prime;
TProfile *hdb_2;
TProfile *hdb_4;
TH1F *hd_2;
TH1F *hv2_4;

//Vectors to store event-level geometric parameters
vector<float> epsilon2;
vector<float> epsilon3;
vector<float> psi2;
vector<float> psi3;

//Participant particles in each event
vector<particle> collisionparticles;

//Final state hadrons in each event
vector<particle> finalparticles;

//Histograms for plotting
TH1F *hPsi2;
TProfile *hv2_pT;
TProfile *hv3_pT;

//Diagnostic histograms
TH1F *hEta;
TH1F *hpT;
TH1F *hpT_EtaCut;
TH1F *hPhi;

//Miscellaneous options
bool parse_verbosity = false;

//---------------------------------------------------------
// Functions
//---------------------------------------------------------

cmpx multiplyComplex(float x, float y, float u, float v)
{
	//z1 = x + iy
	//z2 = u + iv
	cmpx c;
	c.re = x*u -y*v;
	c.im = x*v + y*u;

	return c;
}

float getNormSq(float x, float y)
{
	//z = x + iy
	return x*x + y*y;
}

void draw()
{
  gStyle->SetOptStat(0);
  TCanvas *c1 = new TCanvas("c1","c1",600,600);
  hEta->Draw();

  TCanvas *c2 = new TCanvas("c2","c2",600,600);
  hpT->Draw();

  TCanvas *c3 = new TCanvas("c3","c3",600,600);
  hPhi->Draw();
}

void writeData()
{	
	TFile *fout = new TFile(Form("/direct/phenix+hhj/jdok/UrQMD_Output/urqmd_out_%i.root",process_number),"RECREATE");
	hdb_2prime->Write();
	hdb_4prime->Write();
	hdb_2->Write();
	hdb_4->Write();
	fout->Close();
}

void computeCumulants()
{
	db_2 = hdb_2->GetBinContent(1);
	db_4 = hdb_4->GetBinContent(1);

	cn_2 = db_2;
	cn_4 = db_4 - 2*pow(db_2,2);

	for(int i=1; i<=hpT_Template->GetNbinsX(); i++)
	{
		float cont = (hdb_4prime->GetBinContent(i)) - 2*(hdb_2prime->GetBinContent(i))*db_2;
		hd_2->SetBinContent(i, cont);
	}

	//cout << "<<2>> = " << db_2 << endl;
	//cout << "<<4>> = " << db_4 << endl;
	
	//Finally, compute v2
	float scaleFactor = -1*pow(-1*cn_4,3.0/4.0);
	hd_2->Scale(1.0/scaleFactor);

	//Two-particle cumulant result
	hdb_2prime->Scale(1.0/sqrt(cn_2));
	//hdb_4prime->Draw();

	//Two particle correlation result
	//hv2_4 = (TH1F*) hd_2->Clone("hv2_4");
	//hd_2->Draw();

	//Four particle correlation result
	//hd_2->Draw();
}

void processEventCumulants()
{
	//Initialize variables
	n   = 0;
	q2x = 0;
	q2y = 0;
	q4x = 0;
	q4y = 0;

	//Compute reference flow for the event at hand
	for(int i=0; i<finalparticles.size(); i++)
	{
		float eta = finalparticles[i].eta;
		float phi = finalparticles[i].phi;
		float pT  = finalparticles[i].pT;

		hEta->Fill(eta);
		hpT->Fill(pT);
		hPhi->Fill(phi);

		//Only consider particles within -2 < eta < 2 and 0 < pT [GeV/c] < 5
		if(TMath::Abs(eta) > 2.0 || pT > 5.0) continue;

		hpT_EtaCut->Fill(eta);

		q2x = q2x + TMath::Cos(2*phi);
		q2y = q2y + TMath::Sin(2*phi);

		q4x = q4x + TMath::Cos(4*phi);
		q4y = q4y + TMath::Sin(4*phi);

		cout << finalparticles[i].pz << endl;

		n++;
	}

	float sb_2 = (getNormSq(q2x,q2y) - n)/(n*(n-1));

	cmpx c1 = multiplyComplex(q2x, -1*q2y, q2x, -1*q2y); 
	cmpx c2 = multiplyComplex(q4x, q4y, c1.re, c1.im);
	float sb_4 = (getNormSq(q2x,q2y)*getNormSq(q2x,q2y) + getNormSq(q4x,q4y) - 2*c2.re - 4*(n-2)*getNormSq(q2x,q2y) + 2*n*(n-3))/(n*(n-1)*(n-2)*(n-3));

	cout << "<2> = " << sb_2 << endl;
	cout << "<4> = " << sb_4 << endl << endl;

	hdb_2->Fill(0.5,sb_2,1);
	hdb_4->Fill(0.5,sb_4,1);

	//Compute differential flow for the event at hand
	TH1F *hp2x         = (TH1F*) hpT_Template->Clone("hp2x");
	TH1F *hp2y         = (TH1F*) hpT_Template->Clone("hp2y");
	TH1F *hp4x         = (TH1F*) hpT_Template->Clone("hp4x");
	TH1F *hp4y         = (TH1F*) hpT_Template->Clone("hp4y");
	TH1F *hsb_2prime   = (TH1F*) hpT_Template->Clone("hsb_2prime");
	TH1F *hsb_4prime   = (TH1F*) hpT_Template->Clone("hsb_4prime");

	float m[NOB] = {0};

	for(int i=0; i<finalparticles.size(); i++)
	{
		float eta = finalparticles[i].eta;
		float phi = finalparticles[i].phi;
		float pT = finalparticles[i].pT;

		if(TMath::Abs(eta) > 2.0 || pT > 5.0) continue;

		int bin = hpT_Template->FindBin(pT);
		hp2x->SetBinContent(bin, hp2x->GetBinContent(bin) + TMath::Cos(2*phi));
		hp2y->SetBinContent(bin, hp2y->GetBinContent(bin) + TMath::Sin(2*phi));
		hp4x->SetBinContent(bin, hp4x->GetBinContent(bin) + TMath::Cos(4*phi));
		hp4y->SetBinContent(bin, hp4y->GetBinContent(bin) + TMath::Sin(4*phi));

		m[bin-1]++;
	}
 
	for(int i=1; i<= hsb_2prime->GetNbinsX(); i++)
	{
		cmpx aux1 = multiplyComplex(hp2x->GetBinContent(i), hp2y->GetBinContent(i), q2x, -1*q2y);
		float aux2 = m[i-1]*(n-1); 
		hsb_2prime->SetBinContent(i, (aux1.re-m[i-1])/aux2);
	}

	for(int i=1; i<hdb_4prime->GetNbinsX(); i++)
	{
		cmpx aux1_0 = multiplyComplex(hp2x->GetBinContent(i), hp2y->GetBinContent(i), q2x, q2y);
		cmpx aux1_1 = multiplyComplex(q2x, -1*q2y, q2x, -1*q2y);
		cmpx aux1_2 = multiplyComplex(aux1_0.re, aux1_0.im, aux1_1.re, aux1_1.im);
		float aux1  = aux1_2.re;

		cmpx aux2_0 = multiplyComplex(hp4x->GetBinContent(i), hp4y->GetBinContent(i), q2x, -1*q2y);
		cmpx aux2_1 = multiplyComplex(aux2_0.re, aux2_0.im, q2x, -1*q2y);
		float aux2  = aux2_1.re;

		cmpx aux3_0 = multiplyComplex(hp2x->GetBinContent(i), hp2y->GetBinContent(i), q2x, q2y);
		cmpx aux3_1 = multiplyComplex(aux3_0.re, aux3_0.im, q4x, -1*q4y);
		float aux3  = aux2_1.re;

		cmpx aux4_0 = multiplyComplex(hp2x->GetBinContent(i), hp2y->GetBinContent(i), q2y, -1*q2y);
		float aux4  = 2*n*aux4_0.re;

		float aux5  = 2*m[i]*getNormSq(q2x,q2y);

		cmpx aux6_0 = multiplyComplex(hp2x->GetBinContent(i), hp2y->GetBinContent(i), q2x, -1*q2y);
		float aux6  = 7*aux6_0.re;

		cmpx aux7_0 = multiplyComplex(q2x, q2y, hp2x->GetBinContent(i), -1*hp2y->GetBinContent(i));
		float aux7  = aux7_0.re;

		cmpx aux8_0 = multiplyComplex(hp4x->GetBinContent(i), hp4y->GetBinContent(i), q4x, -1*q4y);
		float aux8  = aux8_0.re;

		cmpx aux9_0 = multiplyComplex(hp2x->GetBinContent(i), hp2y->GetBinContent(i), q2x, -1*q2y);
		float aux9  = 2*aux9_0.re;

		float aux10 = 2*m[i]*n;
		float aux11 = 6*m[i];
		float aux12 = (m[i]*n-3*m[i])*(n-1)*(n-2);

		hsb_4prime->SetBinContent(i, (aux1-aux2-aux3-aux4-aux5+aux6-aux7+aux8+aux9+aux10-aux11)/(aux12));
	}

	for(int i=1; i<=hpT_Template->GetNbinsX(); i++)
	{
		hdb_2prime->Fill(hsb_2prime->GetBinCenter(i), hsb_2prime->GetBinContent(i),1);
		hdb_4prime->Fill(hsb_4prime->GetBinCenter(i), hsb_4prime->GetBinContent(i),1);
	}
}

void processEvent()
{
	if(collisionparticles.size() == 0) return;

	float cmx=0;
	float cmy=0;
	float xsum=0;
	float ysum=0;
	float avercos2=0;
	float avercos3=0;
	float aversin2=0;
	float aversin3=0;
	float aver2=0;

	//Calculate centroid
	for (unsigned int i=0; i<collisionparticles.size(); i++)
	{
		xsum = xsum + collisionparticles[i].x;
		ysum = ysum + collisionparticles[i].y;
	}

	cmx = xsum/collisionparticles.size();
	cmy = ysum/collisionparticles.size();

	//Calculate each useful value of collision particles
	//Store them in vectors
	for (unsigned int i=0; i<collisionparticles.size(); i++)
	{
		//Shift to center of mass frame
		collisionparticles[i].x = collisionparticles[i].x - cmx;
		collisionparticles[i].y = collisionparticles[i].y - cmy;

		collisionparticles[i].phi = TMath::ATan2(collisionparticles[i].y,collisionparticles[i].x);
		collisionparticles[i].rsquare = collisionparticles[i].x*collisionparticles[i].x + collisionparticles[i].y*collisionparticles[i].y;
		collisionparticles[i].pT = TMath::Sqrt(collisionparticles[i].px*collisionparticles[i].px + collisionparticles[i].py*collisionparticles[i].py);
	}

	//Calculate the average values for computing epsilon_2
	for (unsigned int i=0; i<collisionparticles.size(); i++)
	{
		avercos2 = avercos2 + collisionparticles[i].rsquare * TMath::Cos(2*collisionparticles[i].phi);
		aversin2 = aversin2 + collisionparticles[i].rsquare * TMath::Sin(2*collisionparticles[i].phi);
		aver2    = aver2 + collisionparticles[i].rsquare;
	}

	avercos2 = avercos2 / collisionparticles.size();
	aversin2 = aversin2 / collisionparticles.size();
	aver2    = aver2 / collisionparticles.size();

	//Calculate the average values for n=3
	for (unsigned int i=0; i<collisionparticles.size(); i++)
	{
		avercos3 = avercos3 + collisionparticles[i].rsquare * TMath::Cos(3*collisionparticles[i].phi);
		aversin3 = aversin3 + collisionparticles[i].rsquare * TMath::Sin(3*collisionparticles[i].phi);
	}

	avercos3 = avercos3 / collisionparticles.size();
	aversin3 = aversin3 / collisionparticles.size();

	//Calculate epsilon_n and psi_n and put them into vectors
	//Using formulas from arXiv:1501.06880
	float e2;
	float e3;
	float s2;
	float s3;

	e2 = TMath::Sqrt(avercos2*avercos2 + aversin2*aversin2) / aver2;
	e3 = TMath::Sqrt(avercos3*avercos3 + aversin3*aversin3) / aver2;
	s2 = (TMath::ATan2(aversin2,avercos2) + TMath::Pi())/2.0;
	s3 = (TMath::ATan2(aversin3,avercos3) + TMath::Pi())/3.0;

	epsilon2.push_back(e2);
	epsilon3.push_back(e3);
	psi2.push_back(s2);
	psi3.push_back(s3);

	hPsi2->Fill(s2);

	//Calculate v2 and v3
	for(unsigned int i=0; i<finalparticles.size(); i++)
	{
		//Only particles at midrapidity
		if(finalparticles[i].eta < -1 || finalparticles[i].eta > 1) continue;

		float v2_val = TMath::Cos(2*(finalparticles[i].phi - s2));
		float v3_val = TMath::Cos(3*(finalparticles[i].phi - s3));

		hv2_pT->Fill(finalparticles[i].pT,v2_val,1);
		hv3_pT->Fill(finalparticles[i].pT,v3_val,1);
	}

}

void parseGeometry(int proc)
{
	process_number = proc;

	//Initialize histogram
	hPsi2      = new TH1F("hPsi2","hPsi2",50,0,TMath::Pi());
	hv2_pT     = new TProfile("hv2_pT","hv2_pT",NOB,0,5,-2,2);
	hv3_pT     = new TProfile("hv3_pT","hv3_pT",NOB,0,5,-2,2);
	hpT_Template = new TH1F("hpT_Template","hpT_Template",NOB,0,5);
	hdb_2prime = new TProfile("hdb_2prime","hdb_2prime",NOB,0,5,-500,500);
	hdb_4prime = new TProfile("hdb_4prime","hdb_4prime",NOB,0,5,-500,500);
	hd_2       = new TH1F("hd_2","hd_2",NOB,0,5);
	hdb_2      = new TProfile("hdb_2","hdb_2",1,0,1,-100,100);
	hdb_4      = new TProfile("hdb_4","hdb_4",1,0,1,-100,100);

	hEta         = new TH1F("hEta","hEta;#eta;Counts",NOB,-6,6);
  	hpT          = new TH1F("hpT","hpT;p_{T};Counts",NOB,0,10);
  	hpT_EtaCut   = new TH1F("hpT_EtaCut","hpT_EtaCut;p_{T};Counts",NOB,0,20);
  	hPhi         = new TH1F("hPhi","hPhi;#phi;Counts",NOB,-TMath::Pi(),TMath::Pi());

	//Read in test.f20 file
	ifstream dataFile;
	dataFile.open("/direct/phenix+hhj/jdok/UrQMD/urqmd-3.4/test.f20");
	if (!dataFile)
	{
		printf("File does not exist\n");
		return;
	}
	else
	{
		cout << Form("--> Successfully opened file") << endl << endl;
	}

	string linestr;
	vector<string> tokens;
	while(dataFile)
	{
		tokens.clear();
		std::getline(dataFile,linestr);
		istringstream iss(linestr);
		copy(istream_iterator<string>(iss), istream_iterator<string>(), back_inserter(tokens));
		if(parse_verbosity) cout << "*     " << linestr << endl;

  		//Ignore file header starting with #
		if(tokens[0] == "#") continue;

  		//Find new event header
		if(atoi(tokens[0].c_str()) == 0 && atoi(tokens[1].c_str()) == NUCL)
		{
			numevent++;

    		//Skip the following lines
			for(int i=0; i<NUCL; i++)
			{
				std::getline(dataFile,linestr);
				if(parse_verbosity) cout << "**    " << linestr << endl;
			}

			continue;
		}

 		//Find collision header in the current event for 2 -> many scatterings
		if(tokens[0] == "2")
		{
			int numcollproducts = atoi(tokens[1].c_str());

			particle p1;
			particle p2;

			tokens.clear();
			std::getline(dataFile,linestr);
			istringstream iss2(linestr);
			copy(istream_iterator<string>(iss2), istream_iterator<string>(), back_inserter(tokens));
			if(parse_verbosity) cout << "+     " << linestr << endl;


    		//Is this scattering between nucleons?
			if(atoi(tokens[0].c_str()) >= 1 && atoi(tokens[0].c_str()) <= NUCL)
			{
  				//Store px,py,pz,x,y,z for first particle
				p1.px=atof(tokens[3].c_str());
				p1.py=atof(tokens[4].c_str());
				p1.pz=atof(tokens[5].c_str());
				p1.x=atof(tokens[8].c_str());
				p1.y=atof(tokens[9].c_str());
				p1.z=atof(tokens[10].c_str());

				tokens.clear();
				std::getline(dataFile,linestr);
				istringstream iss3(linestr);
				copy(istream_iterator<string>(iss3), istream_iterator<string>(), back_inserter(tokens));
				if(parse_verbosity) cout << "+     " << linestr << endl;

				if(atoi(tokens[0].c_str()) >= 1 && atoi(tokens[0].c_str()) <= NUCL)
				{
    				//Store px,py,pz,x,y,z for second particle
					p2.px=atof(tokens[3].c_str());
					p2.py=atof(tokens[4].c_str());
					p2.pz=atof(tokens[5].c_str());
					p2.x=atof(tokens[8].c_str());
					p2.y=atof(tokens[9].c_str());
					p2.z=atof(tokens[10].c_str());

					collisionparticles.push_back(p1);
					collisionparticles.push_back(p2);

					npartsum = npartsum + 2;
					npart = npart + 2;
				}

				for(int i=0; i<numcollproducts; i++)
				{
					std::getline(dataFile,linestr);
					if(parse_verbosity) cout << "***   " << linestr << endl;
				}
				continue;
			}

			//Skip over the rest of the lines corresponding to collision products
			for(int i=0; i<numcollproducts+1; i++)
			{
				std::getline(dataFile,linestr);	
				if(parse_verbosity) cout << "***   " << linestr << endl;
			}
			continue;
		}

		//Identify final-state particles and store them
		if(atoi(tokens[0].c_str()) > 0 && tokens[1] == "0" && tokens.size() == 2)
		{
			if(parse_verbosity) cout << "END OF EVENT" << endl;
			for(int i=0; i<atoi(tokens[0].c_str()); i++)
			{
				vector<string> finalstate_tokens;
				std::getline(dataFile,linestr);
				istringstream iss4(linestr);
				copy(istream_iterator<string>(iss4), istream_iterator<string>(), back_inserter(finalstate_tokens));				

				particle pf;
				TLorentzVector ev(atof(finalstate_tokens[3].c_str()),atof(finalstate_tokens[4].c_str()),atof(finalstate_tokens[5].c_str()),atof(finalstate_tokens[6].c_str()));
				pf.pT  = ev.Pt();
				pf.phi = ev.Phi();
				pf.eta = ev.Eta();
				pf.px = atof(finalstate_tokens[3].c_str());				
				pf.py = atof(finalstate_tokens[4].c_str());
				pf.pz = atof(finalstate_tokens[5].c_str());

				if(atoi(finalstate_tokens[0].c_str()) <= NUCL) nspectator++;

				finalparticles.push_back(pf);	
				if(parse_verbosity) cout << "****  " << linestr << endl;
			}
		}

  		//Process event
		if(tokens[0] == "0" && tokens[1] == "0")
		{
			//processEvent();
			processEventCumulants();
			collisionparticles.clear();
			finalparticles.clear();
			npart = 0;
		}

		if (!dataFile) break;

	}

	computeCumulants();
	draw();

	//Compute mean epsilon2 and epsilon3
	float ep2avg = 0;
	float ep3avg = 0;

	for(unsigned int i=0; i<epsilon2.size(); i++)
	{
		ep2avg += epsilon2[i];
	}
	ep2avg = ep2avg/epsilon2.size();

	for(unsigned int i=0; i<epsilon3.size(); i++)
	{
		ep3avg += epsilon3[i];
	}
	ep3avg = ep3avg/epsilon3.size();

	cout << "Nevt         = " << numevent << endl;
	cout << "<ep2>        = " << ep2avg << endl;
	cout << "<ep3>        = " << ep3avg << endl;
	cout << "<Npart_NN>   = " << (float) npartsum/numevent << endl;
	cout << "<Nspec>      = " << (float) nspectator/numevent << endl;
	cout << "<Npart>      = " << NUCL - (float) nspectator/numevent << endl;

	//writeData();
}
