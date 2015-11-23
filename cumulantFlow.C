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
#include "TF1.h"
#include "TProfile.h"

using namespace std;

//-------------------------------------------
// Variables
//-------------------------------------------

//Data structure to contain particles
struct particle
{
	int id;
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

//Variables for reference flow calculation
TProfile *h_db_2;     // <<2>>
TProfile *h_db_4;     // <<4>>

TProfile *h_db_2prime; // <<2'>>
TProfile *h_db_4prime; // <<4'>>

//Participant particles in each event
vector<particle> collisionparticles;

//Final state hadrons in each event
vector<particle> finalparticles;

//Spectators determined by looking at final state PIDs
vector<particle> spectatorParticles;

//Participants determined by looking at final state PIDS
vector<particle> participantParticles;

bool parse_verbosity = false;

//-------------------------------------------
// Functions
//-------------------------------------------

void computeFlow()
{
	//REFERENCE FLOW

	//Start by computing q-vector
	float q2x = 0;
	float q2y = 0;

	float q4x = 0;
	float q4y = 0;

	int M = finalparticles.size();

	for(int i=0; i<M; i++)
	{
		particle p = finalparticles[i];
		float phi = p.phi;

		q2x += TMath::Cos(2*phi);
		q2y += TMath::Sin(2*phi);

		q4x += TMath::Cos(4*phi);
		q4y += TMath::Sin(4*phi);	
	}

	//Single event reference cumulants
	float sb_2 = 0;
	float sb_4 = 0;

	sb_2 = ((q2x*q2x + q2y*q2y) - M)/(M*(M-1));

	float aux1 = q2x*q2x + q2y*q2y;
	float aux2 = q4x*q4x + q4y*q4y;
	float aux3 = -2*(q2x*q2x*q4x - q2y*q2y*q4x + 2*q2x*q2y*q4y);
	float aux4 = -4*(M-2)*(q2x*q2x + q2y*q2y);
	float aux5 = 2*M*(M-3);
	sb_4 = (aux1+aux2+aux3+aux4+aux5)/(M*(M-1)*(M-2)*(M-3));

	//Fill TProfile to compute event averages
	h_db_2->Fill(0.5,sb_2,1);
	h_db_4->Fill(0.5,sb_4,1);

	//DIFFERENTIAL FLOW

	//Compute p-vector for particles of interest (each pT bin)
	TH1F *hp2x = new TH1F("hp2","hp2",50,0,5);
	TH1F *hp4x = new TH1F("hp4","hp4",50,0,5);
	TH1F *hp2y = new TH1F("hp2","hp2",50,0,5);
	TH1F *hp4y = new TH1F("hp4","hp4",50,0,5);
	TH1F *hm  = new TH1F("hm","hm",50,0,5);

	//Initialize p-vectors with zeroes
	for(int i=0; i<hp2x->GetNbinsX(); i++)
	{
		hp2x->SetBinContent(i+1,0);
		hp4x->SetBinContent(i+1,0);
		hp2y->SetBinContent(i+1,0);
		hp4y->SetBinContent(i+1,0);
	}

	for(int i=0; i<M; i++)
	{
		particle p = finalparticles[i];
		float pT = sqrt(p.px*p.px + p.py*p.py);
		float phi = p.phi;
		int bin = hm->FindBin(pT);

		//Keep track of how many particles there are per bin
		hm->Fill(pT);

		hp2x->SetBinContent(bin,hp2x->GetBinContent(i) + TMath::Cos(2*phi));
		hp2y->SetBinContent(bin,hp2y->GetBinContent(i) + TMath::Sin(2*phi));

		hp4x->SetBinContent(bin,hp4x->GetBinContent(i) + TMath::Cos(4*phi));
		hp4y->SetBinContent(bin,hp4y->GetBinContent(i) + TMath::Sin(4*phi));
	}

	//Single event differential cumulants
	TH1F *hsb_2prime;
	hsb_2prime = new TH1F("hsb_2prime","hsb_2prime",50,0,5);

	TH1F *hsb_4prime;
	hsb_4prime = new TH1F("hsb_4prime","hsb_4prime",50,0,5);

	for(int i=0; i<hsb_2prime->GetNbinsX(); i++)
	{
		int n = hm->GetBinContent(i+1);
		float val = ((hp2x->GetBinContent(i+1))*(q2x) + (hp2y->GetBinContent(i+1))*(q2y))/(n*M-n);
		hsb_2prime->SetBinContent(i+1,val);
	}

	//Accumulate event averages
	for(int i=0; i<hsb_2prime->GetNbinsX(); i++)
	{
		h_db_2prime->Fill(hsb_2prime->GetBinCenter(i+1),hsb_2prime->GetBinContent(i+1),1);
	}
}

void processEvent()
{
	computeFlow();
}

void parseFile20()
{
	//Read in test.f20 file
	ifstream dataFile;
	//dataFile.open("/direct/phenix+hhj/jdok/UrQMD/urqmd-3.4/test.f20");
	dataFile.open("/direct/phenix+hhj2/jdok/urqmd-hulthen-3.4/test.f20");
	if (!dataFile)
	{
		printf("File does not exist\n");
		return;
	}
	else
	{
		cout << Form("--> Successfully opened file f20!") << endl << endl;
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
				pf.id = atoi(finalstate_tokens[0].c_str());

				if(atoi(finalstate_tokens[0].c_str()) <= NUCL)
				{
					nspectator++;
					spectatorParticles.push_back(pf);
				}

				finalparticles.push_back(pf);
				if(parse_verbosity) cout << "****  " << linestr << endl;
			}
		}

  		//Process event
		if(tokens[0] == "0" && tokens[1] == "0")
		{
			processEvent();
			finalparticles.clear();
			spectatorParticles.clear();
			collisionparticles.clear();
			participantParticles.clear();
			npart = 0;
		}

		if (!dataFile) break;
	}
}

void cumulantFlow(int proc)
{
	//Set label for output file (when running in parallel)
	process_number = proc;

	//Initialize histograms
	h_db_2 = new TProfile("h_db_2","h_db_2",1,0,1,-1000,1000);
	h_db_4 = new TProfile("h_db_4","h_db_4",1,0,1,-1000,1000);

	h_db_2prime = new TProfile("h_db_2prime","h_db_2prime",50,0,5,-1000,1000);
	h_db_4prime = new TProfile("h_db_4prime","h_db_4prime",50,0,5,-1000,1000);

	parseFile20();
	h_db_2prime->Scale(1.0/sqrt((h_db_2->GetBinContent(1))));
	h_db_2prime->Draw();
}