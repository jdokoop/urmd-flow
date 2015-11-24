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

const int NBIN = 50;

//Event characterization variables
int npart          = 0;
int npartsum       = 0;
int nspectator     = 0;
int numevent       = 0;
int process_number = 0;

//Event-wise reference flow
//Must be cleared after every event
float q2x  = 0;
float q2y  = 0;
float q4x  = 0;
float q4y  = 0;
float sb_2 = 0;
float sb_4 = 0;

TH1F *h_sb_2prime;
TH1F *h_sb_4prime;

TH1F *hp2x;
TH1F *hp2y;
TH1F *hp4x;
TH1F *hp4y;

TH1F *hn;

//Average over events
TProfile *h_db_2;
TProfile *h_db_4;
TProfile *h_db_2prime;
TProfile *h_db_4prime;

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
	float c2 = h_db_2->GetBinContent(1);
	h_db_2prime->Scale(1.0/sqrt(c2));
	h_db_2prime->Draw();
}

void processFlow()
{
	//Clear event-wise variables
	q2x  = 0;
	q2y  = 0;
	q4x  = 0;
	q4y  = 0;
	sb_2 = 0;
	sb_4 = 0;
	h_sb_2prime->Reset();
	h_sb_4prime->Reset();
	hp2x->Reset();
	hp2y->Reset();
	hp4x->Reset();
	hp4y->Reset();
	hn->Reset();

	//Loop over particles, computing Q-vector
	int M = finalparticles.size();
	for(int i=0; i<M; i++)
	{
		q2x = q2x + TMath::Cos(2*finalparticles[i].phi);
		q2y = q2y + TMath::Sin(2*finalparticles[i].phi);

		q4x = q4x + TMath::Cos(4*finalparticles[i].phi);
		q4y = q4y + TMath::Sin(4*finalparticles[i].phi);		
	}

	//Reference flow
	sb_2 = ((q2x*q2x + q2y*q2y) - M)/(M*(M - 1));
	h_db_2->Fill(0.5,sb_2,1);

	float sb_4_aux1 = (q2x*q2x + q2y*q2y)*(q2x*q2x + q2y*q2y);
	float sb_4_aux2 = q4x*q4x + q4y*q4y;
	float sb_4_aux3 = -2*(q2x*q2x*q4x - q2y*q2y*q4x + 2*q2x*q2y*q4y);
	float sb_4_aux4 = -4*(M-2)*(q2x*q2x + q2y*q2y);
	float sb_4_aux5 = 2*M*(M-3);
	float sb_4_aux6 = M*(M-1)*(M-2)*(M-3);
	sb_4 = (sb_4_aux1 + sb_4_aux2 + sb_4_aux3 + sb_4_aux4 + sb_4_aux5)/sb_4_aux6;
	h_db_4->Fill(0.5,sb_4,1);

	//Differential cumulants
	for(int i=0; i<M; i++)
	{
		float pT = sqrt(pow(finalparticles[i].px,2) + pow(finalparticles[i].py,2));
		int bin = hp2x->FindBin(pT);

		hp2x->SetBinContent(bin, hp2x->GetBinContent(bin) + TMath::Cos(2*finalparticles[i].phi));
		hp2y->SetBinContent(bin, hp2y->GetBinContent(bin) + TMath::Sin(2*finalparticles[i].phi));
		hp4x->SetBinContent(bin, hp4x->GetBinContent(bin) + TMath::Cos(4*finalparticles[i].phi));
		hp4y->SetBinContent(bin, hp4y->GetBinContent(bin) + TMath::Sin(4*finalparticles[i].phi));

		hn->Fill(pT);
	}

	//Computing <2'>
	for(int i=1; i<NBIN; i++)
	{
		float p2x = hp2x->GetBinContent(i);
		float p2y = hp2y->GetBinContent(i);
		float p4x = hp4x->GetBinContent(i);
		float p4y = hp4y->GetBinContent(i);

		float sb_2prime_aux = p2x*q2x + p2y*q2y;
		int n = hn->GetBinContent(i);

		float sb_2prime_cont = (float) (sb_2prime_aux - n)/(n*M - n);
		h_sb_2prime->SetBinContent(i,sb_2prime_cont);
	}

	//Accumulate in <<2'>>
	for(int i=1; i<NBIN; i++)
	{
		float pT = h_sb_2prime->GetBinCenter(i);
		h_db_2prime->Fill(pT, h_sb_2prime->GetBinContent(i), 1);
	}
}

void processEvent()
{
	processFlow();
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
	h_db_2      = new TProfile("h_db_2","h_db_2",1,0,1,-500,500);
	h_db_4      = new TProfile("h_db_4","h_db_4",1,0,1,-500,500);
	h_db_2prime = new TProfile("h_db_2prime","h_db_2prime",NBIN,0,5,-500,500);
	h_db_4prime = new TProfile("h_db_4prime","h_db_4prime",NBIN,0,5,-500,500);

	h_sb_2prime = new TH1F("h_sb_2prime","h_sb_2prime",NBIN,0,5);
	h_sb_4prime = new TH1F("h_sb_4prime","h_sb_4prime",NBIN,0,5);
	hn          = new TH1F("hn","hn",NBIN,0,5);
	hp2x        = new TH1F("hp2x","hp2x",NBIN,0,5);
	hp2y        = new TH1F("hp2y","hp2y",NBIN,0,5);
	hp4x        = new TH1F("hp4x","hp4x",NBIN,0,5);
	hp4y        = new TH1F("hp4y","hp4y",NBIN,0,5);

	parseFile20();

	computeFlow();
}