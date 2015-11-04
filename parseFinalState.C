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
	float phi;
	float eta;
	float pT;
};

int nevent = 0;

vector<float> psi_2;
vector<float> psi_3;
vector<particle> particles;

TProfile *hv2_pT = new TProfile("hv2_pT","hv2_pT",50,0,5,-2,2);
TProfile *hv3_pT = new TProfile("hv3_pT","hv3_pT",50,0,5,-2,2);

//-------------------------------------------
// Functions
//-------------------------------------------

void writeFile(int process)
{
	TFile *fout = new TFile(Form("/direct/phenix+hhj/jdok/UrQMD_Output/urqmd_out_%i.root",process),"RECREATE");
	hv2_pT->Write();
	hv3_pT->Write();
	fout->Close();
}

void loopParticles()
{
	if(nevent % 100 == 0)
	{
		cout << "--> Processing event No. " << nevent << endl;
	}

	for(int i=0; i<particles.size(); i++)
	{
		if(particles[i].eta > -1 && particles[i].eta < 1) continue;

		float v2 = TMath::Cos(2*(particles[i].phi - psi_2[nevent]));
		hv2_pT->Fill(particles[i].pT, v2, 1);

		float v3 = TMath::Cos(3*(particles[i].phi - psi_2[nevent]));
		hv3_pT->Fill(particles[i].pT, v3, 1);		
	}
}

void loadEventPlane()
{
	//Psi2
	ifstream psi_2_file;
	psi_2_file.open("psi2.txt");

	float psi_2_value = 0;
	while(psi_2_file)
	{
		psi_2_file >> psi_2_value;
		psi_2.push_back(psi_2_value);
	}

	psi_2_file.close();

	//Psi3
	ifstream psi_3_file;
	psi_3_file.open("psi3.txt");

	float psi_3_value = 0;
	while(psi_3_file)
	{
		psi_3_file >> psi_3_value;
		psi_3.push_back(psi_3_value);
	}

	psi_3_file.close();
}

void parseFinalState(int process)
{
	//Load event plane from file
	loadEventPlane();

	cout << psi_2.size() << endl;
	//Open File
	ifstream myFile;
	myFile.open("test.f14");

	if (!myFile) 
	{
		printf("Input file does not exist!\n");
		return;
	}
	else
	{
		cout << "--> Successfully opened file! " << endl;
	}

	string linestr;
	vector<string> tokens;
	while (myFile) 
	{
		//Skip the 17 event header lines
		for(int i=0; i<17; i++)
		{
			std::getline(myFile,linestr);
		}

		//Get the number of particles in the event
		tokens.clear();
		std::getline(myFile,linestr);
		istringstream iss(linestr);
		copy(istream_iterator<string>(iss), istream_iterator<string>(), back_inserter(tokens));
		int npart = atoi(tokens[0].c_str());

		//Skip over collision counter line
		std::getline(myFile,linestr);

		//Loop over the npart particles in the event
		for(int i=0; i<npart; i++)
		{
			tokens.clear();
			std::getline(myFile,linestr);
			istringstream iss(linestr);
			copy(istream_iterator<string>(iss), istream_iterator<string>(), back_inserter(tokens));

			float energy = atof(tokens[4].c_str());
			float px     = atof(tokens[5].c_str());	
			float py     = atof(tokens[6].c_str());
			float pz     = atof(tokens[7].c_str());
			int pid      = atoi(tokens[9].c_str());
			int chg      = atoi(tokens[11].c_str());

			TLorentzVector ev(px,py,pz,energy);

			float pT  = ev.Pt();
			float phi = ev.Phi();
			float eta = ev.Eta();

			particle p;
			p.pT  = pT;
			p.phi = phi;
			p.eta = eta;

			if(chg != 0)
			{
				particles.push_back(p);
			}
		}

		//Process event
		loopParticles();
		nevent++;
		particles.clear();
	}  

	writeFile(process);
}
