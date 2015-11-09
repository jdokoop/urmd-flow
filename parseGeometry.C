
//-----------------------------------------------
// Code to parse the OSCAR199A output file from
// Ultrarelativistic Quantum Molecular Dynamics
// and calculate anisotropic flow coefficients
//   (1) with respect to the event plane
//   (2) with two-particle cumulants
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

//Number of bins in histograms
const int NOB = 100;

int npart = 0;
int npartsum = 0;
int numevent=0;
int process_number = 0;

//Vectors to store event-wise geometric parameters
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

//Miscellaneous options
bool parse_verbosity = false;

//---------------------------------------------------------
// Functions
//---------------------------------------------------------

void writeData()
{
//Write event plane angles to file
/*	
ofstream myfile;
myfile.open ("psi2.txt");
for(unsigned int i=0; i<psi2.size(); i++)
{
	myfile << psi2[i] << endl;
}
myfile.close();

myfile.open ("psi3.txt");
for(unsigned int i=0; i<psi3.size(); i++)
{
	myfile << psi3[i] << endl;
}
myfile.close();
*/

TFile *fout = new TFile(Form("/direct/phenix+hhj/jdok/UrQMD_Output/urqmd_out_%i.root",process_number),"RECREATE");
hv2_pT->Write();
hv3_pT->Write();
fout->Close();
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
	hPsi2  = new TH1F("hPsi2","hPsi2",50,0,TMath::Pi());
	hv2_pT = new TProfile("hv2_pT","hv2_pT",NOB,0,5,-2,2);
	hv3_pT = new TProfile("hv3_pT","hv3_pT",NOB,0,5,-2,2);

	//Read in test.f20 file
	ifstream dataFile;
	dataFile.open("/direct/phenix+hhj/jdok/UrQMD/urqmd-3.4/test_dPb_5TeV.f20");
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
		if(tokens[0] == "0" && tokens[1] == "210")
		{
			numevent++;

    		//Skip the following lines
			for(int i=0; i<210; i++)
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
			if(atoi(tokens[0].c_str()) >= 1 && atoi(tokens[0].c_str()) <= 210)
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

				if(atoi(tokens[0].c_str()) >= 1 && atoi(tokens[0].c_str()) <= 210)
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

				finalparticles.push_back(pf);	
				if(parse_verbosity) cout << "****  " << linestr << endl;
			}
		}

  		//Process event
		if(tokens[0] == "0" && tokens[1] == "0")
		{
			processEvent();
			collisionparticles.clear();
			finalparticles.clear();
			npart = 0;
		}

		if (!dataFile) break;

	}

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

	cout << "Nevt    = " << numevent << endl;
	cout << "<ep2>   = " << ep2avg << endl;
	cout << "<ep3>   = " << ep3avg << endl;
	cout << "<Npart> = " << (float) npartsum/numevent << endl;

	//writeData();
}
