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

//Number of bins in histograms
const int NOB = 100;

//Number of nucleons in system
// --> p+Au = 198
// --> d+Au = 199
// --> d+Pb = 210
// --> p+Pb = 209
const int NUCL = 198;

//Event characterization variables
int npart          = 0;
int npartsum       = 0;
int nspectator     = 0;
int numevent       = 0;
int process_number = 0;

//Vector to store initial nucleon coordinates for every event
vector< vector<float> > r0_x;
vector< vector<float> > r0_y;

//Vectors to store event-level geometric parameters
vector<float> epsilon2;
vector<float> epsilon3;
vector<float> psi2;
vector<float> psi3;

//Participant particles in each event
vector<particle> collisionparticles;

//Final state hadrons in each event
vector<particle> finalparticles;

//Spectators determined by looking at final state PIDs
vector<particle> spectatorParticles;

//Participants determined by looking at final state PIDS
vector<particle> participantParticles;

//Histograms for plotting
TH1F *hPsi2;
TProfile *hv2_pT;
TProfile *hv3_pT;

//Miscellaneous options
bool parse_verbosity = false;

//Smear nucleon positions with Gaussian?
bool smear_nucleons = false;

//---------------------------------------------------------
// Functions
//---------------------------------------------------------

void writeData()
{	
	TFile *fout = new TFile(Form("/direct/phenix+hhj/jdok/UrQMD_Output/urqmd_out_%i.root",process_number),"RECREATE");
	
	fout->Close();
}

bool isParticipant(int identifier)
{
	for(int i=0; i<spectatorParticles.size(); i++)
	{
		int id = spectatorParticles[i].id;

		if(identifier == id) return false;
	}

	return true;
}

void determineParticipants()
{
	//Get the original nucleon positions for the event at hand
	int eventNucleons = r0_x[numevent-1].size();
	vector<float> nucleons_x = r0_x[numevent-1];
	vector<float> nucleons_y = r0_y[numevent-1];

	for(int i=0; i< eventNucleons; i++)
	{
		if(!isParticipant(i+1)) continue;

		particle p;
		p.x = nucleons_x[i];
		p.y = nucleons_y[i];

		participantParticles.push_back(p);
	}
}

void processEvent()
{
	if(participantParticles.size() == 0) return;

	double Count = participantParticles.size();

	float cmx=0;
	float cmy=0;
	float xsum=0;
	float ysum=0;
	float avercos2=0;
	float avercos3=0;
	float aversin2=0;
	float aversin3=0;
	float aver2=0;

	//If Gaussian smearing is enabled for nucleon positions
	double SmearSigma = 0.4; // units of fm
	TF1 *fradius = new TF1("fradius","x*TMath::Exp(-x*x/(2*[0]*[0]))",0.0,2.0);
	fradius->SetParameter(0,SmearSigma);
	TF1 *fphi = new TF1("fphi","1.0",0.0,2.0*TMath::Pi());
	int SmearSamplings = 100;

	//Calculate centroid
	for (unsigned int i=0; i<participantParticles.size(); i++)
	{
		xsum = xsum + participantParticles[i].x;
		ysum = ysum + participantParticles[i].y;
	}

	cmx = xsum/participantParticles.size();
	cmy = ysum/participantParticles.size();

	//Calculate each useful value of collision particles
	//Store them in vectors
	for (unsigned int i=0; i<participantParticles.size(); i++)
	{
		//Shift to center of mass frame
		participantParticles[i].x = participantParticles[i].x - cmx;
		participantParticles[i].y = participantParticles[i].y - cmy;

		participantParticles[i].phi = TMath::ATan2(participantParticles[i].y,participantParticles[i].x);
		participantParticles[i].rsquare = participantParticles[i].x*participantParticles[i].x + participantParticles[i].y*participantParticles[i].y;
		participantParticles[i].pT = TMath::Sqrt(participantParticles[i].px*participantParticles[i].px + participantParticles[i].py*participantParticles[i].py);
	}

	if(smear_nucleons)
	{
		for(int i=0; i<participantParticles.size(); i++)
		{
			for (int is=0;is<SmearSamplings;is++) 
			{
				double rtemp   = fradius->GetRandom();
				double phitemp = fphi->GetRandom();
				double xtemp   = participantParticles[i].x + rtemp*TMath::Sin(phitemp);
				double ytemp   = participantParticles[i].y + rtemp*TMath::Cos(phitemp);

				float r    = TMath::Sqrt(xtemp*xtemp + ytemp*ytemp);
				float phi  = TMath::ATan2(ytemp,xtemp);

				avercos2 += r*r*TMath::Cos(2*phi);
				aversin2 += r*r*TMath::Sin(2*phi);
				avercos3 += r*r*TMath::Cos(3*phi);
				aversin3 += r*r*TMath::Sin(3*phi);
				aver2 += r*r;

			Count = Count + 1;  // noting this is a double
		}
	}
}
else
	{
		//Calculate the average values for computing epsilon_2
		for (unsigned int i=0; i<participantParticles.size(); i++)
		{
			avercos2 = avercos2 + participantParticles[i].rsquare * TMath::Cos(2*participantParticles[i].phi);
			aversin2 = aversin2 + participantParticles[i].rsquare * TMath::Sin(2*participantParticles[i].phi);
			aver2    = aver2 + participantParticles[i].rsquare;
		}

		//Calculate the average values for n=3
		for (unsigned int i=0; i<participantParticles.size(); i++)
		{
			avercos3 = avercos3 + participantParticles[i].rsquare * TMath::Cos(3*participantParticles[i].phi);
			aversin3 = aversin3 + participantParticles[i].rsquare * TMath::Sin(3*participantParticles[i].phi);
		}
	}

	avercos2 = avercos2 / Count;
	aversin2 = aversin2 / Count;

	avercos3 = avercos3 / Count;
	aversin3 = aversin3 / Count;

	aver2    = aver2 / Count;

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

void parseFile14()
{	
	ifstream dataFile14;
	dataFile14.open("/direct/phenix+hhj/jdok/UrQMD/urqmd-3.4/test.f14");
	if (!dataFile14)
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
	vector<float> xvals;
	vector<float> yvals;
	int evtHeaderCounter = 0;
	while(dataFile14)
	{	
		//Skip the 17 event header lines
		if(evtHeaderCounter < 17)
		{
			std::getline(dataFile14,linestr);
			evtHeaderCounter++;
			continue;
		}

		//Skip next 2 lines
		tokens.clear();
		std::getline(dataFile14,linestr);
		istringstream iss1(linestr);
		copy(istream_iterator<string>(iss1), istream_iterator<string>(), back_inserter(tokens));

		//Determine total number of nucleons and do a sanity check
		int nnucleons = atoi(tokens[0].c_str());
		if(nnucleons != NUCL) cout << "Error! Wrong number of nucleons." << endl;

		//Skip next line
		tokens.clear();
		std::getline(dataFile14,linestr);

		//Iterate over nucleons
		for(int i=0; i<nnucleons; i++)
		{
			tokens.clear();
			std::getline(dataFile14,linestr);
			istringstream iss2(linestr);
			copy(istream_iterator<string>(iss2), istream_iterator<string>(), back_inserter(tokens));

			xvals.push_back(atof(tokens[1].c_str()));
			yvals.push_back(atof(tokens[2].c_str()));
		}

		//Finish processing event
		r0_x.push_back(xvals);
		r0_y.push_back(yvals);
		xvals.clear();
		yvals.clear();

		//Skip over remainder of event
		tokens.clear();
		std::getline(dataFile14,linestr);
		istringstream iss3(linestr);
		copy(istream_iterator<string>(iss3), istream_iterator<string>(), back_inserter(tokens));

		int nskip = atoi(tokens[0].c_str()) + 1;
		for(int i=0; i<nskip; i++)
		{		
			std::getline(dataFile14,linestr);
		}

		evtHeaderCounter = 0;
	}

	if(parse_verbosity)
	{
		for(int i=0; i<r0_x.size(); i++)
		{
			cout << "********** EVENT " << i << "*************" << endl;

			vector<float> vx = r0_x[i];
			vector<float> vy = r0_y[i];
			for(int j=0; j<vx.size(); j++)
			{
				cout << "   " << vx[j] << "  " << vy[j] << endl;
			}
			cout << endl;
		}
	}
}

void parseFile20()
{
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
			determineParticipants();
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

void parseEventPlane(int proc)
{
	//Set label for output file (when running in parallel)
	process_number = proc;

	//Initialize histograms
	hPsi2  = new TH1F("hPsi2","hPsi2",50,0,TMath::Pi());
	hv2_pT = new TProfile("hv2_pT","hv2_pT",NOB,0,5,-2,2);
	hv3_pT = new TProfile("hv3_pT","hv3_pT",NOB,0,5,-2,2);

	//Parse geometry files
	parseFile14();
	parseFile20();

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

/*
	writeData();
	*/
}
