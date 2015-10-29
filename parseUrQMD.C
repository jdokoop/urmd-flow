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
  int   count;
  int   partid;
  float x;
  float y;
  float z;
};

int npart;
int numevent=0;

particle p;

//-------------------------------------------
// Functions
//-------------------------------------------

void parseUrQMD()
{
  //Read in test.f20 file
  ifstream dataFile;
  dataFile.open("test200mb.f20");
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

      //Ignore file header starting with #
      if(tokens[0] == "#") continue;

      //Find new event header
      if(tokens[0] == "0" && tokens[1] == "199")
      {
        cout << "NEW EVENT!" << endl;
        numevent++;
        //Skip the following 199 lines
        for(int i=0; i<199; i++)
        {
          std::getline(dataFile,linestr);
        }
      }

      //Find collision header in the current event for 2 -> many scatterings
      if(tokens[0] == "2")
      {
        cout << "  NEW COLLISION" << endl;
        tokens.clear();
        std::getline(dataFile,linestr);
        istringstream iss(linestr);
        copy(istream_iterator<string>(iss), istream_iterator<string>(), back_inserter(tokens));

        //Is this scattering between nucleons?
        if(std::stoi(tokens[0]) >= 0 && std::stoi(tokens[0]) <= 199)
        {
          tokens.clear();
          std::getline(dataFile,linestr);
          istringstream iss(linestr);
          copy(istream_iterator<string>(iss), istream_iterator<string>(), back_inserter(tokens));

          if(std::stoi(tokens[0]) >= 0 && std::stoi(tokens[0]) <= 199)
          {
            cout << "    NUCLEON SCATTERING! " << tokens[0] << endl;
            npart = npart + 2;
          }
          else
          {
            continue;
          }
        }
      }

      if (!dataFile) break;

    }//close while over dataFile

    float mean_npart = (float) npart/numevent;
    cout << "<Npart> = " << npart << " / " << numevent << " =  " << mean_npart << endl;
 
}//close void()
