#include "functions.h"
#include <vector>
#include <fstream>
#include <iostream>
#include <map>
#include <set>

using namespace std;

/////////////////////////////////////////////////////////////////////////////////
/////This code has been adopted and slightly modified from the study of Sung, J. et al. Global metabolic interaction network of 
/////the human gut microbiota for context-specific community-scale analysis. Nat. Commun. 8, 15393  
/////doi:10.1038/ncomms15393 (2017).
////////////////////////////////////
/////STEP1
/////////////////////////////////////////////////////////////////////////////////
void met_microbe_coupling(std::string datafolder, std::string outputfolder)
{
	/////////////////////////////////////////////////
	// READING FILES AND INITIALIZING VARIABLES
	////////////////////////////////////////////////

	ifstream inputF;
	vector<int> microbeIDs; //list of microbe IDs
	vector<double> relab;	//relative abundance of microbes in healthy sample
	vector<int> smallmets;	//list of small molecules to be used in this analysis
	vector<int> macros;		//list of macromolecules to be used in this analysis

	set<int> smallmetset;
	set<int> macroset;

	map<int, int> microbeIDtoInd; //key: ID of microbe, value: corresponding index in "microbeIDs"
	map<int, int> metsIDtoInd; //key: ID of small molecule, value: corresponding index in "smallmets"
	map<int, int> macroIDtoInd; //key: ID of macromolecule, value: corresponding index in "macros"

	vector<vector<int> > producermap; //key: small molecules index in "smallmets", value: direct producers index in "microbeIDs"
	vector<vector<int> > consumermap; //key: small molecules index in "smallmets", value: direct consumers index in "microbeIDs"
	map<int, set<int> > macrodegraders; //key: macromolecules index in "macros", value: degraders index in "microbeIDs"
	vector<vector<int> > metmacromap; //key: small molecules index in "smallmets", value: corresponding macromolecules index in "macros"


	//read relative abundance of microbes
	int N = 0;
	int currID, diffexp;
	double abH, abT2D;
	inputF.open(datafolder+"/microbe_info.txt");
	if(!inputF.is_open())
	{
		cout<<"NOT OPENED"<<endl;
	}

	inputF >> N;
	for (int i=0; i < N; ++i)
	{
		inputF >> currID >> abH;
		microbeIDs.push_back(currID);
		relab.push_back(abH);
		microbeIDtoInd[currID] = i;
	}

	inputF.close();
	
	//read the list of small molecules
	int currMetID;
	
	inputF.open(datafolder+"/small_molecule_info.txt");
	if(!inputF.is_open())
	{
		cout<<"NOT OPENED"<<endl;
	}
	inputF >> N;
	for (int i=0; i < N; ++i)
	{
		inputF >> currMetID;
		smallmets.push_back(currMetID);
		metsIDtoInd[currMetID] = i;
	}
	inputF.close();

	smallmetset.insert(smallmets.begin(), smallmets.end());
	producermap.resize(N, vector<int>());
	consumermap.resize(N, vector<int>());
	metmacromap.resize(N, vector<int>());

	//read the list of macromolecules
	int macroID;

	inputF.open(datafolder+"/macromolecule_info.txt");
	if(!inputF.is_open())
	{
		cout<<"NOT OPENED"<<endl;
	}
	inputF >> N;
	for (int i=0; i < N; ++i)
	{
		inputF >> macroID;
		macros.push_back(macroID);
		macroIDtoInd[macroID] = i;
	}
	macroset.insert(macros.begin(), macros.end());
	inputF.close();

	//read network information
	int ID1, ID2, currtype;

	inputF.open(datafolder+"/network_info.txt");
	if(!inputF.is_open())
	{
		cout<<"NOT OPENED"<<endl;
	}
	inputF >> N;

	for(int i=0; i < N; ++i)
	{
		inputF >> ID1 >> ID2 >> currtype;
		
		if (currtype == 2 || currtype == 4)
		{//ID1: metabolite, ID2: microbe, Import case
			auto metIt = metsIDtoInd.find(ID1);
			auto microbeIt = microbeIDtoInd.find(ID2);
			if(metIt != metsIDtoInd.end() && microbeIt != microbeIDtoInd.end())
				consumermap[metIt->second].push_back(microbeIt->second);
		}
		if (currtype == 3 || currtype == 4)
		{//ID1: metabolite, ID2: microbe, Export case
			auto metIt = metsIDtoInd.find(ID1);
			auto microbeIt = microbeIDtoInd.find(ID2);
			if(metIt != metsIDtoInd.end() && microbeIt != microbeIDtoInd.end())
				producermap[metIt->second].push_back(microbeIt->second);
		}
		if (currtype == 5)
		{//ID1: macromolecule, ID2: species, degrading
			auto macroIt = macroIDtoInd.find(ID1);
			auto microbeIt = microbeIDtoInd.find(ID2);
			if(macroIt != macroIDtoInd.end() && microbeIt != microbeIDtoInd.end())
				macrodegraders[macroIt->second].insert(microbeIt->second);
			
		}
		if (currtype == 7)
		{//ID1: macromolecule, ID2: small molecules, macro-met pairing
			auto macroIt = macroIDtoInd.find(ID1);
			auto metIt = metsIDtoInd.find(ID2);

			if(macroIt != macroIDtoInd.end() &&	metIt != metsIDtoInd.end() )
				metmacromap[metIt->second].push_back(macroIt->second);
		}		
	}

	inputF.close();

	////////////////////////////////////////////////////////////////////////
	// CALCULATING PRODUCER AND CONSUMER TERMS FOR METABOLITE-MICROBE PAIRS
	////////////////////////////////////////////////////////////////////////
	
	double thre1 = 3e-2; //Threshold for microbial population. thre1=1e-3=0.001
						//Equivalent to theta1 in Supplementary Information
	int currB; //1 if a given metabolite comes from degradation of a macromolecule by a given microbe
			   //0 otherwise
	int currmet, currmicrobe;

	//calculating total population of degraders of each macromolecule
	vector<double> degraderpop(macros.size(), 0.0); 
			//key: macromolecules index in "macros", value: total population of degraders
	
	for(auto it = macrodegraders.begin(); it != macrodegraders.end(); ++it)
	{
		int currmacro = it->first;
		set<int> currdegraders = it->second;
		double totalpop = 0;

		for(auto it2 = currdegraders.begin(); it2 != currdegraders.end(); ++it2)
			totalpop += relab[*it2];
		
		degraderpop[currmacro] = totalpop;
	}

	//For each microbe, allocate the list of macromolecules (m's) such that
	//the microbe degrades m and degraderpop[m] >= thre1 
	map<int, vector<int> > degraderToMacro;
	for(auto it = macrodegraders.begin(); it != macrodegraders.end(); ++it)
	{
		int currmacro = it->first;
		set<int> currdegraders = it->second;

		if (degraderpop[currmacro] < thre1)
			continue;
		for(auto it2 = currdegraders.begin(); it2 != currdegraders.end(); ++it2)
			degraderToMacro[*it2].push_back(currmacro);
	}
	
	//open output file

	ofstream met_microbe_F(outputfolder+"/met_microbe_coupling.txt");
	met_microbe_F.precision(16);

	//calculate producer and consumer terms
	for(int i=0; i < smallmets.size(); ++i)
	{
		vector<double> currProdTerm(microbeIDs.size(), 0.0);
		vector<double> currConTerm(microbeIDs.size(), 0.0);

		vector<int> currProducers = producermap[i];
		vector<int> currConsumers = consumermap[i];
		set<int> currMacros;

		double totalconspop = 0.0;
		double totalprodpop = 0.0;

		int currNumMacro = 0; // # of macromolecules whose total population of degraders is >= thre1

		currmet = smallmets[i];
		currB = 0;
		
		for(int j=0; j < currConsumers.size(); ++j)
			totalconspop += relab[currConsumers[j]];
		for(int j=0; j < currProducers.size(); ++j)
			totalprodpop += relab[currProducers[j]];
		for(auto it2 = metmacromap[i].begin(); it2 != metmacromap[i].end(); ++it2)
		{
			if (degraderpop[*it2] >= thre1)
			{
				currMacros.insert(*it2);
				++currNumMacro;
			}
		}

		//calculate consumer term
		for(int j=0; j < currConsumers.size(); ++j)
		{
			int currind = currConsumers[j];

			if (totalconspop < thre1)
				currConTerm[currind] = 0.0;
			else
				currConTerm[currind] = relab[currind]/totalconspop;
			
		}
		//calculate producer term
		if (currNumMacro == 0)
		{//non macromolecule case
			currB = 0;
			
			for(int j=0; j < currProducers.size(); ++j)
			{
				int currind = currProducers[j];

				if (totalprodpop < thre1)
					currProdTerm[currind] = 0.0;
				else
					currProdTerm[currind] = relab[currind]/totalprodpop;
			}
		}
		else
		{//macromolecule case
			currB = 1;
			set<int> degradersK;
			
			for(auto it2 = currMacros.begin(); it2 != currMacros.end(); ++it2)
				degradersK.insert(macrodegraders[*it2].begin(),
										macrodegraders[*it2].end());


			for(auto it2 = degradersK.begin(); it2 != degradersK.end(); ++it2)
			{		
				int currind = *it2;
				vector<int> macroOfCurrDegrader = degraderToMacro[currind];
				vector<int> macrosCurrSp; 
							//macromolecules (m's) such that total population of degraders of m is >= thre1
							//,j is one of the degraders, and i can come out from m.
				for(int k=0; k < macroOfCurrDegrader.size(); ++k)
				{
					if(currMacros.find(macroOfCurrDegrader[k]) != currMacros.end())
					{
						macrosCurrSp.push_back(macroOfCurrDegrader[k]);
					}
				}
				
				if (relab[currind] >= thre1)
				{
					currProdTerm[currind] = static_cast<double>(macrosCurrSp.size())/currNumMacro;

				}
				else
				{
					double prodterm_ij = 0;
					for(int k=0; k < macrosCurrSp.size(); ++k)
					{
						prodterm_ij += relab[currind]/degraderpop[macrosCurrSp[k]];
					}
					currProdTerm[currind] =  prodterm_ij / currNumMacro;
					
				}
				
			}
		}

		//write current information into file
		//Format: B, Producer Terms, Consumer Terms
		met_microbe_F<<currB<<'\t';
		for(int j=0; j < currProdTerm.size(); ++j)
			met_microbe_F<<currProdTerm[j]<<'\t';
		for(int j=0; j < currConTerm.size(); ++j)
			met_microbe_F<<currConTerm[j]<<'\t';
		met_microbe_F<<endl;
	}

	met_microbe_F.close();

}
