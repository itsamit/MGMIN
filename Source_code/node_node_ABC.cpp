#include "functions.h"
#include "matrix.h"
#include <vector>
#include <fstream>
#include <iostream>
#include <map>
#include <set>

using namespace std;

///////////////////////////////////////////////////////////////////////////
//////This code has been adopted and slightly modified from the study of Sung, J. et al. Global metabolic interaction network of 
/////the human gut microbiota for context-specific community-scale analysis. Nat. Commun. 8, 15393  
////doi:10.1038/ncomms15393 (2017)
/////Step2
///////////////////////////////////////////////////////////////////////////
void node_node_ABC(std::string datafolder, std::string outputfolder)
{
	/////////////////////////////////////////////////
	// READING FILES AND INITIALIZING VARIABLES
	////////////////////////////////////////////////

	ifstream inputF;
	vector<int> microbeIDs; //list of microbe IDs
	vector<int> smallmets;	//list of small molecules to be used in this analysis
	vector<int> groupIDs;	//list of microbial group IDs
	vector<set<int> > groupInfo; //list of microbes in each group
	vector<double> relab;	//relative abundance of microbes in healthy sample
	vector<double> grouppop;//relative abundance of microbial groups in healthy sample

	set<int> smallmetset;
	set<int> macroset;

	map<int, int> microbeIDtoInd; //key: ID of microbe, value: corresponding index in "microbeIDs"
	map<int, int> metsIDtoInd; //key: ID of small molecule, value: corresponding index in "smallmets"
	map<int, int> macroIDtoInd; //key: ID of macromolecule, value: corresponding index in "macros"

	vector<vector<int> > consumingmets; //key: small molecules index in "smallmets", value: direct consumers index in "microbeIDs"
	
	vector<int> metBs; // 1: macromolecule case, 0: non macromolecule case
	vector<vector<double> > producerTerms; //row: metabolite, col: microbe
	vector<vector<double> > consumerTerms; //row: metabolite, col: microbe
	

	//read the list of microbes
	int N = 0;
	int currID, diffexp;
	double abH, abT2D;
	inputF.open(datafolder+"/microbe_info.txt");
	if(!inputF.is_open())
	{
		cout<<"MICROBE INFO NOT OPENED"<<endl;
	}
	inputF >> N;
	for (int i=0; i < N; ++i)
	{
		inputF >> currID >> abH;
		microbeIDs.push_back(currID);
		relab.push_back(abH);
		microbeIDtoInd[currID] = i;
	}

	consumingmets.resize(N, vector<int>());

	inputF.close();

	//read group data
	inputF.open(datafolder+"/group_info.txt");
	if(!inputF.is_open())
	{
		cout<<"GROUP INFO NOT OPENED"<<endl;
	}
	inputF >> N;
	groupIDs.resize(N, 0);
	grouppop.resize(N, 0.0);
	groupInfo.resize(N, set<int>());
	for(int i=0; i < N; ++i)
	{
		int groupsize, currmic;
		double currpop = 0.0;

		inputF>>groupIDs[i]>>groupsize;
		
		for(int j=0; j < groupsize; ++j)
		{
			inputF >> currmic;
			groupInfo[i].insert(microbeIDtoInd[currmic]);
			currpop += relab[microbeIDtoInd[currmic]];
		}
		grouppop[i] = currpop;
	}

	inputF.close();

	//read the list of small molecules
	int currMetID;
	
	inputF.open(datafolder+"/small_molecule_info.txt");
	if(!inputF.is_open())
	{
		cout<<"SMALL MOLECULE FILE NOT OPENED"<<endl;
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
	metBs.resize(N, 0);
	producerTerms.resize(N, vector<double>(microbeIDs.size(), 0));
	consumerTerms.resize(N, vector<double>(microbeIDs.size(), 0));

	//read network information
	int ID1, ID2, currtype;

	inputF.open(datafolder+"/network_info.txt");
	if(!inputF.is_open())
	{
		cout<<"NETWORK FILE NOT OPENED"<<endl;
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
				consumingmets[microbeIt->second].push_back(metIt->second);
		}				
	}
	inputF.close();

	

	//read met_microbe_coupling data
	inputF.open(outputfolder+"/met_microbe_coupling.txt");
	if(!inputF.is_open())
	{
		cout<<"MET MIC COUPLING NOT OPENED"<<endl;
	}
	for(int i=0; i < smallmets.size(); ++i)
	{
		inputF >> metBs[i];

		for(int j=0; j < microbeIDs.size(); ++j)
			inputF >> producerTerms[i][j];
		
		for(int j=0; j < microbeIDs.size(); ++j)
			inputF >> consumerTerms[i][j];
	}

	inputF.close();

	////////////////////////////////////////////////////////////////////////
	// CALCULATING A, B, and C where W = A*alpha + B*beta + C
	////////////////////////////////////////////////////////////////////////
	
	ofstream node_node_file(outputfolder+"/node_node_ABC.txt");
	node_node_file.precision(12);
	matrix2D micMicA(microbeIDs.size(), microbeIDs.size()),
				micMicB(microbeIDs.size(), microbeIDs.size()),
				micMicC(microbeIDs.size(), microbeIDs.size());

	//microbe->microbe
	for(int i=0; i < microbeIDs.size(); ++i)
	{
		for(int j=0; j < microbeIDs.size(); ++j)
		{//calculating W_ij
			micMicA(i, j) = 0;
			micMicB(i, j) = 0;
			micMicC(i, j) = 0;

			if(i==j)
				continue;

			for(auto itk = consumingmets[j].begin(); itk != consumingmets[j].end(); ++itk)
			{
				if (metBs[*itk] == 0)
					micMicA(i, j) += producerTerms[*itk][i];
				else
					micMicB(i, j) += producerTerms[*itk][i];
				micMicC(i, j) -= consumerTerms[*itk][i];
			}
		}
	}

	//microbe->microbe or group
	for(int i=0; i < microbeIDs.size(); ++i)
	{
		//write microbe->microbe
		for(int j=0; j < microbeIDs.size(); ++j)
		{
			node_node_file << micMicA(i, j) << '\t'
							<< micMicB(i, j) << '\t'
							<< micMicC(i, j) << '\t';
		}

		for(int j=0; j < groupIDs.size(); ++j)
		{//W_iG
			double currA, currB, currC;
			currA = 0.0;
			currB = 0.0;
			currC = 0.0;

			if (grouppop[j] > 0)
			{
				for(auto itk = groupInfo[j].begin(); itk != groupInfo[j].end(); ++itk)
				{
					currA += relab[*itk] / grouppop[j] * micMicA(i, *itk);
					currB += relab[*itk] / grouppop[j] * micMicB(i, *itk);
					currC += relab[*itk] / grouppop[j] * micMicC(i, *itk);
				}
			}
			else
			{
				for(auto itk = groupInfo[j].begin(); itk != groupInfo[j].end(); ++itk)
				{
					currA += micMicA(i, *itk);
					currB += micMicB(i, *itk);
					currC += micMicC(i, *itk);
				}
				currA /= groupInfo[j].size();
				currB /= groupInfo[j].size();
				currC /= groupInfo[j].size();
			}

			node_node_file << currA << '\t'
							<< currB << '\t'
							<< currC << '\t';
		}

		node_node_file<<endl;
	}

	//write group -> microbe and group -> group
	for(int i=0; i < groupIDs.size(); ++i)
	{//W_Gi
		//write group -> microbe
		for(int j = 0; j < microbeIDs.size(); ++j)
		{
			double currA, currB, currC;
			currA = 0.0;
			currB = 0.0;
			currC = 0.0;

			for(auto itk = groupInfo[i].begin(); itk != groupInfo[i].end(); ++itk)
			{
				currA += micMicA(*itk, j);
				currB += micMicB(*itk, j);
				currC += micMicC(*itk, j);
			}

			node_node_file << currA << '\t'
							<< currB << '\t'
							<< currC << '\t';
		}	


		//write group -> group
		for(int j=0; j < groupIDs.size(); ++j)
		{//W_Gr
			double currA, currB, currC;
			vector<int> currTargetGroup;
			double currTargetPop = 0.0;
			
			currA = 0.0;
			currB = 0.0;
			currC = 0.0;

			if(i == j)
			{
				node_node_file << currA << '\t'
							<< currB << '\t'
							<< currC << '\t';
				continue;
			}
			
			for(auto itk = groupInfo[j].begin(); itk != groupInfo[j].end(); ++itk)
			{
				if(groupInfo[i].find(*itk) == groupInfo[i].end())
					currTargetGroup.push_back(*itk);
			}
			
			for(auto itk = currTargetGroup.begin(); itk != currTargetGroup.end(); ++itk)
			{
				currTargetPop += relab[*itk];
			}


			if (currTargetPop > 0)
			{
				for(auto itk1 = groupInfo[i].begin(); itk1 != groupInfo[i].end(); ++itk1)
				{
					for(auto itk2 = currTargetGroup.begin(); itk2 != currTargetGroup.end(); ++itk2)
					{
						currA += relab[*itk2] / currTargetPop * micMicA(*itk1, *itk2);
						currB += relab[*itk2] / currTargetPop * micMicB(*itk1, *itk2);
						currC += relab[*itk2] / currTargetPop * micMicC(*itk1, *itk2);
					}
				}
			}
			else
			{
				for(auto itk1 = groupInfo[i].begin(); itk1 != groupInfo[i].end(); ++itk1)
				{
					for(auto itk2 = currTargetGroup.begin(); itk2 != currTargetGroup.end(); ++itk2)
					{
						currA += micMicA(*itk1, *itk2);
						currB += micMicB(*itk1, *itk2);
						currC += micMicC(*itk1, *itk2);
					}
				}
				currA /= currTargetGroup.size();
				currB /= currTargetGroup.size();
				currC /= currTargetGroup.size();
			}

			node_node_file << currA << '\t'
							<< currB << '\t'
							<< currC << '\t';
		}
		
		node_node_file<<endl;
	}

	node_node_file.close();

	return;
}
