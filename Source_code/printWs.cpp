#include "functions.h"
#include "matrix.h"
#include <vector>
#include <fstream>
#include <iostream>
#include <map>
#include <cmath>
#include <set>

using namespace std;

///////////////////////////////////////////////////////////////////////////
/////This code has been adopted and slightly modified from the study of Sung, J. et al. Global metabolic interaction network of 
/////the human gut microbiota for context-specific community-scale analysis. Nat. Commun. 8, 15393  
/////doi:10.1038/ncomms15393 (2017).
//////////////////
/////Step3
///////////////////////////////////////////////////////////////////////////
void printWs(std::string datafolder, std::string outputfolder, double alpha, double beta, double thre)
{
	/////////////////////////////////////////////////
	// READING FILES AND INITIALIZING VARIABLES
	////////////////////////////////////////////////

	ifstream inputF;
	vector<int> nodeIDs;	//list of microbe and microbial group IDs (here called nodes)
	vector<int> diffExpNodes;//list of indices of differentially expressed nodes in "nodeIDs" 
	vector<int> diffExp;	//differentially expressed? 1: in Healthy, -1: in T2D, 0: nowhere

	int numMic = 0;

	//read the list of microbes
	int N = 0;
	int currID, diffexp;
	double abH, abT2D;
	
	inputF.open(datafolder+"/microbe_info.txt");
	if(!inputF.is_open())
	{
		cout<<"NOT OPENED"<<endl;
	}
	inputF >> N;
	numMic = N;
	
	for (int i=0; i < N; ++i)
	{
		inputF >> currID >> abH;
		nodeIDs.push_back(currID);
	}
	
	inputF.close();
	
	//read group data
	inputF.open(datafolder+"/group_info.txt");
	if(!inputF.is_open())
	{
		cout<<"NOT OPENED"<<endl;
	}
	inputF >> N;
	
	for(int i=0; i < N; ++i)
	{
		int groupID,groupsize, currmic;
		double currpop = 0.0;

		inputF>>groupID>>groupsize;
		
		for(int j=0; j < groupsize; ++j)
			inputF >> currmic;
		nodeIDs.push_back(groupID);
	}

	inputF.close();

	//read A, B, C
	int numNodes = nodeIDs.size();
	matrix2D nodeA(numNodes, numNodes), nodeB(numNodes, numNodes), 
			nodeC(numNodes, numNodes); 

	inputF.open(outputfolder+"/node_node_ABC.txt");
	if(!inputF.is_open())
	{
		cout<<"NOT OPENED"<<endl;
	}
	
	for(int i=0; i < numNodes; ++i)
	{
		for(int j=0; j < numNodes; ++j)
		{
			inputF >> nodeA(i, j) >> nodeB(i, j) >> nodeC(i, j);
		}
	}
	inputF.close();

	////////////////////////////////////////////////////////////////////////
	// CALCULATING AND PRINTING W's IN MATRIX FORM FOR EVERY NODE PAIR
	//     NOTE: W = A*alpha + B*beta + C
	////////////////////////////////////////////////////////////////////////
	ofstream Wfile;
	double currW_ij;
		
	Wfile.open(outputfolder+"/Wmatrix_a1.0_b1.0_thre1e-4.txt");
	
	for(int i=0; i < nodeIDs.size(); ++i)
	{
		for(int j=0; j < nodeIDs.size(); ++j)
		{	

			currW_ij = nodeA(i, j) * alpha + 
						nodeB(i, j) * beta +
						nodeC(i, j);
			if (fabs(currW_ij) < thre)
				currW_ij = 0.0;

			Wfile << currW_ij <<'\t';
		}
		Wfile<<endl;
	}
	Wfile.close();
}
