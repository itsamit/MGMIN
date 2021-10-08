#pragma once

#include <string>

//Functions for each steps

//step1
void met_microbe_coupling(std::string datafolder, std::string outputfolder);

//step2
void node_node_ABC(std::string datafolder, std::string outputfolder);

//step3
void printWs(std::string datafolder, std::string outputfolder, double alpha, double beta, double thre);
