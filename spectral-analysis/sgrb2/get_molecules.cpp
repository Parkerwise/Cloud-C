/////////////////////////////////////////////////////////////////////////////
// Name:        get_molecules.cpp
// Purpose:		output python list with molecules present within sgrB2
// Author:      Parker Wise
// Created:     2024-06-25
/////////////////////////////////////////////////////////////////////////////
#include <iostream>
#include <algorithm>
#include <ostream>
#include <vector>
#include <fstream>
#include <string>

double GetFrequencies(std::string line){ 
	std::vector<char> bytes(line.begin(), line.end());
	bytes.push_back('\0');
	std::string FrequencyString(bytes.begin()+60,bytes.begin()+70);
	double Frequency = std::stod(FrequencyString,0);
	return Frequency;
}
std::string GetTransitions(std::string line){
	std::vector<char> bytes(line.begin(), line.end());
	bytes.push_back('\0');
	std::string Transition(bytes.begin(),bytes.begin()+21);
	std::string Molecule = Transition.substr(0, Transition.find(",", 0));
	Molecule.erase(std::remove_if(Molecule.begin(), Molecule.end(), isspace), Molecule.end());
	return Molecule;
}
int main(){
std::string DataFile;
std::cout<<"Name of Input File \n";
std::cin>> DataFile;
float StartFreq;
std::cout<<"Starting Freq (MHz)?: \n";
std::cin >> StartFreq;
float EndFreq;
std::cout<<"Ending Freq (MHz)?: \n";
std::cin >> EndFreq;
std::ifstream MoleculeList(DataFile);
std::string line;
std::vector<double> AllFrequencies;
std::vector<std::string> AllMolecules;
if (!MoleculeList.is_open()) { 
	std::cerr << "Error opening the file!" << std::endl; 
	return 1;
	}
while (getline(MoleculeList,line)){
	float IndidualFrequency=GetFrequencies(line);
	if (IndidualFrequency > StartFreq && IndidualFrequency < EndFreq){
	std::string IndividualMolecule=GetTransitions(line);
	AllMolecules.push_back(IndividualMolecule);
	}
sort(AllMolecules.begin(),AllMolecules.end());
AllMolecules.erase(std::unique(AllMolecules.begin(),AllMolecules.end()),AllMolecules.end());
}
MoleculeList.close();
std::string SelectedMolecules = "[";
std::string FinalMolecule=AllMolecules.back();
for (std::string i : AllMolecules){
	SelectedMolecules.append(i);
	if (i != FinalMolecule){
	SelectedMolecules.append(",");
	}
}
SelectedMolecules.append("]");
std::cout<<SelectedMolecules;
return 0;
}

