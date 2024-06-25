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
	//Values are specified by bytes in dat files
	//so our string must be converted to bytes
	std::vector<char> bytes(line.begin(), line.end());
	bytes.push_back('\0');
	//Frequencies are from byte 61-70
	std::string FrequencyString(bytes.begin()+60,bytes.begin()+70);
	double Frequency = std::stod(FrequencyString,0);
	return Frequency;
}
//The same process is repeated to get transitions (bytes 1-21)
std::string GetTransitions(std::string line){
	std::vector<char> bytes(line.begin(), line.end());
	bytes.push_back('\0');
	std::string Transition(bytes.begin(),bytes.begin()+21);
	//removes vibrational states since they are specified after a comma
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
while (getline(MoleculeList,line)){ //populates vector with molecules in freq range
	float IndidualFrequency=GetFrequencies(line);
	if (IndidualFrequency > StartFreq && IndidualFrequency < EndFreq){
	std::string IndividualMolecule=GetTransitions(line);
	AllMolecules.push_back(IndividualMolecule);
	}
//Sort and erase remove duplicates (different transitions, same molecule)
sort(AllMolecules.begin(),AllMolecules.end());
AllMolecules.erase(std::unique(AllMolecules.begin(),AllMolecules.end()),AllMolecules.end());
}
MoleculeList.close();
//outputs string formatted as a python list for use in xclass
std::string SelectedMolecules = "[";
std::string FinalMolecule=AllMolecules.back();
for (std::string i : AllMolecules){
	SelectedMolecules.append(i);
	if (i != FinalMolecule){ //avoids comma after final entry
	SelectedMolecules.append(",");
	}
}
SelectedMolecules.append("]");
std::cout<<SelectedMolecules;
return 0;
}

