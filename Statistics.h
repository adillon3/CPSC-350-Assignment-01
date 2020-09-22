/*******************************
 * Andrew Dillon
 * 2382400
 * CPSC 350
 * Assignment 01
 *******************************/

#ifndef STATISTICS_H_
#define STATISTICS_H_

#include <iostream>
using namespace std;
#include <string>
#include <fstream>
#include<cmath>
#include <iomanip>


class Statistics
{
public:
	//cosntructor and destructor
	Statistics();
	~Statistics();

	void ProcessFile(string fileName, ofstream& oFile);
	int GetNumLines(ifstream& inFile);
	int SumNucleotideCounts();
	float Mean(int sum, int numLines);
	float Variance(ifstream& inFile, int sum, float mean);
	float StandardDeviation(float variance) throw(exception);
	void CountNucleotides(ifstream& inFile);
	void PrintProbabilityOfNucleotide(ofstream& oFile, int sum);
	string ToUpperString(string myString);
	void PrintProbabilityOfNucleotideBigram(ifstream& inFile, ofstream& oFile, int sum);
	void InitilizePrivateVariables();
	void Append1000Strings(ofstream& oFile, float standardDev, float mean);
	bool CheckLineForInvalidNucleotides(string myString);
	
private:
	int aCount;
	int cCount;
	int gCount;
	int tCount;

	float aProb;
	float cProb;
	float gProb;
	float tProb;
};






#endif /* STATISTICS_H_ */
