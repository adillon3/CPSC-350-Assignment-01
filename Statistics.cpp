/*******************************
 * Andrew Dillon
 * 2382400
 * CPSC 350
 * Assignment 01
 *******************************/

#include "Statistics.h"


Statistics::Statistics()
{
	InitilizePrivateVariables();
}

Statistics::~Statistics()
{

}

void Statistics::InitilizePrivateVariables()
{
	aCount = 0;
	cCount = 0;
	gCount = 0;
	tCount = 0;

	aProb = 0.0;
	cProb = 0.0;
	gProb = 0.0;
	tProb = 0.0;
}

void Statistics::ProcessFile(string fileName, ofstream& oFile)
{
	ifstream inFile;
	inFile.open(fileName);

	int numLines = GetNumLines(inFile);

	InitilizePrivateVariables();

	if(numLines <= 0)
	{
		oFile << "No DNA strings found in the file" << endl;

		return;
	}

	//DO THE OTHER FUNCTIONS
	//resetting to front of file
	inFile.seekg (0, inFile.beg);
	CountNucleotides(inFile);

	int sum 	 	= SumNucleotideCounts();
	float mean 	 	= Mean(sum, numLines);

	try
	{
		//resetting to front of file
		inFile.seekg (0, inFile.beg);
		float varaince 	= Variance(inFile, sum, mean);

		//resetting to front of file
		inFile.seekg (0, inFile.beg);
		float standardDev	= StandardDeviation(varaince);


		oFile << fixed << showpoint << setprecision(2);

		//This data to file
		oFile << "Basic Statistics\n";
		oFile << " - Sum               : " << sum << endl;
		oFile << " - Mean              : " << mean << endl;
		oFile << " - Variance          : " << varaince << endl;
		oFile << " - Standard Deviation: " << standardDev << endl << endl;

		//resetting output to standard
		oFile.unsetf(ios::fixed);
		oFile << noshowpoint << setprecision(6);

		PrintProbabilityOfNucleotide(oFile, sum);
		PrintProbabilityOfNucleotideBigram(inFile, oFile, sum);
		Append1000Strings(oFile, standardDev, mean);

		oFile << endl << endl;

		inFile.close();
	}
	catch(const char* msg)
	{
		return;
	}
}

int Statistics::GetNumLines(ifstream& inFile)
{
	int numLines = 0;

	string temp;
	while(inFile && !inFile.eof())
	{
		getline(inFile, temp);

		if(CheckLineForInvalidNucleotides(temp))
		{
			++numLines;
		}
	}

	return numLines;
}

int Statistics::SumNucleotideCounts()
{
	return aCount + cCount + gCount + tCount;
}

float Statistics::Mean(int sum, int numLines)
{
	if(numLines == 0)
	{
		cout << "Mean could not be found." << endl;
	}

	return (float)sum / numLines;
}

float Statistics::Variance(ifstream& inFile, int sum, float mean)
{
	// sum cannot be 1 because to find the variance we must divide by the sum - 1
	//	if sum is 1, this would result in a divid by 0 error
	if(sum == 1)
	{
		throw "Variance could not be found.\n";
	}


	//The summation of each item in the data set minus the mean squared
	//The sum of (x - mean)^2
	float summation = 0;
	string stringTemp;
	int intTemp;

	while(inFile && !inFile.eof())
	{
		getline(inFile, stringTemp);

		if(CheckLineForInvalidNucleotides(stringTemp))
		{
			intTemp = stringTemp.length() - mean;
			intTemp = intTemp * intTemp;

			summation += intTemp;
		}

	}

	return summation / (sum - 1);
}

float Statistics::StandardDeviation(float variance) throw(exception)
{
	if(variance < 0)
	{
		throw "Standard Deviation could not be found.\n";
	}

	return sqrt(variance);
}

void Statistics::CountNucleotides(ifstream& inFile)
{
	string temp;

	//count number of each Nucleotide
	while(inFile && !inFile.eof())
	{
		getline(inFile, temp);

		if(CheckLineForInvalidNucleotides(temp))
		{
			int length = temp.length();

			for(int i = 0; i < length; ++i)
			{
				if(toupper(temp[i]) == 'A')
				{
					++aCount;
				}
				if(toupper(temp[i]) == 'C')
				{
					++cCount;
				}
				if(toupper(temp[i]) == 'G')
				{
					++gCount;
				}
				if(toupper(temp[i]) == 'T')
				{
					++tCount;
				}
			}//end for(int i = 0; i < length; ++i)
		}//end if(CheckLineForInvalidNucleotides(temp))
	}//end while(inFile && !inFile.eof())
}//end function

void Statistics::PrintProbabilityOfNucleotide(ofstream& oFile, int sum)
{
	oFile << "Nucleotide Probabilities" << endl;

	//formatting output
	oFile << fixed << showpoint << setprecision(2);

	aProb = (float)aCount / sum * 100;
	cProb = (float)cCount / sum * 100;
	gProb = (float)gCount / sum * 100;
	tProb = (float)tCount / sum * 100;


	if(sum != 0)
	{
		oFile << " - adenine : " << setw(6) << aProb << "%\n";
		oFile << " - cytosine: " << setw(6) << cProb << "%\n";
		oFile << " - guanine : " << setw(6) << gProb << "%\n";
		oFile << " - thymine : " << setw(6) << tProb << "%\n";
	}
	else
	{
		oFile << " - adenine : " << "Could not calculate" << endl;
		oFile << " - cytosine: " << "Could not calculate" << endl;
		oFile << " - guanine : " << "Could not calculate" << endl;
		oFile << " - thymine : " << "Could not calculate" << endl;
	}

	//resetting output to standard
	oFile.unsetf(ios::fixed);
	oFile << noshowpoint << setprecision(6);

	oFile << endl;
}

string Statistics::ToUpperString(string myString)
{
	int length = myString.length();


	for(int i = 0; i < length; ++i)
	{
		if (islower(myString[i]))
		{
			myString[i] = myString[i] - 32;
		}
	}

	return myString;
}






void Statistics::PrintProbabilityOfNucleotideBigram(ifstream& inFile, ofstream& oFile, int sum)
{
	if(sum == 0)
	{
		oFile << "Nucleotide Bigrams could not be calculated\n";

		return;
	}

	int aaCount = 0;
	int acCount = 0;
	int agCount = 0;
	int atCount = 0;
	int caCount = 0;
	int ccCount = 0;
	int cgCount = 0;
	int ctCount = 0;
	int gaCount = 0;
	int gcCount = 0;
	int ggCount = 0;
	int gtCount = 0;
	int taCount = 0;
	int tcCount = 0;
	int tgCount = 0;
	int ttCount = 0;

	string temp;

	while(inFile && !inFile.eof())
	{
		getline(inFile, temp);
		if(CheckLineForInvalidNucleotides(temp))
		{
			temp = ToUpperString(temp);

			int length = temp.length();
			int j = 1;
			//counting each type of bigram
			for(int i = 0; j < length; ++i && ++j)
			{
				if(temp[i] == 'A' && temp[j] == 'A')
				{
					++aaCount;
				}
				else if(temp[i] == 'A' && temp[j] == 'C')
				{
					++acCount;
				}
				else if(temp[i] == 'A' && temp[j] == 'G')
				{
					++agCount;
				}
				else if(temp[i] == 'A' && temp[j] == 'T')
				{
					++atCount;
				}
				else if(temp[i] == 'C' && temp[j] == 'A')
				{
					++caCount;
				}
				else if(temp[i] == 'C' && temp[j] == 'C')
				{
					++ccCount;
				}
				else if(temp[i] == 'C' && temp[j] == 'G')
				{
					++cgCount;
				}
				else if(temp[i] == 'C' && temp[j] == 'T')
				{
					++ctCount;
				}
				else if(temp[i] == 'G' && temp[j] == 'A')
				{
					++gaCount;
				}
				else if(temp[i] == 'G' && temp[j] == 'C')
				{
					++gcCount;
				}
				else if(temp[i] == 'G' && temp[j] == 'G')
				{
					++ggCount;
				}
				else if(temp[i] == 'G' && temp[j] == 'T')
				{
					++gtCount;
				}
				else if(temp[i] == 'T' && temp[j] == 'A')
				{
					++taCount;
				}
				else if(temp[i] == 'T' && temp[j] == 'C')
				{
					++tcCount;
				}
				else if(temp[i] == 'T' && temp[j] == 'G')
				{
					++tgCount;
				}
				else if(temp[i] == 'T' && temp[j] == 'T')
				{
					++ttCount;
				}
			}//end if(CheckLineForInvalidNucleotides(temp))
		}//end (CheckLineForInvalidNucleotides(temp))
	}//end while(inFile && !inFile.eof())


	//sum up all the bigram counts
	int bigramSum = aaCount + acCount + agCount + atCount +
					caCount + ccCount + cgCount + ctCount +
					gaCount + gcCount + ggCount + gtCount +
					taCount + tcCount + tgCount + ttCount;

	//formatting output
	oFile << fixed << showpoint << setprecision(2);

	oFile << "Nucleotide Bigram Probabilities\n";
	oFile << " - adenine-adenine  : ";
	oFile << setw(6) << (float)aaCount / bigramSum * 100 << "%\n";
	oFile << " - adenine-cytosine : ";
	oFile << setw(6) << (float)acCount / bigramSum * 100 << "%\n";
	oFile << " - adenine-guanine  : ";
	oFile << setw(6) << (float)agCount / bigramSum * 100 << "%\n";
	oFile << " - adenine-thymine  : ";
	oFile << setw(6) << (float)atCount / bigramSum * 100 << "%\n";

	oFile << " - cytosine-adenine : ";
	oFile << setw(6) << (float)caCount / bigramSum * 100 << "%\n";
	oFile << " - cytosine-cytosine: ";
	oFile << setw(6) << (float)ccCount / bigramSum * 100 << "%\n";
	oFile << " - cytosine-guanine : ";
	oFile << setw(6) << (float)cgCount / bigramSum * 100 << "%\n";
	oFile << " - cytosine-thymine : ";
	oFile << setw(6) << (float)ctCount / bigramSum * 100 << "%\n";

	oFile << " - guanine-adenine  : ";
	oFile << setw(6) << (float)gaCount / bigramSum * 100 << "%\n";
	oFile << " - guanine-cytosine : ";
	oFile << setw(6) << (float)gcCount / bigramSum * 100 << "%\n";;
	oFile << " - guanine-guanine  : ";
	oFile << setw(6) << (float)ggCount / bigramSum * 100 << "%\n";
	oFile << " - guanine-thymine  : ";
	oFile << setw(6) << (float)gtCount / bigramSum * 100 << "%\n";

	oFile << " - thymine-adenine  : ";
	oFile << setw(6) << (float)taCount / bigramSum * 100 << "%\n";
	oFile << " - thymine-cytosine : ";
	oFile << setw(6) << (float)tcCount / bigramSum * 100 << "%\n";
	oFile << " - thymine-guanine  : ";
	oFile << setw(6) << (float)tgCount / bigramSum * 100 << "%\n";
	oFile << " - thymine-thymine  : ";
	oFile << setw(6) << (float)ttCount / bigramSum * 100 << "%\n";

	//resetting output to standard
	oFile.unsetf(ios::fixed);
	oFile << noshowpoint << setprecision(6);

	oFile << endl;
}

bool Statistics::CheckLineForInvalidNucleotides(string myString)
{
	int length = myString.length();
	if(length == 0)
	{
		return false;
	}

	for(int i = 0; i < length -1; ++i)
	{

		if(toupper(myString[i]) != 'A' && toupper(myString[i]) != 'C'
		 &&toupper(myString[i]) != 'G' && toupper(myString[i]) != 'T')
		{
			return false;
		}
	}
	return true;
}

void Statistics::Append1000Strings(ofstream& oFile, float standardDev, float mean)
{
	if(standardDev == 0 || mean == 0)
	{
		oFile << "DNA Strings could not be added:" << endl <<
				 " - Normal distribution could not be calculated." << endl;

	}

	float randomValueA;
	float randomValueB;
	float c;
	float d;
	//int randomNucleotideNum;
	float randomNucleotideNum;


	oFile << "1000 DNA Strings\n";
	for(int i = 0; i < 1000; ++i)
	{
		//getting random values between [0, 1)
		randomValueA = (float)rand() / RAND_MAX;
		randomValueB = (float)rand() / RAND_MAX;

		//Calculating C and D
		c = sqrt((-2) * (log(randomValueA))) * cos((2) * M_PI * randomValueB);
		d = standardDev * c + mean;

		for(int j = 0; j < d; ++j)
		{
			//reating random percet value
			randomNucleotideNum = (rand() % 10000);

			randomNucleotideNum = randomNucleotideNum / 100;

			//creating an proportinal number of each nucleotide based on its own percentage
			//if the random value is between 0 and aProb, output A
			if(randomNucleotideNum < aProb)
			{
				oFile << 'A';
			}
			//If the random number is between aProb and aProb + cProb output c
			else if(randomNucleotideNum < aProb + cProb)
			{
				oFile << 'C';
			}
			else if(randomNucleotideNum < aProb + cProb + gProb)
			{
				oFile << 'G';
			}
			else if(randomNucleotideNum < aProb + cProb + gProb + tProb)
			{
				oFile << 'T';
			}
			else
			{
				oFile << "Error";
			}






			/*
			randomNucleotideNum = (rand() % 4) + 1;

			switch(randomNucleotideNum)
			{
				case 1:
					oFile << 'A';
					break;
				case 2:
					oFile << 'C';
					break;
				case 3:
					oFile << 'G';
					break;
				case 4:
					oFile << 'T';
					break;
				default :
					break;
			}//switch(randomNucleotideNum)
			*/
		}//end for(int j = 0; j < d; ++j)
		oFile << endl;
	}//for(int i = 0; i < 1000; ++i)
}//end function
