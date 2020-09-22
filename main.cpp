/*******************************
 * Andrew Dillon
 * 2382400
 * CPSC 350
 * Assignment 01
 *******************************/


#include "Statistics.h"


int main(int argc, char* argv[])
{

	string fileName;
	char continueChar = 'Y';

	ifstream inFile;
	ofstream oFile;
	int fileNum = 1;

	Statistics myStatistics;

	if(argc < 2)
	{
		cout << "Please enter a file name: " ;
		getline(cin, fileName);
	}
	else
	{
		fileName = argv[1];
	}

	//Creating header in file
	oFile.open("AndrewDillon.out");
	oFile << "/*******************************\n"
			 " * Andrew Dillon\n"
			 " * 2382400\n"
			 " * CPSC 350\n"
			 " * Assignment 01\n"
			 " *******************************/\n\n";






	do
	{
		oFile << "File #" << fileNum << ": " << fileName << endl;

		//Process File
		myStatistics.ProcessFile(fileName, oFile);



		//Ask if user wants to use another list
		cout << "Would you like to process another list? (y/n): ";

		do
		{

			cin.get(continueChar);
			cin.ignore(10000000,'\n');
			continueChar = toupper(continueChar);


			if(continueChar != 'Y' && continueChar != 'N')
			{
				cout << "Sorry, please enter either a /'y/' to continue or an /'n/' to exit: ";
			}
		}while(continueChar != 'Y' && continueChar != 'N');

		if(continueChar == 'Y')
		{
			cout << "Please enter the name of the next file: ";
			getline(cin, fileName);
			++fileNum;
		}

	}while(continueChar == 'Y');

	oFile << "END OF FILE\n";

	oFile.close();

	cout << "\n\n\nThank you for using my program.  Have a good day\n\n\n";




	return 0;
}
