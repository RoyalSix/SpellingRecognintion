#include <cstdlib>
#include <ctime>
#include <cstring>
#include <iostream>
#include <fstream>
#include <conio.h>
#include <array>
#include <vector>
#include "models.h"
#include <fstream>
#include <math.h> 
ofstream fsOut;
ifstream fsIn;
char inLine[] = "";
char outputLine[101];

using namespace std;



//TODO: FINISH COMMENTING

//******************************************************************************
//******************************************************************************
//
//  Parameter setting for the storage capacity
//
//******************************************************************************
//******************************************************************************

//The maximum capacity (maximum number of characters allowed) 
//		for the storing the vocabulary set
#define VocabularyStorageLimit 20000

//The maximum capacity (maximum number of characters allowed) 
//		for storing the corrupted words during spelling recognition
#define NoisyWordsStorageLimit 15000





//******************************************************************************
//******************************************************************************
//
//  Parameter settings for the Spelling model
//
//******************************************************************************
//******************************************************************************
double prSpRepeat = 0.2;
//The probability of repeating the current cognitive state again as the next state


double prSpMoveOn = 0.8;
//The probability of moving from the current cognitive state to other states
//	as the next state

//********************************************************
//prSpRepeat + prSpMoveon should always equal 1
//********************************************************

double spDegenerateTransitionDistancePower = 2;
//The likelihood of moving from the cognitive state of typing some character in a word 
//to the next cognitive state is proportion to the inverse of 
//(spDegenerateTransitionDistancePower) to the 
//(the distance between the current state to the next state)th power.
//In the setting of the original spelling model in the project,


double spDegenerateInitialStateDistancePower = 2;
//The likelihood of some character in a word as the initial cognitive state
//is proportion to the inverse of 
//(spDegenerateInitialStateDistancePower) to the 
//(the position of the character in the word)th power.
//In the setting of the original spelling model in the project,
// spDegenerateInitialStateDistancePower and spDegenerateTransitionDistancePower
//have the same value, but you can make them different to see the effects



/**
  * @desc Lets user set the spelling model parameters, these control
  *		the probability that an agent will type the same word twice
  * 	or move on.
*/
void setParametersSpellingModel() {
	cout << endl
		<< "Reset the parameters of the spelling model:" << endl << endl;

	cout << "Reset P_moveOn, the probability of moving on" << endl
		<< "   from the current cognitive state to other states:" << endl
		<< "P_moveOn = ";
	cin >> prSpMoveOn;

	prSpRepeat = 1 - prSpMoveOn;
	cout << endl
		<< "Automatically reset P_repeat to 1-P_moveOn" << endl
		<< "P_repeat = " << prSpRepeat << endl;

	cout << endl
		<< "Do you want to change the deg_sp? (y or n)" << endl;

	if (_getch() == 'y')
	{
		cout << "Reset deg_sp, the power coefficient governing the probability of " << endl
			<< "   skipping from the current cognitive state to other states:" << endl
			<< "deg_sp = ";

		cin >> spDegenerateTransitionDistancePower;

		spDegenerateInitialStateDistancePower = spDegenerateTransitionDistancePower;
	}
}


void displayParametersSpellingModel() {
	cout << endl
		<< "Parameter values of the spelling model:" << endl
		<< "   P_repeat  = " << prSpRepeat << endl
		<< "   P_moveOn  = " << prSpMoveOn << endl
		<< "   deg_sp = " << spDegenerateTransitionDistancePower << endl << endl;
}

//******************************************************************************
//******************************************************************************
//
//  Parameter settings for the keyboard model
//
//******************************************************************************
//******************************************************************************

double prKbHit = 0.6;
//The probability that you want to type a character and you do successfully make it
//In the setting of the original keyboard model in the project,

double prKbMiss = 0.4;
//The sum of probabilities that you want to type a character but end in touching 
//a different character.

//*******************************************************
//prKbHit + prKbMiss should always equal 1
//*******************************************************



double kbDegenerateDistancePower = 2;
//The likelihood you want to type a character but end in touching a different character
//is proportion to the inverse of 
//(kbDegenerateDistancePower) raised to the (distance between them) th power
//In the setting of the original keyboard model in the project,

double scaleFactor = 0;



/**
  * @desc Lets user set the keyboard model parameters, these control
  *		the probability that an agent will try to type a character
  * 	in their mind and they will hit it.
*/
void setParametersKbModel() {
	cout << endl
		<< "Reset the parameters of the keyboard model:" << endl << endl;

	cout << "Reset P_hit, the probability of hitting" << endl
		<< "   the right character wanted on the keyboard:" << endl
		<< "P_hit = ";
	cin >> prKbHit;

	prKbMiss = 1 - prKbHit;
	cout << endl
		<< "Automatically reset P_miss to 1-P_hit" << endl
		<< "P_miss = " << prKbMiss << endl;

	cout << endl
		<< "Do you want to change the deg_kb? (y or n)" << endl;

	if (_getch() == 'y')
	{
		cout << "Reset deg_kb, the power coefficient governing the probability of " << endl
			<< "   skipping from the current cognitive state to other states:" << endl
			<< "deg_kb = ";

		cin >> kbDegenerateDistancePower;
	}
}

void displayParametersKbModel(){
	cout << endl
		<< "Parameter values of the keyboard model:" << endl
		<< "   P_hit  = " << prKbHit << endl
		<< "   P_miss = " << prKbMiss << endl
		<< "   deg_kb = " << kbDegenerateDistancePower << endl << endl;
}

bool traceON = true;   // output tracing of state transitions

/**
  * @desc Calculate and return the probability of charGenerated actually typed
  *		given charOfTheState as the underlying cognitive state. 
  * @param char charGenerated - What we actually touched (typed)
  * @param char charOfTheState - What we want to type in our mind (our cognitive state)
  * @param double scaleFactor - the pre computated scale factor
  * @return double - Pr(charOfTheState | charGenerated)
*/
double prCharGivenCharOfState(char charGenerated, char charOfTheState, double scaleFactor) {
	if (charGenerated == charOfTheState) return prKbHit;
	else {
		double bottom = pow(kbDegenerateDistancePower, distanceOfLetters(charGenerated, charOfTheState));
		double theNum = (prKbMiss / bottom);
		return  theNum / scaleFactor;
	}
}


/**
  * @desc find the distance of two observed letters based
  * 	on a keyboard model that loops around in a 360
  * @param string $msg - the message to be displayed
  * @return int dist - the distance of two given 
*/
int distanceOfLetters(char charOne, char charTwo) {
	int dist = abs(charTwo - charOne);
	if (dist > 13) {
		return abs(dist - 26);
	}
	return dist;
}

/**
* @desc - determines the scale factor for the model based on the
* 	amount of possible observations
* @param int size - size of the table
* @return double scale - the computed scale
*/
double getScaleFactor(int size) {
	double scale = 0;
	for (int i = 1; i < size; i++) {
		if (i > size / 2) {
			int j = abs(i - size);
			scale = scale + (1.0 / pow(2.0, j));
		}
		else {
			scale = scale + (1.0 / pow(2.0, i));
		}
	}
	return scale;
}

/**
* @desc - Calculate for each cognitive state excluding the special states I and F,
*	the probability of that cognitive state being the first cognitive state
*	after the transition out of the special state I.
* @param double * prTable - The transition table to
* 	store state Pr for inital states based on amount of possible states
* @param int sizeOfTable - size of the table
*/
void getPrTableForPossibleInitialStates(double prTable[], int sizeOfTable) {
	//It is a word of sizeOfTable characters:
	//     i.e. the number of character states is sizeOfTable.
	//     these characters are indexed from 0 to sizeOfTable-1.
	//

	double scale = 0;
	for (int i = 0; i < sizeOfTable; i++) {
		scale += pow(2, i + 1);
		//finds scale by adding sum of each element in array (state)
	}
	for (int j = 0; j < sizeOfTable; j++) {
		prTable[abs(j - sizeOfTable) - 1] = pow(2, j + 1) / scale;
	}
}




/**
* @desc Calculate for each actual cognitive state for a word
* 	(excluding the special I state), the probability of 
* 	that cognitive state being the next cognitive state
* 	given that currentState is the index of the current state.
* 	The value of the parameter sizeOfTable should be
* 	1 + the number characters in the word,
* @param double * transitionPrTable - The transition table to
*	store state Pr for next possible states (not including
* 	initial states)
* @param int sizeOfTable - size of the table
* @param int currentState - current state index before the transition
*/
void getPrTableForPossibleNextStates (double transitionPrTable[], int sizeOfTable, int currentState) {
	double scale = 0;
	for (int i = 0; i < sizeOfTable - 1 - currentState; i++) {
		scale += pow(2, i + 1);
	}

	for (int j = 0; j <= sizeOfTable; j++) {
		if (abs(j - sizeOfTable) != sizeOfTable) {
			transitionPrTable[abs(j - sizeOfTable)] = (pow(2, j) / scale) * prSpMoveOn;
			if (abs(j - sizeOfTable) < currentState) {
				transitionPrTable[abs(j - sizeOfTable)] = 0;
			}
			if (abs(j - sizeOfTable) == currentState) {
				transitionPrTable[abs(j - sizeOfTable)] = prSpRepeat;
			}
		}
	}

}

/**
* @desc Given the probabilities (of sizeOfTable elements) stored in prTable,
*	randomly takes a sample out of sizeOfTable elements
*	according to the probabilities of sizeOfTable elements;
* @param double * prTable - Pr table that includes Pr for all observations
* 	in model regardless of current state
* @param int sizeOfTable - the size of the HMM table
* @return int - the index of the non-deterministically sampled element.
*/
int take1SampleFrom1PrSpace(double prTable[], int sizeOfTable) {
	int i;
	double prSum = 0;
	for (i = 0; i<sizeOfTable; i++)
		prSum += prTable[i];
	if (prSum < 0.999 || prSum > 1.001)
		cout << "Something is wrong with a random sampleing" << endl
		<< "The sum of all probabilities in the table is " << prSum << endl;

	//Calculate the probability intervals of the characters
	double prAccumulated = 0;
	double * prIntervals;
	prIntervals = new double[sizeOfTable];
	for (i = 0; i<sizeOfTable; i++)
	{
		prAccumulated = prAccumulated + prTable[i];
		prIntervals[i] = prAccumulated;
	}

	/*
	cout << endl;
	for (i=0; i<sizeOfTable; i++)
	{ if (i>0)
	cout << "state " << i << " lower bounded by " << prIntervals[i-1]<< endl;
	cout << "state " << i << " upper bounded by " << prIntervals[i]<< endl;
	}
	*/

	// Generate a random number in [0,1]
	i = rand() % 1001;
	double temp = i / 1000.0;
	//cout << "The random number pair generated is " << i << endl << temp << endl;

	bool sampleTaken = false;
	for (i = 0; i < sizeOfTable && !sampleTaken; i++)
		if (temp <= prIntervals[i])
		{
			delete[] prIntervals;
			//cout << "The random interval id is " << i << endl;
			return i;
		}
	return sizeOfTable - 1;

}




/**
  * @desc Given the character to type (charToType) calculate 
  * 	pr(charGenerated = a observation | charToType) for 
  * 	each observation, and store these probabilities in prTable.
  * @param char charToType - character in the state the agent is
  * 	trying to type
  * @param double prTable - The transition table to
*	store observation Pr's in account of current charToType
*/
void getKeyboardProbabilityTable(char charToType, double prTable[]) {
	int currentState = 0;
	int sizeOfTable = 26;
	char asciiMap[] = "abcdefghijklmnopqrstuvwxyz";
	double scaleFactor = getScaleFactor(26);
	for (int i = 0; i < sizeOfTable; i++) {
		if (prTable[i] == charToType) {
			currentState = i;
		}
	}
	for (int i = 0; i < 26; i++) {
		prTable[i] = prCharGivenCharOfState(asciiMap[i], charToType, scaleFactor);
	}
}

/**
  * @desc Given the charToType, simulates what character may actually
  * 	be typed and return it as the result.
  * @param char charToType - The character in the state being attempted
  *		to type
  * @return char - the randomly chosen character based on the spelling
  * 	model and the keyboard model
*/
char typeOneChar(char charToType) {
	char asciiMap[] = "abcdefghijklmnopqrstuvwxyz";
	double prTable[26];
	getKeyboardProbabilityTable(charToType, prTable);
	return asciiMap[take1SampleFrom1PrSpace(prTable, 26)];
}


/**
  * @desc Given a word stored in the word array, simulates what may actually
  * be typed and store the result in the output array.
*/
void typeOneWord(char word[], char output[], bool traceON, int maxOutput) {
	int size = 0;
	int currentState = 0;
	char charInState;
	for (int max = 0; max < 50; max++) {
		output[max] = '\0';
	}
	while (word[size]) {
		size++;
	}
	double* wordArray = new double[size];
	getPrTableForPossibleInitialStates(wordArray, size);
	currentState = take1SampleFrom1PrSpace(wordArray, size);
	charInState = word[currentState];
	output[0] = typeOneChar(charInState);
	if (traceON) cout << "'" << output[0] << "'" << " pressed, current state: "
		<< word[currentState] << ", " << currentState;
	if (traceON) {
		cout << " the next state: ";
	}
	delete[] wordArray;

	//main loop
	size++;
	bool finalState = false;
	for (int i = 1; !finalState; i++) {
		double* wordArray = new double[size];
		getPrTableForPossibleNextStates(wordArray, size, currentState);
		currentState = take1SampleFrom1PrSpace(wordArray, size);
		if (currentState >= size - 1) {
			if (traceON) {
				finalState = true;
				cout << "final state." << endl << endl;
			}
			break;

		}
		cout << word[currentState] << ", " << currentState << endl;
		charInState = word[currentState];
		output[i] = typeOneChar(charInState);
		if (traceON) {
			cout << "'" << output[i] << "'" << " pressed, current state: " << word[currentState] << ", " << currentState << " the next state: ";
		}
		delete[] wordArray;
	}
}


void typeOneArticle(char * corruptedMessageFile, char * sourceArticle, bool trace) {
	fsIn.open(sourceArticle);
	fsOut.open(corruptedMessageFile);
	while (fsIn >> inLine) {
		typeOneWord(inLine, outputLine, true);
		fsOut << outputLine << endl;
	}
	fsOut.close();
	fsIn.close();
}


double prOf1CharSeriesWhenTyping1Word(char observedString[], char wordString[]) {
	double scaleFactor = getScaleFactor(26);
	int sizeOfStates = strlen(wordString);
	int sizeOfObservedString = strlen(observedString);
	int sizeOfTableStates = sizeOfStates + 2;
	int sizeOfTableOberservedString = sizeOfObservedString + 2;

	double sum = 0;
	vector < vector <double> > fowardTable(sizeOfTableOberservedString, vector <double>(sizeOfTableStates, 0));
	fowardTable[0][0] = 1;
	for (int i = 1; i < sizeOfTableOberservedString; i++) {
		//cout << i + 1 << " Column: " << endl;
		for (int j = 0; j < sizeOfTableStates; j++) {
			//cout << j + 1 << " Row: " << endl;
			sum = 0;
			for (int k = 0; k < sizeOfTableStates; k++) {
				if (k == 0) {
					double* initialStates = new double[sizeOfTableStates];
					getPrTableForPossibleInitialStates(initialStates + 1, sizeOfStates);
					initialStates[sizeOfTableStates - 1] = 0;
					initialStates[0] = 0;
					//cout << fowardTable[i - 1][k] << " * " << initialStates[j] << " = " << fowardTable[i - 1][k] * initialStates[j] << endl;
					sum += fowardTable[i - 1][k] * initialStates[j];
					delete[] initialStates;
				}
				else {
					double* currentStates = new double[sizeOfTableStates];
					getPrTableForPossibleNextStates(currentStates, sizeOfStates + 1, k - 1);
					if (j < 1) {
						//Initial State
						//cout << fowardTable[i - 1][k] << " * " << 0 << " = 0" << endl;
						sum += fowardTable[i - 1][k] * 0;
					}
					else {
						//cout << fowardTable[i - 1][k] << " * " << currentStates[j - 1] << " = " << fowardTable[i - 1][k] * currentStates[j - 1] << endl;
						sum += fowardTable[i - 1][k] * currentStates[j - 1];
					}
					delete[] currentStates;
				}

			}
			double prOfTyping = prCharGivenCharOfState(observedString[i - 1], wordString[j - 1], scaleFactor);
			if (j > sizeOfStates || j < 1 || i > sizeOfObservedString)
			{
				prOfTyping = 0;
				if (j > sizeOfStates && i > sizeOfObservedString)
					prOfTyping = 1;
			}
			//cout << "total is " << sum << " * " << prOfTyping << " = ";
			sum *= prOfTyping;
			//cout << sum << endl << endl;
			fowardTable[i][j] = sum;
		}
	}
	return fowardTable[sizeOfTableOberservedString - 1][sizeOfTableStates - 1];
}


void testPrOfTypingDocument() {
	bool cont = true;
	while (cont) {
		string doc1;
		char person;
		cout << "Enter the letter of the person to test: " << endl;
		cout << "W.		Winnie" << endl;
		cout << "M.		Manny" << endl;
		cout << "C.		Cathy" << endl;
		cout << "J.		Johnny" << endl;
		cin >> person;
		switch (person) {
		case 'W': case 'w':
			prKbHit = 0.7;
			prKbMiss = 0.3;
			kbDegenerateDistancePower = 2;
			spDegenerateTransitionDistancePower = 2;
			prSpRepeat = 0.1;
			prSpMoveOn = 0.9;

			break;
		case 'M': case 'm':
			prKbHit = 0.9;
			prKbMiss = 0.1;
			kbDegenerateDistancePower = 2;
			spDegenerateTransitionDistancePower = 2;
			prSpRepeat = 0.3;
			prSpMoveOn = 0.7;
			break;
		case 'C': case 'c':
			prKbHit = 0.7;
			prKbMiss = 0.3;
			kbDegenerateDistancePower = 2;
			spDegenerateTransitionDistancePower = 2;
			prSpRepeat = 0.3;
			prSpMoveOn = 0.7;
			break;
		case 'J': case 'j':
			prKbHit = 0.9;
			prKbMiss = 0.1;
			kbDegenerateDistancePower = 2;
			spDegenerateTransitionDistancePower = 2;
			prSpRepeat = 0.1;
			prSpMoveOn = 0.9;
			break;
		}
		cout << "Got the parameters: " << endl;
		displayParametersKbModel();
		displayParametersSpellingModel();
		cout << "Enter the letter of the document you want to test: " << endl;
		cin >> doc1;
		cout << "The probability of seeing this document when trying to type the biola statement is: " <<
			logPrOfGettingDocument1WhenTypingDocument2("biola.txt", doc1 + ".txt") << endl;
		cout << "Hit any key to go again or press '.' to end." << endl;
		if (_getch() == '.') cont = false;
	}
 }


double logPrOfGettingDocument1WhenTypingDocument2(string document1, string document2) {
	double sum = 0;
	ifstream fsInDoc1;
	ifstream fsInDoc2;
	string inLineDoc1;
	string inLineDoc2;
	fsInDoc2.open(document2);
	fsInDoc1.open(document1);
	while (fsInDoc2 >> inLineDoc2 && fsInDoc1 >> inLineDoc1) {
		char* doc2Char = new char[inLineDoc2.length() + 1];
		strcpy(doc2Char, inLineDoc2.c_str());
		char* doc1Char = new char[inLineDoc1.length() + 1];
		strcpy(doc1Char, inLineDoc1.c_str());
		//cout << "Pr(" << inLineDoc1 << "|" << inLineDoc2 << ") = " <<
			//prOf1CharSeriesWhenTyping1Word(doc2Char, doc1Char) << " + " << sum << " +" << endl;
		sum += log(prOf1CharSeriesWhenTyping1Word(doc2Char, doc1Char));
		delete[] doc1Char;
		delete[] doc2Char;
	}
	return sum;
}


void getPrFromDocument(string document1, string document2) {
	const double INITIAL_INC = 0.1;
	prKbHit = INITIAL_INC;
	prKbMiss = 1 - INITIAL_INC;
	prSpRepeat = INITIAL_INC;
	prSpMoveOn = 1 - INITIAL_INC;
	double highestPr = logPrOfGettingDocument1WhenTypingDocument2(document1, document2);;
	double highestPrHit;
	double highestPrRepeat;
	for (double i = INITIAL_INC; i < 1; i += INITIAL_INC) {
		for (double j = INITIAL_INC; j < 1; j += INITIAL_INC) {
			i = floor(i * 100.00 + 0.5) / 100.00;
			j = floor(j * 100.00 + 0.5) / 100.00;
			prKbHit = i;
			prKbMiss = 1 - i;
			prSpRepeat = j;
			prSpMoveOn = 1 - j;
			double currPr = logPrOfGettingDocument1WhenTypingDocument2(document1, document2);
			if (currPr > highestPr) {
				highestPr = currPr;
				highestPrHit = i;
				highestPrRepeat = j;
			}
		}
	}
	prKbHit = highestPrHit;
	prKbMiss = 1 - highestPrHit;
	prSpRepeat = highestPrRepeat;
	prSpMoveOn = 1 - highestPrRepeat;
}

void testgetPrFromDocument() {
	bool cont = true;
	string CORRUPTED = "corruptedBiolaVision1";
	while (cont) {
		cout << "To test " << CORRUPTED << " press . or enter name of document" << endl;
		if (_getch() == '.') {
			cout << "Learning from the document: " << CORRUPTED << endl;
		}
		else {
			cin >> CORRUPTED;
			cout << "Learning from the document: " << CORRUPTED << ".txt" << endl;
		}
		getPrFromDocument("biola.txt", CORRUPTED + ".txt");
		cout << "Got most likely parameters as: "; 
		displayParametersKbModel();
		displayParametersSpellingModel();
		cout << "Hit any key to go again or press '.' to end." << endl;
		if (_getch() == '.') cont = false;
	}
}

void testRecoverOriginal() {
	bool cont = true;
	string doc1 = "vocabulary.txt"; 
	string doc2 = "corruptedMessage1.txt";
	string doc3 = "recoveredMessage_V1.txt";
	cout << "Parameters for this test are: " << endl;
	displayParametersKbModel();
	displayParametersSpellingModel();
	while (cont) {
		cout << "Got document " << doc2 << endl;
		cout << "Attempting to recover message..." << endl;
		recoverOriginal(doc1, doc2, doc3);
		cout << "Message recovered, results in " << doc3 << endl;
		cout << "Hit any key to go again or press '.' to end." << endl;
		if (_getch() == '.') cont = false;
	}
}

void recoverOriginal(string doc1, string doc2, string doc3) {
	const int AMOUNT_OF_DATA = 4;
	ifstream fsInDoc1;
	ifstream fsInDoc2;
	ofstream fsOutDoc;
	string inLineDoc1;
	string inLineDoc2;
	fsInDoc2.open(doc2);
	fsOutDoc.open(doc3);
	while(fsInDoc2 >> inLineDoc2){
		//observed words
		vector<double> highestPr;
		vector<string> highestPrString;
		fsInDoc1.open(doc1);
		while (fsInDoc1 >> inLineDoc1) {
			//vocab words
			char* doc2Char = new char[inLineDoc2.length() + 1];
			strcpy(doc2Char, inLineDoc2.c_str());
			char* doc1Char = new char[inLineDoc1.length() + 1];
			strcpy(doc1Char, inLineDoc1.c_str());
			//cout << "Pr(" << inLineDoc1 << "|" << inLineDoc2 << ") = " <<
			//prOf1CharSeriesWhenTyping1Word(doc2Char, doc1Char) << " + " << sum << " +" << endl;
			double currPr = log(prOf1CharSeriesWhenTyping1Word(doc2Char, doc1Char));
			if (highestPr.size() == 0) {
				highestPr.insert(highestPr.end(), currPr);
				highestPrString.insert(highestPrString.end(), doc1Char);
			}
			else {
				for (int i = 0; i < AMOUNT_OF_DATA && i < highestPr.size(); i++) {
					if (currPr > highestPr[i]) {
						highestPr.insert(highestPr.begin() + i, currPr);
						highestPrString.insert(highestPrString.begin() + i, doc1Char);
						break;
					}
					if (i == highestPr.size() - 1) {
						highestPr.insert(highestPr.end(), currPr);
						highestPrString.insert(highestPrString.end(), doc1Char);
						break;
					}
				}
			}
			delete[] doc1Char;
			delete[] doc2Char;
		}
		for (int i = 0; i < AMOUNT_OF_DATA; i++) {
			fsOutDoc << highestPrString[i] << " ";
		}
		fsOutDoc << endl;
		fsInDoc1.close();
	}
	fsOutDoc.close();
}

void testRecoverOriginal2() {
	bool cont = true;
	string doc1 = "vocabulary.txt";
	string doc2 = "corruptedMessage1.txt";
	string doc3 = "corruptedMessage2.txt";
	string doc4 = "recoveredMessage_V2.txt";
	cout << "Got the parameters: " << endl;
	displayParametersKbModel();
	displayParametersSpellingModel();
	while (cont) {
		cout << "Got document " << doc2 << " and " << doc3 << endl;
		cout << "Attempting to recover message..." << endl;
		recoverOriginal2(doc1, doc2, doc3, doc4);
		cout << "Message recovered, results in " << doc3 << endl;
		cout << "Hit any key to go again or press '.' to end." << endl;
		if (_getch() == '.') cont = false;
	}
}

void recoverOriginal2(string doc1, string doc2, string doc3, string doc4) {
	const int AMOUNT_OF_DATA = 4;
	vector<vector <string> > listOfHighestPrsString;
	ifstream fsInDoc1;
	ifstream fsInDoc2;
	ifstream fsInDoc3;
	ofstream fsOutDoc;
	string inLineDoc1;
	string inLineDoc2;
	string inLineDoc3;
	fsInDoc2.open(doc2);
	fsInDoc3.open(doc3);
	fsOutDoc.open(doc4);
	while (fsInDoc2 >> inLineDoc2 && fsInDoc3 >> inLineDoc3) {
		//observed words
		vector<double> highestPr;
		vector<string> highestPrString;
		fsInDoc1.open(doc1);
		while (fsInDoc1 >> inLineDoc1) {
			//vocab words
			char* doc3Char = new char[inLineDoc3.length() + 1];
			strcpy(doc3Char, inLineDoc3.c_str());
			char* doc2Char = new char[inLineDoc2.length() + 1];
			strcpy(doc2Char, inLineDoc2.c_str());
			char* doc1Char = new char[inLineDoc1.length() + 1];
			strcpy(doc1Char, inLineDoc1.c_str());
			//cout << "Pr(" << inLineDoc1 << "|" << inLineDoc2 << ") = " <<
			//prOf1CharSeriesWhenTyping1Word(doc2Char, doc1Char) << " + " << sum << " +" << endl;
			double currPrDoc2 = log(prOf1CharSeriesWhenTyping1Word(doc2Char, doc1Char));
			double currPrDoc3 = log(prOf1CharSeriesWhenTyping1Word(doc3Char, doc1Char));
			double currPr = currPrDoc2 * currPrDoc3 * - 1;
			if (highestPr.size() == 0) {
				highestPr.insert(highestPr.end(), currPr);
				highestPrString.insert(highestPrString.end(), doc1Char);
			}
			else {
				for (int i = 0; i < AMOUNT_OF_DATA && i < highestPr.size(); i++) {
					if (currPr > highestPr[i]) {
						highestPr.insert(highestPr.begin() + i, currPr);
						highestPrString.insert(highestPrString.begin() + i, doc1Char);
						break;
					}
					if (i == highestPr.size() - 1) {
						highestPr.insert(highestPr.end(), currPr);
						highestPrString.insert(highestPrString.end(), doc1Char);
						break;
					}
				}
			}
			delete[] doc1Char;
			delete[] doc2Char;
			delete[] doc3Char;
		}
		for (int i = 0; i < AMOUNT_OF_DATA; i++) {
			fsOutDoc << highestPrString[i] << " ";
		}
		fsOutDoc << endl;
		fsInDoc1.close();
	}
	fsOutDoc.close();
}

void testPrecisionOfRecovered() {
	bool cont = true;
	string doc1 = "message.txt";
	string doc2 = "recoveredMessage_V1.txt";
	string doc3 = "recoveredMessage_V2.txt";
	while (cont) {
		cout << "To test " << doc2 << " and " << doc3 << " press . or enter name of document" << endl;
		if (_getch() == '.') {
			cout << "Getting precision from the document: " << doc3 << endl;
			cout << "Got precision of " << doc2 << ":" << endl << getPrecisionOfRecovered(doc2, doc1) * 100 << " percent." << endl;
			cout << "Got precision of: " << doc3 << ":" << endl << getPrecisionOfRecovered(doc3, doc1) * 100 << " percent." << endl;
		}
		else {
			cin >> doc2;
			cout << "Getting precision from the document: " << doc2 << ".txt" << endl;
		}
		cout << "Hit any key to go again or press '.' to end." << endl;
		if (_getch() == '.') cont = false;
	}
}

double getPrecisionOfRecovered(string observedDoc, string originalDoc) {
	const double TOTAL_TRYS_GIVEN = 4;
	int totalLines = 0;
	int misses = 0;
	ifstream fsInDoc1;
	ifstream fsInDoc2;
	string inLineDoc1;
	string inLineDoc2;
	fsInDoc1.open(observedDoc);
	fsInDoc2.open(originalDoc);
	while (fsInDoc2 >> inLineDoc2 ) {
		totalLines++;
		int lineCounter = 0;
		while (lineCounter < 4 && fsInDoc1 >> inLineDoc1) {
			if (inLineDoc1 == inLineDoc2) {
				getline(fsInDoc1, inLineDoc1);
				break;
			}
			lineCounter++;
			misses++;
		}
	}
	return 1 - (misses/ (totalLines * TOTAL_TRYS_GIVEN));
}
