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
ofstream fsOut;
ifstream fsIn;
char inLine[] = "";
char outputLine[101];

using namespace std;





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
//we make it 0.2 initially, but you can try different values to see the effects.

double prSpMoveOn = 0.8;
//The probability of moving from the current cognitive state to other states
//	as the next state
//we make it 0.8 initially, but you can try different values to see the effects.

//********************************************************
//Note that prSpRepeat + prSpMoveon should always equal 1
//********************************************************

double spDegenerateTransitionDistancePower = 2;
//The likelihood of moving from the cognitive state of typing some character in a word 
//to the next cognitive state is proportion to the inverse of 
//(spDegenerateTransitionDistancePower) to the 
//(the distance between the current state to the next state)th power.
//In the setting of the original spelling model in the project,
//we make it 2, but you can try different values to see the effects.

double spDegenerateInitialStateDistancePower = 2;
//The likelihood of some character in a word as the initial cognitive state
//is proportion to the inverse of 
//(spDegenerateInitialStateDistancePower) to the 
//(the position of the character in the word)th power.
//In the setting of the original spelling model in the project,
// spDegenerateInitialStateDistancePower and spDegenerateTransitionDistancePower
//have the same value, but you can make them different to see the effects
//By default, we make it 2, but you can try different values to see the effects.


void setParametersSpellingModel()
{
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

void displayParametersSpellingModel()
{
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
//we make it 0.9, but you can try different values to see the effects.

double prKbMiss = 0.4;
//The sum of probabilities that you want to type a character but end in touching 
//a different character.
//we make it 0.1, but you can try different values to see the effects.

//*******************************************************
//Note that prKbHit + prKbMiss should always equal 1
//*******************************************************

//In the setting of the original keyboard model in the project,
//we make it 0.2, but you can try different values to see the effects.


double kbDegenerateDistancePower = 2;
//The likelihood you want to type a character but end in touching a different character
//is proportion to the inverse of 
//(kbDegenerateDistancePower) raised to the (distance between them) th power
//In the setting of the original keyboard model in the project,
//we make it 2, but you can try different constants to see the effects.

double scaleFactor = 0;


void setParametersKbModel()
{
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

void displayParametersKbModel()
{
	cout << endl
		<< "Parameter values of the keyboard model:" << endl
		<< "   P_hit  = " << prKbHit << endl
		<< "   P_miss = " << prKbMiss << endl
		<< "   deg_kb = " << kbDegenerateDistancePower << endl << endl;
}


//******************************************************************************
//******************************************************************************
//
//  Trace flag and the alphabet table
//
//******************************************************************************
//******************************************************************************
bool traceON = true;   // output tracing messages


					   /************************************************************************/
					   //Calculate and return the probability of charGenerated actually typed
					   //given charOfTheState as the underlying cognitive state. 
					   /************************************************************************/
double prCharGivenCharOfState(char charGenerated, char charOfTheState, double scaleFactor)
//1/2 + 1/4 + 1/8
{   // KEYBOARD STATE
	// CharGenerated = What we actually touched (typed)
	// CharOfTheState = What we want to type in our mind (our cognitive state)
	if (charGenerated == charOfTheState) return prKbHit;
	else {
		double bottom = pow(kbDegenerateDistancePower, distanceOfLetters(charGenerated, charOfTheState));
		double theNum = (prKbMiss / bottom);
		return  theNum / scaleFactor;
	}
}

int distanceOfLetters(char charOne, char charTwo) {
	int dist = abs(charTwo - charOne);
	if (dist > 13) {
		return abs(dist - 26);
	}
	return dist;
}


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


/************************************************************************/
//Calculate for each cognitive state excluding the special states I and F,
//     the probability of that cognitive state being the first cognitive state
//     after the transition out of the special state I.
//Note that the value of the parameter sizeOfTable should be
//     exactly the number characters in the word.
//Store these prbabilities in the prTable array.
/************************************************************************/
void getPrTableForPossibleInitialStates(double prTable[], int sizeOfTable)
{
	//It is a word of sizeOfTable characters:
	//     i.e. the number of character states is sizeOfTable.
	//     let's index these characters from 0 to sizeOfTable-1.
	//

	//First calculate the sum of ratios of probabilities of
	//     going from the special I state into these character states.
	//This allows you to calculate the scaling factor to determine the probabilities.


	//Second, for each character state calculate the probability
	//     transitioning from the special I state into the character state.
	double scale = 0;
	for (int i = 0; i < sizeOfTable; i++) {
		scale += pow(2, i + 1);
		//finds scale by adding sum of each element in array (state)
	}
	for (int j = 0; j < sizeOfTable; j++) {
		prTable[abs(j - sizeOfTable) - 1] = pow(2, j + 1) / scale;
	}
}




/************************************************************************/
//Calculate for each actual cognitive state for a word
//     (excluding the special I state),
//     the probability of that cognitive state being the next cognitive state
//     given that currentState is the index of the current state.
//Note that the value of the parameter sizeOfTable should be
//     1 + the number characters in the word,
//Store these probabilities in the transitionPrTable[] array.
/************************************************************************/
void getPrTableForPossibleNextStates
(double transitionPrTable[], int sizeOfTable, int currentState)
{
	//We are working on a word of sizeOfTable-1 characters:
	//     i.e. the number of character states is sizeOfTable-1.
	//Index these character states from 0 to sizeOfTable-2 respectively, while
	//     the index of the special final state F is sizeOfTable-1.
	//currentState is the index of the current state in the word

	//First calculate the sum of ratios of probabilities of
	//     going from the current state into the other down-stream states down in word
	//     including all the down-stream character states and the
	//     special F final state.
	//This allows you to calculate the scaling factor to determine the probabilities.

	//Second, for each state (excluding the special I state)
	//     calculate the probability of
	//     transitioning from the current state into the character state
	//     and store the probability into the table.

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

/************************************************************************/
//  Programming #2 
//
//  List below are function prototypes of functions given  or required 
//		to be implemented in Programming #2 
//
/************************************************************************/

/************************************************************************/
//Given the probabilities (of sizeOfTable elements) stored in prTable,
//	try to randomly take a sample out of sizeOfTable elements
//	according to the probabilities of sizeOfTable elements;
//Return the index of the non-deterministically sampled element.
/************************************************************************/
int take1SampleFrom1PrSpace(double prTable[], int sizeOfTable)
{
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





/************************************************************************/
//
//Given the character to type (charToType) 
//	(assuming that the 1D keyboard of 26 keys is used),
//	(assuming that prTable[] for storing 26 prbabilities),
//	calculate pr(charGenerated = 'a' | charToType),
//	calculate pr(charGenerated = 'b' | charToType), 
//	...
//	calculate pr(charGenerated = 'z' | charToType), and
//	store these probabilities in prTable.
/************************************************************************/
void getKeyboardProbabilityTable(char charToType, double prTable[])
{
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


/************************************************************************/
//Simulate the keyboard model:
//Given the charToTye, simulate what character may actually
//	be typed and return it as the result.
/************************************************************************/
char typeOneChar(char charToType)
{
	char asciiMap[] = "abcdefghijklmnopqrstuvwxyz";
	double prTable[26];
	getKeyboardProbabilityTable(charToType, prTable);
	return asciiMap[take1SampleFrom1PrSpace(prTable, 26)];
}



/************************************************************************/
//Simulate the combination of the spelling model and the keyboard model:
//Given a word stored in the word array, simulate what may actually
//	be typed and store the result in the output array.
//maxOutput specifies the capacity limit of the output array, by default it is 100.
//When traceON is true (by default it is false), extra outputs are provided as traces.
/************************************************************************/
void typeOneWord(char word[], char output[], bool traceON, int maxOutput)
{
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
}// end of the function
 /*******************************************************************/
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
		cout << i + 1 << " Column: " << endl;
		for (int j = 0; j < sizeOfTableStates; j++) {
			cout << j + 1 << " Row: " << endl;
			sum = 0;
			for (int k = 0; k < sizeOfTableStates; k++) {
				if (k == 0) {
					double* initialStates = new double[sizeOfTableStates];
					getPrTableForPossibleInitialStates(initialStates + 1, sizeOfStates);
					initialStates[sizeOfTableStates - 1] = 0;
					initialStates[0] = 0;
					cout << fowardTable[i - 1][k] << " * " << initialStates[j] << " = " << fowardTable[i - 1][k] * initialStates[j] << endl;
					sum += fowardTable[i - 1][k] * initialStates[j];
				}
				else {
					double* currentStates = new double[sizeOfTableStates];
					getPrTableForPossibleNextStates(currentStates, sizeOfStates + 1, k - 1);
					if (j < 1) {
						//Initial State
						cout << fowardTable[i - 1][k] << " * " << 0 << " = 0" << endl;
						sum += fowardTable[i - 1][k] * 0;
					}
					else {
						cout << fowardTable[i - 1][k] << " * " << currentStates[j - 1] << " = " << fowardTable[i - 1][k] * currentStates[j - 1] << endl;
						sum += fowardTable[i - 1][k] * currentStates[j - 1];
					}
				}
				
			}
			double prOfTyping = prCharGivenCharOfState(observedString[i - 1], wordString[j - 1], scaleFactor);
			//prOfTyping is giving wrong calculations
			if (j > sizeOfStates || j < 1 || i > sizeOfObservedString)
			{
				prOfTyping = 0;
				if (j > sizeOfStates && i > sizeOfObservedString)
					prOfTyping = 1;
			}
			cout << "total is " << sum << " * " << prOfTyping << " = ";
			sum *= prOfTyping;
			cout << sum << endl << endl;
			fowardTable[i][j] = sum;
		}
	}
	return fowardTable[sizeOfTableOberservedString - 1][sizeOfTableStates - 1];
}

