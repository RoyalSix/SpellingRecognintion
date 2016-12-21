
#include <string>
#define ALPHABET "abcdefghijklmnopqrstuvwxyz"
#define ALPHABET_SIZE 26

using namespace std;

double prCharGivenCharOfState(char charGenerated, char charOfTheState, double scaleFactor);
/************************************************************************/
//Calculate and return the probability of charGenerated actually typed
//given charOfTheState as the underlying cognitive state. 
/************************************************************************/

void getPrTableForPossibleInitialStates(double prTable[], int sizeOfTable);
/************************************************************************/
//Calculate for each actual cognitive state in a word of sizeOfTable characters,
//	the probability of that cognitive state being the actual first cognitive state
//	when the user types the word.
//Store these prbabilities in the prTable array.
/************************************************************************/

void getPrTableForPossibleNextStates (double transitionPrTable[], int sizeOfTable, int currentState);
void setParametersSpellingModel();
/************************************************************************/
//Reset the parameters of the spelling model
/************************************************************************/

void displayParametersSpellingModel();
/************************************************************************/
//Display the parameters of the spelling model
/************************************************************************/

void setParametersKbModel();
/************************************************************************/
//Reset the parameters of the keyboard model
/************************************************************************/


void displayParametersKbModel();
/************************************************************************/
//Display the parameters of the keyboard model
/************************************************************************/



int distanceOfLetters(char charOne, char charTwo);

double getScaleFactor(int num);

void testPrOfTypingDocument();

void testgetPrFromDocument();

void testRecoverOriginal();

void testRecoverOriginal2();

void testPrecisionOfRecovered();

int take1SampleFrom1PrSpace(double prTable[], int sizeOfTable);

void getKeyboardProbabilityTable(char charToType, double prTable[]);

char typeOneChar(char charToType);

void typeOneWord(char word[], char output[], bool traceON = false, int maxOutput = 100);

void typeOneArticle(char * corruptedMessageFile, char * sourceArticle, bool trace = false);

double prOf1CharSeriesWhenTyping1Word(char observedString[], char wordString[]);

double logPrOfGettingDocument1WhenTypingDocument2(string document1, string document2);

void getPrFromDocument(string document1, string document2);

void recoverOriginal(string doc1, string doc2, string doc3);

void recoverOriginal2(string doc1, string doc2, string doc3, string doc4);

double getPrecisionOfRecovered(string doc1, string doc2);
