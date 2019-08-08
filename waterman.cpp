#include <iostream>
#include <unordered_map>
#include <string>
#include <fstream>
#include <cmath>
#include <cstdlib>
#include <algorithm>
#include <list>
#include <sstream>
#include <stdio.h>
#include <cstring>
#include "sequence.h"
using namespace std;

double similarity_score(char a , char b, double match, double mu)
{
	if (a == b)
		return match;
	return mu;
}

int find_array_max(int array[], int length, int& index)
{
	int max = array[0];            // start with max = first element
	index = 0;
	for(int i = 1; i < length; i++)   // finds the max value (and index) in an array
	{
		if(array[i] >= max)
		{
			max = array[i];
			index = i;
		}
	}
	return max;                    // return highest value in array
}


int aligner(const string& sequence1, const string& sequence2, int lowerbound)
{
	double delta = -3;            // indel penalty

	string seq_a = sequence1;
	string seq_b = sequence2;

	int N_a = seq_a.length();     // get the actual lengths of the sequences
	int N_b = seq_b.length();

	// cout << N_a << " " << N_b << endl;

	// int H[N_a + 1][N_b + 1];   // Initialize dynamic programming matrix H
	// short I_i[N_a + 1][N_b + 1],
	//       I_j[N_a + 1][N_b + 1];    // Index matrices to remember the 'path' for backtracking

	int** H = new int* [N_a+1];
	short** I_i = new short* [N_a+1];
	short** I_j = new short* [N_a+1];
	for (int i=0; i<= N_a; i++){
		H[i] = new int[N_b+1];
		I_i[i] = new short[N_b+1];
		I_j[i] = new short[N_b+1];
	}


	H[0][0] = 0;
	I_i[0][0] = I_j[0][0] = 0;

	for (int i=1; i <= N_a; i++)
	{
		H[i][0] = 0;
		I_i[i][0] = 0;
	}
	for (int j=1; j <= N_b; j++)
	{
		H[0][j] = 0;
		I_i[0][j] = 0;
	}

	//Scoring Matrix
	//                    A   C   G   T   R   Y   K   M   S   W   B   D   H   V   N
	char nuc44[15][15]= {{1, -2, -2, -2,1,  -2, -2,1,  -2,1,  -2,1, 1, 1, 1},   // A
			    { -2,1,  -2, -2, -2,1,  -2,1, 1,  -2,1,  -2,1, 1, 1},   // C    (Candidate strand)
			    { -2, -2,1,  -2,1,  -2,1,  -2,1,  -2,1, 1,  -2,1, 1},   // G
			    { -2, -2, -2,1,  -2,1, 1,  -2, -2,1, 1, 1, 1,  -2,1},   // T
			    {1,  -2,1,  -2,1,  -2, -2, -2, -2, -2, -2,1,  -2,1, 1},   // R
			    { -2,1,  -2, -2, -2,1,  -2, -2, -2, -2,1,  -2,1,  -2,1},   // Y    (Candidate strand)
			    { -2, -2,1,  -2, -2, -2,1,  -2, -2, -2,1, 1,  -2, -2,1},   // K
			    { -2, -2, -2,1,  -2, -2, -2,1,  -2, -2, -2, -2,1, 1, 1},   // M
			    {1,  -2, -2, -2, -2, -2, -2, -2,1,  -2,1,  -2, -2,1, 1},   // S
			    { -2,1,  -2, -2, -2, -2, -2, -2, -2,1,  -2,1, 1,  -2,1},   // W    (Candidate strand)
			    { -2,1, 1, 1,  -2,1, 1,  -2,1,  -2,1,  -2, -2, -2,1},   // B
			    {1,  -2,1, 1, 1,  -2,1,  -2, -2,1,  -2,1,  -2, -2,1},   // D
			    {1, 1,  -2,1,  -2,1, 1,  -2, -2,1,  -2, -2,1,  -2,1},   // H
			    {1, 1, 1,  -2,1,  -2, -2,1, 1,  -2, -2, -2, -2,1, 1},   // V    (Candidate strand)
			    { 0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0} }; // N

	unordered_map<char,int> nucToInt;
	nucToInt['A'] = 0;
	nucToInt['C'] = 1;
	nucToInt['G'] = 2;
	nucToInt['T'] = 3;
	nucToInt['R'] = 4;
	nucToInt['Y'] = 5;
	nucToInt['K'] = 6;
	nucToInt['M'] = 7;
	nucToInt['S'] = 8;
	nucToInt['W'] = 9;
	nucToInt['B'] = 10;
	nucToInt['D'] = 11;
	nucToInt['H'] = 12;
	nucToInt['V'] = 13;
	nucToInt['N'] = 14;

	int temp[4];                                      // holds the four possible values for each matrix position

	// Smith-Waterman algorithm

	for(int i = 1; i <= N_a; i++)
	{
		for(int j = 1; j <= N_b; j++)
		{

			temp[0] = H[i-1][j-1] + nuc44[ nucToInt[ seq_a[i-1] ] ][ nucToInt[ seq_b[j-1] ] ];
			temp[1] = H[i-1][j] + delta;
			temp[2] = H[i][j-1] + delta;
			temp[3] = 0;

			int index;
			H[i][j] = find_array_max(temp, 4, index);

			switch(index)
			{
				case 0:                                  // score in (i,j) stems from a match/mismatch
					I_i[i][j] = i-1;
					I_j[i][j] = j-1;
					break;
				case 1:                                  // score in (i,j) stems from a deletion in sequence A
					I_i[i][j] = i-1;
					I_j[i][j] = j;
					break;
				case 2:                                  // score in (i,j) stems from a deletion in sequence B
					I_i[i][j] = i;
					I_j[i][j] = j-1;
					break;
				case 3:                                  // (i,j) is the beginning of a subsequence
					I_i[i][j] = i;
					I_j[i][j] = j;
					break;
			}
		}
	}

	int matrix_max = 0;
	int i_max=0, j_max=0;
	for(int i=1;i<=N_a;i++)
	{
		for(int j=1;j<=N_b;j++)
		{
			if(H[i][j]>matrix_max)
			{
				matrix_max = H[i][j];
				i_max=i;
				j_max=j;
			}
		}
	}

	int current_i = i_max,
	    current_j = j_max;
	int next_i = I_i[current_i][current_j];
	int next_j = I_j[current_i][current_j];
	int tick = 0;

	char consensus_a[N_a + N_b + 2],     // two strings hold the alignment
	     consensus_b[N_a + N_b + 2];

	while(((current_i != next_i) || (current_j != next_j)) && (current_i != 0) && (current_j != 0))
	{
		if(next_i == current_i)
			consensus_a[tick] = '-';                    // deletion in A
		else
			consensus_a[tick] = seq_a[current_i - 1];   // match/mismatch in A

		if(next_j == current_j)
			consensus_b[tick] = '-';                    // deletion in B
		else
			consensus_b[tick] = seq_b[current_j - 1];   // match/mismatch in B


		H[current_i][current_j] = 0;  // reset position in H so not found again

		current_i = next_i;
		current_j = next_j;
		next_i = I_i[current_i][current_j];
		next_j = I_j[current_i][current_j];
		tick++;
	}
	consensus_a[tick] = '\0';
	consensus_b[tick] = '\0';
	cout << consensus_a << endl << consensus_b << endl;

	for(int i=0; i<=N_a; i++){
		delete [] H[i];
		delete  [] I_i[i];
		delete [] I_j[i];
	}
	delete [] H;
	delete [] I_i;
	delete [] I_j;

	return matrix_max;

}

const string tab = "      ";

string readfasta(string filename){
	string line;
	int lineNum = 0;
	ifstream myfile(filename);
	string sequence = "";
	if (myfile.is_open())
	{
		while ( getline (myfile,line) )
		{
			lineNum++;
			if (lineNum >1)
			 sequence += line;
		}
		myfile.close();
 }
 return sequence;
}

unordered_map<int,int> minVal(vector<int> vec1){
    vector<int>::iterator it;
    unordered_map<int, int> map1;
    it= std::min_element(vec1.begin(), vec1.end());
    map1[*it] = std::distance(vec1.begin(), it);
    return map1;
}

vector<int> topTen(int num, vector<int>& veccheck){

    unordered_map<int, int> map2= minVal(veccheck);
    unordered_map<int, int>::iterator itt;
    itt = map2.begin();
    while (itt != map2.end()){
      if (num > itt->first)
        veccheck[itt->second]= num;
      itt++;
    }

  return veccheck;
}

string translate(string sequence, int offset){
	const int codonNum = 3;
	string new_str="";
	string subb = "";
	int i=0;
	while (offset <sequence.length()){
		subb = sequence.substr(offset,codonNum);
		if (subb.length()==3)
			new_str+=subb;
		offset+=3;
	}
	return new_str;
}

void processLine(char line[10000], Sequence& obj) {
  char * ptr_ch;
  ptr_ch = strtok (line,"\t");
  int index = 0; //index of the element in the sequence
  while (ptr_ch != NULL) //parse the string line by tab
  {
    obj.fillInObject(index, ptr_ch);
    ptr_ch = strtok (NULL, "\t");
    index++;
  }
  /*cout << obj.id << endl;
  cout << obj.refName << endl;
  cout << obj.pos << endl;
  cout << obj.seq << endl;*/
  return;
}

std::vector<Sequence> returnVector(string filename){
  std::vector<Sequence> allInfo;
  Sequence obj;

  string line1;
  ifstream myfile(filename);
  if (myfile.is_open())
  {
    while ( getline (myfile,line1) )
    {
        char cstr[line1.size() + 1];
        strcpy(cstr, line1.c_str());
        processLine(cstr, obj); //process the line to fill in objects with actual values
        allInfo.push_back(obj); //add each object to the vector}
    }
    myfile.close();
 }

  return allInfo;
}

int main(){

	string readfile = readfasta("sequence.txt");
	vector<Sequence> vec1 = returnVector("viewVirus.txt");
	// string seqcheck1 = "ATGAAAAACCCAAAGAAGAAATTCGGAGGATTCCGGATTGTCAATATGCTAAAACGGGGAGTAGCCC";
	// string seqcheck2 = "GGGG";
	int threshold = 10;
	vector<int> topten;
	string sequencing;

	int myArray[10] = {-1,-2,-3,-4,-5,-6,-7,-8,-9,-10};
  vector<int> topn(myArray, myArray + sizeof(myArray)/sizeof(myArray[0]));

	for (int i=0; i<vec1.size(); i++){
		int scoring= aligner(readfile, vec1[i].seq, threshold);
		topten= topTen(scoring,topn);
		cout << scoring << endl;}
	for (int j=0; j<10; j++)
		cout << topten[j] << endl;


	return 0;
}
