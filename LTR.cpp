#include <string>
#include <unordered_map>
#include <fstream>
#include <iostream>
#include <vector>
#include <algorithm>

using namespace std;

void checkltr(string filename){
  string line;
  ifstream myfile(filename);
  unordered_map<string, vector<int>> map1;
  string text;
  int count1 =0;
  while (getline(myfile, line)){
    count1++;
    if (count1 > 1)
      text+=line;
  }
  int count= 20;
  int i=0;
  int j=0;
  while (i+20 < text.size()){
    string seqcheck = text.substr(i,20);
    map1[seqcheck].push_back(i);
    i++;
  }

  unordered_map<string, vector<int>>::iterator it;
  it= map1.begin();

  unordered_map<string,vector<int>> map2;
  while (it != map1.end()){
    int left=1;
    int right=1;
    string new_left= "";
    string new_right="";
    if (it->second.size() > 1){
      int j=0;
      for (int i=0; i<it->second.size()-1; i++){
        j= i+1;
        if ((it->second[j]-it->second[i] >= 10000) && (it->second[j]-it->second[i] <= 20000)){
            while ((text[it->second[i]-left] == text[it->second[j]-left])){
              new_left=text[it->second[i]-left] + new_left;
              left++;
            }
    		    while ((text[it->second[i]+19+right] == text[it->second[j]+19+right])){
              new_right+=text[it->second[i]+19+right];
              right++;
            }
            string ltr= new_left+it->first+new_right;
            if (map2.find(ltr) == map2.end()){
              map2[ltr].push_back(it->second[i]-left);
              map2[ltr].push_back(it->second[j]-left);
            }
        }

      }

    }
    it++;
  }

  unordered_map<string, vector<int>>::iterator itt;
  itt= map2.begin();
  while (itt != map2.end()){
      cout << itt->first << ", "<< endl;
      cout << itt->first.length();
      for (int i=0; i<itt->second.size();i++){
        cout << itt->second[i] << ",";
      }
      cout << endl;
    itt++;
  }
}

int main(){
  checkltr("sequence.fasta");
  return 0;
}
