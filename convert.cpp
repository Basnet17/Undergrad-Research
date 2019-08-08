#include <string>
#include <fstream>
#include <iostream>

using namespace std;

void convert(string filename){
  string line;
  ifstream myfile(filename);
  int count=0;
  while (getline(myfile, line)){
    count++;
    if (count <=2){
      if (line[0] == '@'){
        line[0] = '>';
      }
      cout << line << endl;
    }
    if (count ==4)
      count=0;
  }

}

int main(){
  convert("example1.fastq");
  return 0;
}
