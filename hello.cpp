#include <iostream>
#include <string>
#include <fstream>


using namespace std;



int main () {
  //
  srand(time(0));
  ofstream output;
  output.open("hello.dat",ios::app);

  int num;



  // output
  cout << "hello world " << endl;
  num = rand()%10;
  cout << "generate a random number between 0 and 9: "<<num<<endl;


  output << "hello world" << endl;
  output <<num<<endl;

  return 0;
}
