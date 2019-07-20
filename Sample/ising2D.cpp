// to realize the Monte Carlo Simulation of the 2D Ising model
// by using the Metropolis algorithm

#include <cmath> //need for abs equal to fabs
#include <ctime>
#include <iostream>//needed for io
#include <string>
#include <cstdlib>
#include <iomanip>
#include <fstream>
#include <map>
#include <sstream> // needed for internal io
#include <stdexcept>
#include <vector>
#include <cstdio>
#include <iomanip>
#include <stdlib.h>
#include <algorithm>
#include <random>
#include <vector>

using namespace std;

const int IsingRow = 50;  // the size of 2D square lattice
const int IsingCol = 50;
const long steps = 1000000;  // the steps of the MC


void IsingInitial(int (&latt)[IsingRow][IsingCol])
{
  // initialization of a 2D Ising lattice
  int i,j;

  for (i=0;i<IsingRow;i++)
    for (j=0;j<IsingCol;j++)
      {
        //
        latt[i][j]=rand()%2;
        latt[i][j]=latt[i][j]*2-1;
      }
}


void IsingPrint(int (&latt)[IsingRow][IsingCol])
{
  // print out the 2D Ising lattice
  int i,j;

  for (i=0;i<IsingRow;++i)
    {
      for (j=0;j<IsingCol;++j)
        if (latt[i][j]==1)
          cout<<"+ ";
        else if (latt[i][j]==-1)
          cout<<"- ";
      cout<<endl;
    }
}

int IsingSummation(int (&latt)[IsingRow][IsingCol])
{
  // obtain the summation of the 2D Ising lattice
  int i1,i2,sum=0;
  int j1,j2,k1,k2;

  for (size_t i1=0;i1<IsingRow;++i1)
    for (size_t i2=0;i2<IsingCol;++i2)
    {
      // define the indices for the nearest neighbor sites
      j1=i1-1;
      j2=i2-1;
      k1=i1+1;
      k2=i2+1;

      // impose the periodic boundary conditions
      if (j1<0)
        j1=IsingRow-1;
      if (j2<0)
        j2=IsingCol-1;
      if (k1>=IsingRow)
        k1=0;
      if (k2>=IsingCol)
        k2=0;

      sum+=latt[i1][i2]*(latt[j1][i2]+latt[i1][j2]+latt[i1][k2]+latt[k1][i2]);
    }
  sum=sum/2;
  return sum;
}


void IsingFlip(int (&latt)[IsingRow][IsingCol], double tt)
{
  // to select a spin in the 2D Ising lattice randomly
  // determine to flip the spin or not, by checking the energy change
  int i1,i2;
  int j1,j2,k1,k2;
  double dE=0;

  i1=rand()%IsingRow;
  i2=rand()%IsingCol;

  j1=i1-1;
  j2=i2-1;
  k1=i1+1;
  k2=i2+1;
  if (j1<0)
    j1=IsingRow-1;
  if (j2<0)
    j2=IsingCol-1;
  if (k1>=IsingRow)
    k1=0;
  if (k2>=IsingCol)
    k2=0;

  if (latt[i1][i2]==1)
    dE=2*(latt[j1][i2]+latt[i1][j2]+latt[i1][k2]+latt[k1][i2]);
  else
    dE=-2*(latt[j1][i2]+latt[i1][j2]+latt[i1][k2]+latt[k1][i2]);

  // determine to flip the spin or not
  if (dE<0)
    latt[i1][i2]=-latt[i1][i2];
  else
    {
      //
      bool flip=(rand()%10000)<10000*exp(-dE/tt);
      if (flip==true)
        latt[i1][i2]=-latt[i1][i2];
    }

}


double IsingMag(int (&latt)[IsingRow][IsingCol])
{
  // obtain the summation of the 2D Ising lattice
  int i,j;
  double ave=0;

  for (size_t i=0;i<IsingRow;++i)
    for (size_t j=0;j<IsingCol;++j)
      ave+=latt[i][j];

  ave=1.0*ave/IsingRow/IsingCol;
  return ave;
}


int main()
{

  clock_t counting;
  counting=clock();

  ofstream output;
  output.open("ising2D_L50L50.dat",ios::app);
  srand(time(0));


  int Isinglattice[IsingRow][IsingCol];
  int IsinglatticeInit[IsingRow][IsingCol];
  int IsinglatticeFinal[IsingRow][IsingCol];
  int ii,jj;
  long trys,Tlen;
  double IsingSum=0;
  double DeltaE;
  double J=1.0; // the couplings of the nearest site interactions
  double Temp; // the temperature of the system

  // initialization
  IsingInitial(IsinglatticeInit);

  // loops over different temperatures
  for (Tlen=0;Tlen<1000;Tlen++)
    {
      //srand(time(0));
      Temp=0.01+Tlen*0.01;

      cout<<"========================================================="<<endl;
      // return the initial value of the nearest neighbor summation of the 2D Ising lattice
      cout<<"the system temperature is: "<<Temp<<endl;
      cout<<"print out the initial random 2D Ising lattice as follows:"<<endl;
      IsingPrint(IsinglatticeInit);
      IsingSum = -J*IsingSummation(IsinglatticeInit) ;
      cout<<"the total summation over the nearest neighbor sites of the 2D Ising lattice is: "<<IsingSum<<endl;


      for (ii=0;ii<IsingRow;++ii)
        for (jj=0;jj<IsingCol;++jj)
          IsinglatticeFinal[ii][jj]=IsinglatticeInit[ii][jj];

      // iterations of the random initial 2D Ising lattice
      for (trys=0;trys<steps;++trys)
        {
          //
          IsingFlip(IsinglatticeFinal,Temp);
          IsingSum = -J*IsingSummation(IsinglatticeFinal);
        }
      cout<<"print out the 2D Ising lattice after the iteration:"<<endl;
      IsingPrint(IsinglatticeFinal);
      cout<<"After the iteration, the total summation over the nearest neighbor sites of the 2D Ising lattice is: "<<IsingSum<<endl;
      output<<Temp<<" "<<IsingSum<<" "<<IsingMag(IsinglatticeFinal)<<endl;
    }

  //IsingPrint(Isinglattice);

  output.close();
  counting=clock()-counting;
  cout<<"This run took "<<((float)counting)/CLOCKS_PER_SEC<<" seconds."<<endl;

  return 1;
}
