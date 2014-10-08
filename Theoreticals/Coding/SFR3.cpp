#include<iostream>
#include<iomanip>
#include<cmath>
using namespace std;
double f (double M, double phi, double Mstar, double alpha)
{
  double i = pow(10.0, (-0.4*(M-41.598)))*0.4*log(10.0)*phi*pow(10.0, (0.4*(Mstar-M)*(alpha+1)))*exp(-pow(10.0, (0.4*(Mstar-M))));
  return (i);
}
int main ()
{
  double phi, Mstar, alpha, M1, M2, z1, z2;
  cout << "Input upper limit of redshift" << endl;
  cin >> z2;
  cout << "Input lower limit of redshift" << endl;
  cin >> z1;
  cout << "Input normalisation" << endl;
  cin >> phi;
  cout << "Input alpha" << endl;
  cin >> alpha;
  cout << "Input upper limit magnitude" << endl;
  cin >> M1;
  cout << "Input lower limit magnitude" <<endl;
  cin >> M2;
  double n, k;
  double rap2, sum2, dens;
  cout << "Input number of intervals" << endl;
  cin >> n;
      Mstar = -20.11;
      sum2 = 0.0;
      //Integrate apparent magnitude schechter function to get number density
      for (k=1; k<=(n-1); k++)
	{
	  rap2 = f((M1+(k*((M2-M1)/n))), phi, Mstar, alpha);
	  sum2 = sum2 + rap2;
	}
      dens = ((M2-M1)/n)* (((f (M1, phi, Mstar, alpha)+ f (M2, phi, Mstar, alpha))/2.0)+ sum2);
      cout << dens;
  return 0;
}
