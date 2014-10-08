#include<iostream>
#include<iomanip>
#include<cmath>
using namespace std;
double f (double M, double phi, double Mstar, double alpha, double DL)
{
  double i = 0.4*log(10.0)*phi*pow(10.0, (0.4*(M-Mstar-(5*(log10(DL)-1)))*(alpha+1)))*exp(-pow(10.0, (0.4*(M-Mstar-(5*(log10(DL)-1))))));
	return (i);
}
double g (double z)
{
  double j = ((3*pow(10.0, 5.0))/71)*(1.0/(sqrt(0.27*pow((1.0+z), 3.0)+8.24*pow(10.0, -5.0)*pow((1.0+z), 4.0)+0.73)));
	return (j);
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
  cout << "Input M*" << endl;
  cin >> Mstar;
  cout << "Input alpha" << endl;
  cin >> alpha;
  cout << "Input upper limit magnitude" << endl;
  cin >> M1;
  cout << "Input lower limit magnitude" <<endl;
  cin >> M2;
  double n, k, l, z, o, p ,q;
  double rap2, sum2, trap, tsum, sum1, rap1, DL, dens, rap3, sum3, rap4, sum4, z3, z4, r3, r4, DV;
  cout << "Input number of intervals" << endl;
  cin >> n;
  tsum=0.0;
  for (l=0; l<=n; l++)
  {
    //split into shells of redshift
    z= z1+(((z2-z1)*l)/n);
    sum1 = 0.0;
    //for each redshift shell calculate luminosity distance for conversion to apparent magnitude
    for (o=1; o<=(n-1); o++)
    {
      rap1 = g(o *(z/n));
      sum1 = sum1 + rap1;
    }
    DL = (1+z)*((z/n)*(((g(0) + g(z))/2.0) + sum1))*pow(10.0, 6.0);
    sum2 = 0.0;
    //Integrate apparent magnitude schechter function to get number density
    for (k=1; k<=(n-1); k++)
    {
      rap2 = f((M1+(k*((M2-M1)/n))), phi, Mstar, alpha, DL);
      sum2 = sum2 + rap2;
    }
    dens = ((M2-M1)/n)* (((f (M1, phi, Mstar, alpha, DL)+ f (M2, phi, Mstar, alpha, DL))/2.0)+ sum2);
    sum3 = 0.0;
    //Calculate a comoving volume for each redshift shell
    if (l==n)
      {
	z3 = z;
      }
    else
      {
	z3 = z + ((z2-z1)/(2.0*n));
      }
    for (p=1; p<=(n-1); p++)
      {
	rap3 = g(p*(z3/n));
	sum3 = sum3 + rap3;
      }
    r3=(z3/n)*(((g(0) + g(z3))/2.0) + sum3);
  sum4= 0.0;
  if (l==0)
    {
      z4 = z;
    }
  else
    {     
      z4 = z - ((z2-z1)/(2.0*n));
    }
  for (q=1; q<=(n-1); q++)
    {
      rap4 = g(q*(z4/n));
      sum4 = sum4 + rap4;
    }
  r4=(z4/n)*(((g(0) + g(z4))/2.0) + sum4);
  DV = (1/3.0)*(pow(r3, 3.0)-pow(r4, 3.0));
  //Times the comoving volume by density to get number of galaxies in redshift shell and then sum these shells
trap = DV * dens;
    tsum = tsum + trap;
    cout << "lower limit z"  << "   " << "upper limit radius"<< "   " << "Co-moving Volume" << "   " << "no. of galaxies" << "   " << "summed number" << endl;
    cout << "     " << z4 << "             " << z3 << "               " << DV << "            " << trap << "        " << tsum << endl;
  }
  return 0;
}
