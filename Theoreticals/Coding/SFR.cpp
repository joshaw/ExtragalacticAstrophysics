#include<iostream>
#include<cmath>
using namespace std;
double f (double z)
{
  double i = pow(10.0, 53.5)*0.2*(0.58*exp(-0.135*z)+0.013);
	return i;
}
double g (double z)
{
  double j, critdens, Hcritdens;
  critdens = (3*(70000/(3.086e22))*(70000/(3.086e22)))/(8*3.141592*6.67398e-11)*(0.27*pow((1.0+z), 3.0)+8.24e-5*pow((1.0+z), 4.0)+0.73);
  Hcritdens = (1/(1+6.0e-5+3e-5+2.4e-1+3.0e-9))*critdens;
  j =  Hcritdens/1.66053886e-27;
  return (j);
}
int main ()
{
  double acc, z2;
   cout << "Input upper limit of redshift" << endl;
   cin >> z2;
  // cout << "Input lower limit of redshift" << endl;
  // cin >> z1;
  cout << "Input accuracy" << endl;
  cin >> acc;
  double l, k, z, n, t, tdash, zdash;
  double nion, sum2, rap2, nH, frac;
  cout << "Input number of intervals" << endl;
  cin >> n;  
  l=0.0;
      while (frac < acc)
      // for (l=0; l<=n; l++)
    {
      //split into shells of redshift
      //z= z2-(((z2-z1)*l)/n);
      z = z2-(l/10.0);
      l++;
	t = (((13.77e9*3.1557e7)/(pow((1+z), 1.5))));
	sum2 = 0.0;
      for (k=1; k<=(n-1); k++)
        {
	  tdash = k*(t/n);
	  zdash = pow(((13.77e9*3.1557e7)/tdash), (2/3))-1.0;
          rap2 = f(zdash);
          sum2 = sum2 + rap2;
        }
      nion = (t/n)* ((f (t)/2.0)+ sum2);
         nH = g(z)*pow(3.086e22, 3.0);
         frac = nion/nH;
	cout << z << " " << nion << " " << nH << " " << frac << endl;
    }
  return 0;
}
