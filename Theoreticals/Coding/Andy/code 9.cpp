#include <iostream>
#include <iomanip>
#include <vector>
#include <cmath>
#include <fstream>

using namespace std;


/*//////////////////////////////////////////////////////////////////////////////////////////////////
				variables
//////////////////////////////////////////////////////////////////////////////////////////////////*/

const double  pi = 3.1415926,
			  c  = 3e8,
			  e  = 2.7182818,

			  omegaR = 8.24e-5,
			  omegaM = 0.27,
			  omegaL = 0.73,

			  H0 = 71;

struct datapoint;

/*//////////////////////////////////////////////////////////////////////////////////////////////////
				function declarations
//////////////////////////////////////////////////////////////////////////////////////////////////*/

double to_absolute_magnitude(double m, double z);
	// converts apparent to absolute magnitude in ab system

double to_apparent_magnitude(double M, double z);
	// converts absolute to apparent magnitude in ab system

double hubble(double z);
	// calculates hubble parameter at given z in km/s/Mpc

double alpha(double z, unsigned short i);
	// calculates schechter parameter alpha at redshift z
	//		( i=1 for underestimate, i=2 for best estimate, i=3 for overestimate )

double MX(double z, unsigned short i);
	// calculates schechter parameter M* at redshift z
	//		( i=1 for underestimate, i=2 for best estimate, i=3 for overestimate )

double phiX(double z, unsigned short i);
	// calculates schechter parameter phi* at redshift z
	//		( i=1 for underestimate, i=2 for best estimate, i=3 for overestimate )

double comoving_distance(double z1, double z2, int N=80);
	// calculates comoving distance in Mpc between two z's.

double luminosity_distance(double z, int N=80);
	// calculates luminosity distance in Mpc between z and z=0.

double schechterM(double M, double z, unsigned short i, unsigned short j, unsigned short k);
    // calculates number density /Mpc3 /magnitude at input absolute magnitude and z

double schechterm(double m, double z, unsigned short i, unsigned short j, unsigned short k);
	// calculates number density /Mpc3 /magnitude at input apparent magnitude and z

double schechterMparam(double M, double z, double phiX, double alpha, double MX);
	// calculates schechter function value with provided parameters

double schechtermparam(double m, double z, double phiX, double alpha, double MX);
    // calculates schechter function value with provided parameters

datapoint total_N(double z1, double z2, int shellN, double m1, double m2, int binN);
	// calculates the number of galaxies in a specified z interval and m interval per unit solid angle.

void m_z_distribution(double z1, double z2, double Dz, int shellN, double m1, double m2, double Dm, int binN);
	// prints table to output file giving number of galaxies per square degree for each magnitude and z interval

void N_input();
	// runs input system to find total_N()

void m_z_distribution_input();
	// runs input system to find m_z_distribution()

double cosmic_variance(double z1, double z2, double theta, double phi, double N);
    // cosmic variance as a percentage, angles are input in arcseconds, N is the number of pointings

void variance_input();
    // cosmic variance input system

/*//////////////////////////////////////////////////////////////////////////////////////////////////
				class definitions
//////////////////////////////////////////////////////////////////////////////////////////////////*/

struct datapoint
{
	double N, V, upper, lower;

	datapoint();

	datapoint operator+ (datapoint A);
	datapoint operator+= (datapoint A);
};

datapoint::datapoint() : N(0), V(0), upper(0), lower(0)
{

}

datapoint datapoint::operator+ (datapoint A)
{
	datapoint B;

	B.N = A.N + N;
	B.V = A.V + V;
	B.upper = A.upper + upper;
	B.lower = A.lower + lower;

	return B;
}

datapoint datapoint::operator+= (datapoint A)
{
	N += A.N;
	V += A.V;
	upper += A.upper;
	lower += A.lower;

	return *this;
}



// set the default output-behaviour of a datapoint object to simply output parameter N
std::ostream& operator<< (std::ostream &o, const datapoint d)
{
    o << d.N;
    return o;
}

/*//////////////////////////////////////////////////////////////////////////////////////////////////
				function implementation
//////////////////////////////////////////////////////////////////////////////////////////////////*/


double to_absolute_magnitude(double m, double z)
    { return m - 5*( log10(luminosity_distance(z)*1000000) - 1 ); }

double to_apparent_magnitude(double M, double z)
    { return 5*( log10(luminosity_distance(z)*1000000) - 1 ) + M; }

double hubble(double z)
{
    return H0*sqrt( omegaM*pow((1+z),3) + omegaR*pow((1+z),4) + omegaL);
}

double alpha(double z, unsigned short i)
{
	if      (i==1)	{   if (( -0.0261637*z - 1.60201 )>-2) return ( -0.0261637*z - 1.60201 );
						else return -2;
					}
	else if (i==2)	return ( -0.0146437*z - 1.66423 );
	else if (i==3)	return ( -0.0031237*z - 1.72634 );
	else            return 0;
}


double MX(double z, unsigned short i)
{
	if      (i==1)	return ( 0.18700*z - 21.45773 );
	else if (i==2)	return ( 0.22116*z - 21.64204 );
	else if (i==3)	return ( 0.25532*z - 21.82641 );
	else            return 0;
}


double phiX(double z, unsigned short i)
{
	double value = 0.0521538*pow(e,-1.48655*z) + 0.00105869,

		   error = value * sqrt( 1.84*1.84 + 0.5292*0.5292 + 0.1268*0.1268 );

	//if      (i==1)	  if (value-error>0) return ( value - error );
	//	  	  	  	  	  else return 0;
	//	  	  	  	  }
	if      (i==1)  return ( value - error );
	else if (i==2)	return ( value );
	else if (i==3)	return ( value + error );
	else            return 0;
}

double comoving_distance(double z1, double z2, int N)
{
    // use trapezium rule to solve comoving distance integral

	double T=0, h=(z2-z1)/N;

	for (int i=1; i<N-1; i++)
	{
		T += c/hubble( z1 + h*i );
	}

	return (0.5*h*( c/hubble(z1) + c/hubble(z2) ) + h*T) / 1000;
}


double luminosity_distance(double z, int N)
{
	return (1+z)*comoving_distance(0,z,N);
}


double schechterM(double M, double z, unsigned short i, unsigned short j, unsigned short k)
{
    // calculate schechter paramaters at z
	double phiXi  = phiX(z,i),
	   	   MXi	  = MX(z,j),
	   	   alphai = alpha(z,k),

    //calculate schecter result
	   	   ten = pow( 10, (0.4*(MXi-M)) ),
	   	   phi = 0.4 * phiXi * log(10) * pow( ten, alphai+1 ) * exp( -ten );

	return phi;
}

double schechterm(double m, double z, unsigned short i, unsigned short j, unsigned short k)
{
    double M = to_absolute_magnitude(m,z);

	return schechterM(M,z,i,j,k);
}


double schechterMparam(double M, double phiX, double alpha, double MX)
{
	double ten = pow( 10, (0.4*(MX-M)) ),
   		   phi = 0.4 * phiX * log(10) * pow( ten, alpha+1 ) * exp( -ten );

	return phi;
}


double schechtermparam(double m, double z, double phiX, double alpha, double MX)
{
    double M = to_absolute_magnitude(m,z);

	return schechterMparam(M,phiX,alpha,MX);
}

datapoint total_N(double z1, double z2, int shellN, double m1, double m2, int binN)
{
    // ensure the second paramater is the larger for m's and z's
    if (z1>z2) swap(z1,z2);
    if (m1>m2) swap(m1,m2);

    // calculate pillar widths
	double dz = (z2-z1)/shellN,
		   dm = (m2-m1)/binN,

		   N1 = 0,
		   N2 = 0,
		   N3 = 0,
		   V = 0;

    // calculate total number of galaxies by summing over every z and m interval
	for ( double z=z1 ; z<z2 ; z+=dz )
	{
	    //calculate shell volume
		double dV = (4/3) * pi * ( pow( comoving_distance(0,z+dz), 3 ) - pow( comoving_distance(0,z), 3 ) );

        V += dV;

		for ( double m=m1 ; m<m2 ; m+=dm )
		{
			double dN1 = schechterm(m,z,1,3,3) * dm * dV;
			double dN2 = schechterm(m,z,2,2,2) * dm * dV;
			double dN3 = schechterm(m,z,3,1,1) * dm * dV;

			N1 += dN1;
			N2 += dN2;
			N3 += dN3;
		}
	}

	datapoint d;

	d.V = V;
	d.lower = N1 / 41252.96;
	d.N     = N2 / 41252.96;
	d.upper = N3 / 41252.96;

    //return # of galaxies per unit solid angle (degree2)
    // 4pi steradians = 41252.96 degrees2 = 148510660 minutes2 = 5.346384e11 seconds2
	return d;
}


void m_z_distribution(double z1, double z2, double Dz, int shellN, double m1, double m2, double Dm, int binN)
{
    if (z1>z2) swap(z1,z2);
    if (m1>m2) swap(m1,m2);

    vector< vector<datapoint> > table;		// make a table to store all data

	for ( double m=m1 ; m<m2 ; m+=Dm )		// use total_N() function on a range of z and m intervals and output all results as a table in an output file
	{
		vector<datapoint> row;		// vector to store row data

        cout << endl << "  " << m << '\t';

		for ( double z=z1 ; z<z2 ; z+=Dz )
		{
		    row.push_back( total_N(z,z+Dz,shellN,m,m+Dm,binN) );

            cout << '.';
		}

		table.push_back(row);
	}

    /////////////////////////////////// print

	ofstream fout("output.txt");

	fout << "     ";

	for ( double z=z1 ; z<z2 ; z+=Dz )
	{
		fout << setw(11) << z;
	}

	fout << endl << endl;
	fout << setprecision(3);
	fout.setf(ios::fixed);

	for ( unsigned int i=0; i<table.size(); i++ )
	{
		fout << setw(4) << m1+(i*Dm) << " |";

		for ( unsigned int j=0; j<table[0].size(); j++ )
		{
			fout << setw(11) << table[i][j];
		}

        datapoint row;
        for ( unsigned int a=0; a<table[i].size(); a++) // calculate row total
            row += table[i][a];

        fout << "  |" << setw(11) << row << endl;
	}

	fout << endl <<  "totals |";

	for ( unsigned int a=0; a<table[0].size(); a++ ) // calculate column totals
	{
		datapoint column;

		for ( unsigned int b=0; b<table.size(); b++ )
            column += table[b][a];

		fout << setw(11) << column;
	}

	datapoint total;

	for ( unsigned int i=0; i<table.size(); i++ )   // calculate total total
        for ( unsigned int j=0; j<table[0].size(); j++ )
            total += table[i][j];

    fout << "  |" << setw(11) << total << endl << endl;

	fout.close();

    ////////////////////////////////////
}


void N_input()
{

	double z1, z2, m1, m2, shellN, binN;

    while (1)
    {
        cout << endl;
        cout << " z1: ";    cin  >> z1;
        cout << " z2: ";    cin  >> z2;
        cout << " m1: ";    cin  >> m1;
        cout << " m2: ";    cin  >> m2;
        cout << " shellN: ";cin  >> shellN;
        cout << " binN: ";  cin  >> binN;

        cout << endl << " searching for galaxies..." << endl;

        cout << endl
             << " N = " << total_N(z1,z2,400,m1,m2,400) // 7 to 11 with 16700
             << "  / degree2" << endl;                  // 5 to 25 with 6000000

        cout << endl << endl << endl;
    }
}


void m_z_distribution_input()
{
	double z1, z2, Dz, m1, m2, Dm, shellN, binN;

    while (1)
    {
	    cout << endl;
	    cout << " z1: ";    cin  >> z1;
	    cout << " z2: ";    cin  >> z2;
	    cout << " Dz: ";    cin  >> Dz;
	    cout << " m1: ";    cin  >> m1;
	    cout << " m2: ";    cin  >> m2;
	    cout << " Dm: ";    cin  >> Dm;
	    cout << " shellN: ";cin  >> shellN;
	    cout << " binN: ";  cin  >> binN;

	    cout << endl << " searching for galaxies..." << endl << endl;

		m_z_distribution(z1,z2,Dz,shellN,m1,m2,Dm,binN);

	    cout << endl << endl << " calculation complete.";
		cout << endl << endl << endl;
	}

}

double cosmic_variance(double z1, double z2, double theta, double phi, double N)
{
	if (phi>theta) swap(theta,phi);
	if (z1>z2)     swap(z1,z2);

	double X = theta/phi,
		   C = comoving_distance(z1,z2,1000),
		   r = comoving_distance(z1,z1+0.5*(z2-z1),1000),
		   area = theta*phi * pow(pi/180,4),
		   log = log10( area * 291 ),
		   zeta = ( 1-0.03*sqrt(X-1) ) * ( 219.7 - 52.4*log + 3.21*log*log ) / sqrt(N*C/291); // /sqrt(N)

	return zeta;
}

void variance_input()
{
	double z1, z2, theta, phi, N;

    while (1)
    {
	    cout << endl;
	    cout << " z1: ";     cin  >> z1;
	    cout << " z2: ";     cin  >> z2;
	    cout << " angle1: "; cin  >> theta;
	    cout << " angle2: ";   cin  >> phi;
	    cout << " N: ";     cin  >> N;

	    cout << endl;
		cout << " cosmic variance = " << cosmic_variance(z1,z2,theta,phi,N) << " %";
		cout << endl << endl << endl;
	}
}

/*//////////////////////////////////////////////////////////////////////////////////////////////////
				main
//////////////////////////////////////////////////////////////////////////////////////////////////*/


int main()
{
	//variance_input();

	//m_z_distribution(6,15,1,50,25,30,0.2,50);

	m_z_distribution_input();

	//ofstream fout("shellconverge.txt");

	//for (int N=1; N<500; N++)
	//{
	//	  cout << N << '\t' << total_N(6,7,100,25,30,N) << endl;
	//	  fout << N << '\t' << total_N(6,7,100,25,30,N) << endl;
	//}
//
	//fout.close();

//	  double phiXi = phiX(7,2),
//	     	 alphai= alpha(7,2);
//
//	  for (double M=-25; M<-15; M+=0.1)
//	  {
//	  	  fout << endl << M << ' ';
//
//	  	  for (double MXi=-21; MXi<-18; MXi+=0.2)
//
//	  	  	  fout << log10(schechterMparam(M,phiXi,alphai,MXi)) << ' ';
//	  }




//	  for (double M=-30; M<-10; M+=0.2)
//	  {
//	  	  fout << endl << M << ' ';
//
//	  	  for (double z = 5; z<30; z+=1)
//	  	  {
//	  	  	  double phiXi = phiX(z,2),
//	  	  	    	 alphai= alpha(z,2),
//	  	  	     	 MXi   = MX(z,2);
//
//	  	  	   fout << log10(schechterMparam(M,phiXi,alphai,MXi)) << ' ';
//	  	  }
//	  }

//	  vector<double> data(30,0);
//
//	  for (double m=20; m<40; m+=0.2)
//	  {
//	  	  fout << endl << m << ' ';
//
//	  	  for (double z=5; z<30; z+=1)
//	  	  {
//	  	  	  double phiXi = phiX(z,2),
//	  	  	    	 alphai= alpha(z,2),
//	  	  	     	 MXi   = MX(z,2);
//
//	  	  	  data[z] += total_N(z,z+1,20,m,m+0.2,20).N;
//
//	  	  	  fout << log10(data[z]) << ' ';
//	  	  }
//	  }
//
//	  fout.close();



	//N_input();

	//m_z_distribution_input();

	//m_z_distribution(6,15,1,100,25,30,0.1,100);

    //for (double m=-30; m<-10; m+=0.5)
    //{
    //    total_N(10,11,300,m,m+0.5,300);
    //}

	//ofstream fout("converge.txt");

	//for (int N=1; N<1000; N++)
	//	  fout << N << '\t' << total_N(5,10,N,25,30,500) << endl;

	//fout.close();

	return 0;
}

/*//////////////////////////////////////////////////////////////////////////////////////////////////
				plans
//////////////////////////////////////////////////////////////////////////////////////////////////*/
/*

	fix cosmic variance
	implement upper and lower bounds
    converging analysis to optimise accuracy for given computation time

    3D plots for log and non-log

    report:

    overloading operators

*/

