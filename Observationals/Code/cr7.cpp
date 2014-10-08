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

double comoving_distance(double z1, double z2, int N=150);
	// calculates comoving distance in Mpc between two z's.

double luminosity_distance(double z, int N=150);
	// calculates luminosity distance in Mpc between z and z=0.

double schechterM(double M, double z);
    // calculates number density /Mpc3 /magnitude at input absolute magnitude and z

double schechterm(double m, double z);
	// calculates number density /Mpc3 /magnitude at input apparent magnitude and z

double logschechterM(double M, double z);
    // calculates phi in log space.

double logschechterm(double m, double z);
	// calculates phi in log space.

double luminosity_density(double m, double z);
    // i dont even know

datapoint total_N(double z1, double z2, int shellN, double m1, double m2, int binN);
	// calculates the number of galaxies in a specified z interval and m interval per unit solid angle.

void z_distribution(double z1, double z2, double Dz, int shellN, double m1, double m2, int binN);
	//

void m_distribution(double z1, double z2, int shellN, double m1, double m2, double Dm, int binN);
	//

void m_z_distribution(double z1, double z2, double Dz, int shellN, double m1, double m2, double Dm, int binN);
	// prints table to output file giving number of galaxies per square degree for each magnitude and z interval

void N_input();
	// runs input system to find total_N()

double cosmic_variance(double V, double X, double N);
	//

/*//////////////////////////////////////////////////////////////////////////////////////////////////
				class definitions
//////////////////////////////////////////////////////////////////////////////////////////////////*/

struct datapoint
{
	double N, V, upper, lower, error;

	datapoint();

	datapoint operator+ (datapoint A);
	datapoint operator+= (datapoint A);
};

datapoint::datapoint() : N(0), V(0), upper(0), lower(0), error(0)
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
    o << d.N;// << " % " << d.error*d.N*0.01;
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


double schechterM(double M, double z)
{
    // calculate schechter paramaters at z
	double phiX  = 0.0027*pow(e,-0.161*z),
	   	   MX	 = 0.2097*z - 21.59,
	   	   alpha = -0.0113*z - 1.6771,

    //calculate schecter result
	   	   ten = pow( 10, (0.4*(MX-M)) ),
	   	   phi = 0.4 * phiX * log(10) * pow( ten, alpha+1 ) * exp( -ten );

	return phi;
}

double schechterm(double m, double z)
{
    double M = to_absolute_magnitude(m,z);

	return schechterM(M, z);
}

double logschechterM(double M, double z)
{
    // calculate paramaters at z
	double phiX  = 0.0027*pow(e,-0.161*z),
	   	   MX	 = 0.2097*z - 21.59,
	   	   alpha = -0.0113*z - 1.6771,

    // calculate schechter result in log space
	   	   index = (0.4*(MX-M)),
	   	   constant = log10( 0.4 * phiX * log(10) ),

	   	   logphi = constant + index*(alpha+1) - ( log10(e)*pow( 10, index ) );

	return logphi;
}

double logschechterm(double m, double z)
{
    double M = to_absolute_magnitude(m,z);

	return logschechterM(M,z);
}

datapoint total_N(double z1, double z2, int shellN, double m1, double m2, int binN)
{
    // ensure the second paramater is the larger for m's and z's
    if (z1>z2) swap(z1,z2);
    if (m1>m2) swap(m1,m2);

    // calculate pillar widths
	double dz = (z2-z1)/shellN,
		   dm = (m2-m1)/binN,

		   N = 0,
		   V = 0;

    // calculate total number of galaxies by summing over every z and m interval
	for ( double z=z1 ; z<z2 ; z+=dz )
	{
	    //calculate shell volume
		double dV = (4/3) * pi * ( pow( comoving_distance(0,z+dz), 3 ) - pow( comoving_distance(0,z), 3 ) );

        V += dV;

		for ( double m=m1 ; m<m2 ; m+=dm )
		{
			double dN = schechterm(m,z) * dm * dV;

			N += dN;
		}
	}

	datapoint d;

	d.V = V;
	d.N = N / 41252.96;
	d.error = cosmic_variance(V,1,1);

    //return # of galaxies per unit solid angle (degree2)
    // 4pi steradians = 41252.96 degrees2 = 148510660 minutes2 = 5.346384e11 seconds2
	return d;
}

void z_distribution(double z1, double z2, double Dz, int shellN, double m1, double m2, int binN)
{
    if (z1>z2) swap(z1,z2);

    // print histogram of galaxy numbers across z
	for ( double z=z1 ; z<z2 ; z+=Dz )
	{
		cout << z << '\t' << total_N(z,z+Dz,shellN,m1,m2,binN) << endl;
	}
}


void m_distribution(double z1, double z2, int shellN, double m1, double m2, double Dm, int binN)
{
    if (m1>m2) swap(m1,m2);

    // print histogram of galaxy numbers across m
	for ( double m=m1 ; m<m2 ; m+=Dm )
	{
		cout << m << '\t' << total_N(z1,z2,shellN,m,m+Dm,binN) << endl;
	}
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

double cosmic_variance(double V, double X, double N)
{
	double log = log10(V),
		   zeta = ( 1-0.03*sqrt(X-1) ) * ( 219.7 - 52.4*log + 3.21*log*log ); // /sqrt(N)

	return zeta;
}

/*//////////////////////////////////////////////////////////////////////////////////////////////////
				main
//////////////////////////////////////////////////////////////////////////////////////////////////*/


int main()
{


	//N_input();

	m_z_distribution_input();

	//m_z_distribution(6,15,1,10,25,30,0.1,10);

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



*/

