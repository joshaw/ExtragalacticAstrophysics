#include <iostream>
#include <cmath>
#include <iomanip>
#include <fstream>
#include <string>

using namespace std;

double f_comove (double x);
double f_shech (double x, double psi_st, double alpha, double M_star);
double shechter_integration (double a, double b, double c, double m1, double m2, double divisions);
double volume_calculation (double a, double d1, double d2);
double magnitude_conversion(double M, double D);
double redshift_to_comoving_distance (double z, double divisions);
void savedata (double results[16]);


int main ()
{
    ///DECLARATION OF VARIABLES
    double psi_st, alpha, M_star, M_1, M_2, d_1, d_2, z_1, z_2, m_1, m_2, m_star;
    double shechter_fn_result, volume_result, arcsecs, divisions, results[16];
    int state=0;

while (state!=2)
{
    ///USER INPUTS
    //Asks for
    cout << "Insert least bright APPARENT magnitude" << endl;
    cin >> m_1;
    cout << "Insert most bright APPARENT magnitude" << endl;
    cin >> m_2;
    cout << "Insert upper z limit "<< endl;
    cin >> z_2;
    cout << "Insert lower z limit " << endl;
    cin >> z_1;

    ///CONSTANTS
    arcsecs = 6291456;
    psi_st=0.0007;
    M_star=-20.11;
    alpha = -1.72;
    divisions = 1000;

    ///CONVERSIONS
    d_1 = redshift_to_comoving_distance(z_1, divisions);
    d_2 = redshift_to_comoving_distance(z_2, divisions);
    M_1 = magnitude_conversion(m_1, (3.24*pow(10, (-17)))*(d_1+d_2)/2);
    M_2 = magnitude_conversion(m_2, (3.24*pow(10, (-17)))*(d_1+d_2)/2);
//    M_1 = -49.0;
//    M_2 = -50.0;

//    ofstream god;
//    god.open("damn2.txt");
//    ///MAIN CALCULATIONS
//    for (int g = 0; g<1000; g++)
//    {
//    M_1+=0.1;
//    M_2+=0.1;
    shechter_fn_result = shechter_integration(psi_st, alpha, M_star, M_1, M_2, divisions);
    volume_result = volume_calculation(arcsecs, d_1, d_2);
//    god << log10(M_2) << "\t" << log10(shechter_fn_result*volume_result) << endl;
//    }
//    god.close();
//    state=2;

    ///OUTPUTTING RESULTS
    cout << endl << "\tRedshift 1  \t\t\t" << z_1 << endl;
    cout << "\tRedshift 2  \t\t\t" << z_2;
    cout << endl << "\tComoving Distance 1  \t\t" << d_1 << "m\n" << "\tComoving Distance 2  \t\t" << d_2 << "m\n";
    cout << "\tShechter  \t\t\t" << shechter_fn_result << endl << "\tVolume  \t\t\t" << volume_result << "Mpc^3" << endl;
    cout << "\tAbsolute Magnitude Range  \t" << M_1 << " to " << M_2 << endl;
    cout << "\tNumber = "<< shechter_fn_result*volume_result << " (per " << arcsecs << " Arcsec squared)" << endl;


    ///SAVING DATA
    cout << endl << "0. Restart, 1. Save, 2. Quit" << endl;
    cin >> state;
    if (state==1)
    {results[1] = psi_st;
    results[2] = M_star;
    results[3] = alpha;
    results[4] = arcsecs;
    results[5] = divisions;
    results[6] = M_1;
    results[7] = M_2;
    results[8] = m_1;
    results[9] = m_2;
    results[10] = z_1;
    results[11] = z_2;
    results[12] = d_1;
    results[13] = d_2;
    results[14] = shechter_fn_result;
    results[15] = volume_result;
    results[16] = shechter_fn_result*volume_result;
    savedata(results);}
}
    return 0;
}


/// DOES NOT WORK
double shechter_integration(double psi_st, double alpha, double M_star, double m1, double m2, double divisions)
{
    double k;
	double rap, sum, trap;
	sum = 0.0;
	for (k=1; k<=(divisions-1); k++)
	{
	    rap = f_shech (m1+k*((m2-m1)/divisions), psi_st, alpha, M_star);
		sum = sum + rap;
    }
    trap = ((m1-m2)/divisions)* (((f_shech (m1, psi_st, alpha, M_star)+ f_shech (m2, psi_st, alpha, M_star))/2)+ sum);
	//trap = ((m2-m1)/divisions)* (((f_shech (m1, psi_st, alpha, M_star)+ f_shech (m2, psi_st, alpha, M_star))/2)+ sum);
	return trap;
}

/// THIS WORKS - UNTESTED
double volume_calculation (double a, double d1, double d2)
{
    double steradianconversion=pow((3.14159/(180*3600)), 2), meters3, Mpc3;
    meters3 = (a*steradianconversion/3*(pow(d2, 3)-pow(d1, 3)));
    Mpc3= meters3*(pow (3.241*pow(10,(-23)) , 3));
    return Mpc3;
}

/// THIS WORKS
double redshift_to_comoving_distance(double z, double divisions)
{
    double k;
	double rap, sum, trap;
    sum=0.0;

	for (k=1; k<=(divisions-1); k++)
	{
		rap = f_comove (k*(z/divisions));
		sum = sum + rap;
	}
	trap = (z/divisions)* (((f_comove (0)+ f_comove (z))/2)+ sum);
	return trap;
}

/// THIS WORKS
double f_comove (double x)
{
	double i;
	i=(3*pow(10, 8))/(2.2683*pow(10, -18))*(1/(sqrt(0.3*pow(1+x, 3)+8.24*pow(10,-5)*pow(1+x, 4)+0.7)));
	return (i);
}

/// DOES NOT WORK
double f_shech (double x, double psi_st, double alpha, double M_star)
{
	double i;
	//i = (log10(10)/2.5)*psi_st*(pow(10.0, (0.4*(alpha+1)*(M_star-x))))*exp(-pow(10.0, (0.4*(M_star-x))));
	i = (log(10)*0.4)*psi_st*exp(((-0.4*log(10))*(alpha+1)*(x-M_star))-exp(-0.4*log(10)*(x-M_star)));
	//cout << ((-0.4*log(10))*(alpha+1)*(x-M_star))-exp(-0.4*(x-M_star))<< endl;
	return (i);
}

/// THIS WORKS
double magnitude_conversion(double m, double D)
{
    double M=m+5*(log10(10)-log10(D));
    return M;
}

/// THIS WORKS
void savedata(double results[16])
{
    string savename;
    ofstream file;
    cout << "enter save name: " << endl;
    cin >>savename;
    savename += ".txt";
    file.open ( savename.c_str()  );
    file << "Variable\tValue\t\tUnit" << endl << endl;
    file << "Psi* \t\t" << results[1] << "\t\t" << "Megaparsecs^(-3)" << endl;
    file << "M* \t\t" << results[2] << "\t\t" << "Solar Masses" << endl;
    file << "Alpha \t\t" <<  results[3] << endl;
    file << "Area \t\t" << results[4] << "\t\tArcsec Squared" << endl;
    file << "Iterations\t" << results[5] << endl;
    file << "Absolute M\t(" << results[6] << ")-(" << results[7] << ")" << endl;
    file << "Apparent M\t(" << results[8] << ")-(" << results[9] << ")" << endl;
    file << "z range\t\t" << results[10] << "-" << results[11] << endl;
    file << "Comoving D1\t" << results[12] << "\tMeters" << endl;
    file << "Comoving D2\t" << results[13] << "\tMeters" << endl << endl;
    file << "Volume\t\t" << results[14] << "\tMegaparsec Cubed" << endl;
    file << "No. Density\t" << results [15] << "\t\tMegaparsecs^(-3)" << endl << endl;
    file << "Number Of Galaxies = " << results [16] << " per " << results[4] << " Arcsec Squared";
    file.close();
    cout << endl << "saved!" << endl << endl;
}




///PROBLEMS
//Its still giving me higher values for working with stupid-bright magnitudes, should these go toward 0 in the shecter fn?
//No idea how it evolves with time, but that needs to be plotted still
//Need to clean up the code, so that everyone wont laugh at me and call me a silly.






            //M_1 = -18; //26 Apparent
            //M_2 = -21; //23 Apparent
/*          PUT IN for FINAL RELEASE
            cout << "Insert values of psi_st, alpha, M_star, M_1, M_2, z_1, z_2" << endl;
            cin >> psi_st >> alpha >> M_star >> M_1 >> M_2 >> z_1 >> z_2;
            cout << endl << "Insert number of divisions" << endl;
            cin >> divisions;
*/
