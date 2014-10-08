#include <iostream>
#include <cmath>
#include <iomanip>

using namespace std;

double f_comove (double x);
double f_shech (double x, double psi_st, double alpha, double M_star);
double shechter_integration_function (double a, double b, double c, double m1, double m2, double divisions);
double volume_calculation_function (double a, double d1, double d2);
double redshift_to_comoving_distance_function (double z, double divisions);


int main ()
{
    double psi_st, alpha, M_star, M_1, M_2, d_1, d_2, z_1, z_2, L_1, L_2;
    double shechter_fn_result, volume_result, arcsecs, divisions;

    cout << fixed << setprecision (30);
    cout << endl << "Insert number of arcseconds^2" << endl;
    cin >> arcsecs;


    psi_st=0.0007; // Mpc-3
    //L_star=138600000000000000000000000000000000000000000 erg s^(-1)
    M_star=-20.11;
    alpha=-1.72;
    z_1=5.0;
    z_2=15.0;
    divisions = 10;
    M_1 = -28;
    M_2 = -8;
    d_1 = redshift_to_comoving_distance_function(z_1, divisions);
    d_2 = redshift_to_comoving_distance_function(z_2, divisions);
    cout << endl << "distance 1 = " << d_1 << "m\n" << "distance 2 = " << d_2 << "m\n";
    // M_1 = luminosity_to_magnitude()
    shechter_fn_result = shechter_integration_function(psi_st, alpha, M_star, M_1, M_2, divisions);
    volume_result = volume_calculation_function(arcsecs, d_1, d_2);
    cout << "shechter = " << shechter_fn_result << endl << "volume = " << volume_result << "Mpc^3" << endl << endl;
    cout << endl << " Number Density per Arcsec squared = " << shechter_fn_result*volume_result << endl;
    return 0;
}



double shechter_integration_function(double psi_st, double alpha, double M_star, double m1, double m2, double divisions)
{
    double k;
	double rap, sum, trap;
	sum = 0.0;
	//calculate sum
	for (k=1; k<=(divisions-1); k++)
	{
		rap = f_shech (k*(m2-m1/divisions), psi_st, alpha, M_star);
		sum = sum + rap;
    }
	//integrate using trapezium rule
	trap = (-m2+m1/divisions)* (((f_shech (m1, psi_st, alpha, M_star)+ f_shech (m2, psi_st, alpha, M_star))/2)+ sum);
	//cout << (-m2m1/divisions);
	return trap;
}




double volume_calculation_function (double a, double d1, double d2)
{
    double steradianconversion=pow((3.14159/(180*3600)), 2), meters3, Mpc3;
    meters3 = (a*steradianconversion/3*(pow(d2, 3)-pow(d1, 3)));
    Mpc3= meters3*(pow (3.241*pow(10,(-23)) , 3));
    return Mpc3;
}


double redshift_to_comoving_distance_function(double z, double divisions)
{
    double k;
	double rap, sum, trap;
	sum = 0.0;
	//calculate sum
	for (k=1; k<=(divisions-1); k++)
	{
		rap = f_comove (k*(z/divisions));
		sum = sum + rap;
	}
	//integrate using trapezium rule
	trap = (z/divisions)* (((f_comove (0)+ f_comove (z))/2)+ sum);
	return trap;
}




double f_comove (double x)
{
	double i;
	i=(3*pow(10, 8))/(2.2683*pow(10, -18))*(1/(sqrt(0.3*pow(1+x, 3)+8.24*pow(10,-5)*pow(1+x, 4)+0.7)));
	return (i);
}

double f_shech (double x, double psi_st, double alpha, double M_star)
{
	double i;
	i = (log(10)/2.5)*psi_st*(pow(10.0, (0.4*(alpha+1)*(x-M_star))))*exp(-pow(10.0, (0.4*(x-M_star))));
	//i=(3*pow(10, 8))/(2.2683*pow(10, -18))*(1/(sqrt(0.3*pow(1+x, 3)+8.24*pow(10,-5)*pow(1+x, 4)+0.7)));
	//cout << "this = " << i << endl;
	return (i);
}









    //typical value for
    //psi_st = 8e-3 Mpc^-3
    //L star = 1.4e10 times 3.9e33 erg s^-1 (solar luminosity) = 1.386e44
    //alpha = -0.7
    /* PUT IN for FINAL RELEASE
    cout << "Insert values of psi_st, alpha, M_star, M_1, M_2, z_1, z_2" << endl;
    cin >> psi_st >> alpha >> M_star >> M_1 >> M_2 >> z_1 >> z_2;
    cout << endl << "Insert number of divisions" << endl;
    cin >> divisions;
    */
