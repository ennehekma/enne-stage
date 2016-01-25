#include <iostream>
#include <fstream>
#include <vector>
#include"/usr/local/astro/include/Astro/orbitalElementConversions.hpp"
#include"/usr/local/sml/include/SML/basicFunctions.hpp"
#include"/usr/local/pykep/src/core_functions/par2ic.h"
#include"/usr/local/pykep/src/core_functions/ic2par.h"
#include"/usr/local/atom/include/Atom/convertCartesianStateToTwoLineElements.hpp"
#include"/usr/local/pykep/src/lambert_problem.h"

#include"/usr/local/pykep/src/lambert_problem.cpp"

#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/io.hpp>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>
#include <boost/geometry.hpp>
#include <boost/geometry/geometries/linestring.hpp>
#include <boost/geometry/geometries/adapted/boost_array.hpp>
BOOST_GEOMETRY_REGISTER_BOOST_ARRAY_CS(cs::cartesian)


// #include <CoordTopocentric.h>
// #include <CoordGeodetic.h>
// #include <Observer.h>
// #include <SGP4.h>

#include <gsl/gsl_multiroots.h>
#include <gsl/gsl_vector.h>
 	
#include <libsgp4/DateTime.h>
#include <libsgp4/Eci.h>
#include <libsgp4/Globals.h>
#include <libsgp4/SGP4.h>
#include <libsgp4/Tle.h>

#include <Astro/astro.hpp>
#include <SML/sml.hpp>
#include <enne-stage/constants.hpp>

#include <Atom/printFunctions.hpp>

// typedef double Real;
// typedef std::vector< Real > Vector6;

int main()
{
	// std::cout << "test OK enne" << std::endl;
 // std::cout.precision(10); //
 // double tmp[] = {7096137.00,0.0011219, sml::convertDegreesToRadians(92.0316), sml::convertDegreesToRadians(120.6878),sml::convertDegreesToRadians(296.1384), sml::convertDegreesToRadians(239.6546)};
 // boost::array< double, 3> v( tmp, tmp+6 );
 // double c = 398600441000000; 
 // boost::array< double, 3> f = astro::convertKeplerianToCartesianElements(v,c,0.00001);
 // for (int i=0; i<6; i++) std::cout << f[i] << std::endl; 
 // Vector6 cartesianState( 6 );
 //     cartesianState[ 0 ] = -7.1e3;
 //     cartesianState[ 1 ] = 2.7e3;
 //     cartesianState[ 2 ] = 1.3e3;
 //     cartesianState[ 3 ] = -2.5;
 //     cartesianState[ 4 ] = -5.5;
 //     cartesianState[ 5 ] = 5.5;
 // int i;
 //   gsl_vector * fdasfasd = gsl_vector_alloc (3);
 // Tle tle = Tle("UK-DMC 2                ",
 //     "1 35683U 09041C   12289.23158813  .00000484  00000-0  89219-4 0  5863",
 //     "2 35683  98.0221 185.3682 0001499 100.5295 259.6088 14.69819587172294");
 // Tle convertedTle = atom::convertCartesianStateToTwoLineElements< Real, Vector6 >(
 //  	cartesianState, DateTime( ) );
 // std::cout << convertedTle << std::endl; 
 // // Tle g = atom::convertCartesianStateToTwoLineElements(f, DateTime( ) );







// double r_e = enne_stage::kXKMPER *1000;
// std::cout << r_e << std::endl;
double MU_EARTH = enne_stage::ASTRO_ENNE_GM_EARTH;// m^3/s^2

// boost::array< double, 3> r1 = {3126974.99,-6374445.74,28673.59};
boost::array< double, 3> r1 = {6778136,0,0};	
// boost::array< double, 3> r2 = {-3126974.99,6374445.74,-28673.58};
boost::array< double, 3> r2 = {-6978136,0.00007,0};

double r1a = std::sqrt(std::pow(r1[0],2.0)+std::pow(r1[1],2.0)+std::pow(r1[2],2.0));
// double r2a = std::sqrt(std::pow(r2[0],2.0)+std::pow(r2[1],2.0)+std::pow(r2[2],2.0));
// for (int i=0; i<3; i++) std::cout << r1a[i] << std::endl; 

// boost::array< double, 3> v0 = {-254.91197,-83.30107,7485.70674};
boost::array< double, 3> v0 = {0,7668.558733,0};
// boost::array< double, 3> v3 = {254.91197,83.30107,-7485.70673};
boost::array< double, 3> v3 = {0.000000008,-7557.86574,0};

// double v0a = std::sqrt(std::pow(v0[0],2.0)+std::pow(v0[1],2.0)+std::pow(v0[2],2.0));
// double v3a = std::sqrt(std::pow(v3[0],2.0)+std::pow(v3[1],2.0)+std::pow(v3[2],2.0));

double dt_orbit =  2*enne_stage::kPI * std::sqrt(std::pow((r1a+100000),3.0)/MU_EARTH);
// std::cout << dt_orbit << std::endl;

// std::vector<boost::array< double, 3> > dvtotarray;
// boost::array< double, 800> toff;

double start = 500;
const int amount = 100000;
// boost::array< double, 3 > dvtotarray;
// std::vector<boost::array< double, 3> > dvtotarray;
// std::vector<boost::array< double, 3> > dvtotarray;
boost::array< double, amount> toff;

// boost::array< double, amount> dvtotarray;
// std::vector<boost::array< double, amount>> output;

// double dt = dt_orbit*5;
// kep_toolbox::lambert_problem test( r1,r2, dt_orbit, MU_EARTH);
// dvtotarray = test.lambert_problem::get_v1();

for (int j = 0; j < 5; ++j)
{
	boost::array< double, amount> dvtotarray;
	boost::array< double, amount> toff;

	for (int i = start; i < start + amount; ++i)
	{
		const double a = i / 10000.0;
		double dt = dt_orbit * a;
		kep_toolbox::lambert_problem test( r1,r2, dt, MU_EARTH);
		
		int	nmax = test.lambert_problem::get_Nmax();
		// std::cout << nmax << std::endl;
		if (nmax == 0)
		{
		// std::cout << j << std::endl;
		std::vector<boost::array< double, 3> > v1;
		v1 = test.lambert_problem::get_v1();
		
	    double dv1 = boost::geometry::distance(v1[0], v0);

		std::vector<boost::array< double, 3> > v2;
		v2 = test.lambert_problem::get_v2();
	    double dv2 = boost::geometry::distance(v3, v2[0]);

	    double dvtot = dv1+dv2;

		toff[i-start] = a;
		dvtotarray[i-start] = dvtot;
		}
		if (nmax > 0)
		{
		std::vector<boost::array< double, 3> > v1;
		v1 = test.lambert_problem::get_v1();
		
	    double dv1 = boost::geometry::distance(v1[nmax*2-1], v0);

		std::vector<boost::array< double, 3> > v2;
		v2 = test.lambert_problem::get_v2();
	    double dv2 = boost::geometry::distance(v3, v2[nmax*2-1]);

	    double dvtot = dv1+dv2;

		toff[i-start] = a;
		dvtotarray[i-start] = dvtot;

		}
		// else break;
	}
	std::string filename;
	std::ostringstream convert;
	convert << j;
	filename = "../data/dv_" + convert.str();
	filename += ".txt";
  	std::ofstream myfile1;
    myfile1.open (filename);
	for (int i=0; i<amount; i++) myfile1 << dvtotarray[i] << std::endl;
  	myfile1.close();

}
  std::ofstream myfile;
  myfile.open ("../data/i.txt");
  for (int i=0; i<amount; i++) myfile << toff[i] << std::endl;
  myfile.close();
  return 0;
}