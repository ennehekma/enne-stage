#include <iostream>
#include <fstream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>
#include <boost/geometry.hpp>
#include <boost/geometry/geometries/linestring.hpp>
#include <boost/geometry/geometries/adapted/boost_array.hpp>
#include <boost/numeric/ublas/io.hpp>
#include <boost/numeric/ublas/vector.hpp>
BOOST_GEOMETRY_REGISTER_BOOST_ARRAY_CS(cs::cartesian)

#include <enne-stage/constants.hpp>

#include"/usr/local/pykep/src/lambert_problem.cpp"


int main()
{
double MU_EARTH = enne_stage::ASTRO_ENNE_GM_EARTH;// m^3/s^2

// 200 km, 30deg, quarter orbit
boost::array< double, 3> r1 = {6778136,0,0};	
boost::array< double, 3> r2 = {0,6043243.047,3489068};
boost::array< double, 3> v0 = {0,7668.558733,0};
boost::array< double, 3> v3 = {-7557.86574,0,0};

// // Hohmann transfer initial conditions.
// boost::array< double, 3> r1 = {6778136,0,0};	
// boost::array< double, 3> r2 = {-6978136,0.00007,0};
// boost::array< double, 3> v0 = {0,7668.558733,0};
// boost::array< double, 3> v3 = {0.000000008,-7557.86574,0};

double r1a = std::sqrt(std::pow(r1[0],2.0)+std::pow(r1[1],2.0)+std::pow(r1[2],2.0));
double dt_orbit =  2*enne_stage::kPI * std::sqrt(std::pow((r1a),3.0)/MU_EARTH);

double start = 500;
const int amount = 100000;

boost::array< double, amount> toff;

std::ofstream myfile;
myfile.open ("../data/i.txt");
for (int i=start; i<start+amount; i++) myfile << i / 10000.0 << std::endl;
myfile.close();

for (int j = 0; j < 5; ++j)
{
	boost::array< double, amount> dvtotarray;
	boost::array< double, amount> toff;

	for (int i = start; i < start + amount; ++i)
	{
		const double a = i / 10000.0;
		double dt = dt_orbit * a;
		kep_toolbox::lambert_problem test( r1,r2, dt, MU_EARTH,0,j);
		int	nmax = test.lambert_problem::get_Nmax();

		if (nmax > 0 && j>0)
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

		if (nmax == 0 && j==0)
		{
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

		if (nmax == 0 && j>0)
		{
		toff[i-start] = a;
		dvtotarray[i-start] = 100000;
		}
		if (nmax == 1 && j>1)
		{
		toff[i-start] = a;
		dvtotarray[i-start] = 100000;
		}
		if (nmax == 2 && j>2)
		{
		toff[i-start] = a;
		dvtotarray[i-start] = 100000;
		}
		if (nmax == 3 && j>3)
		{
		toff[i-start] = a;
		dvtotarray[i-start] = 100000;
		}
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
  
  return 0;
}