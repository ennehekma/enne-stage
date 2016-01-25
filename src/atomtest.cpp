#include <iostream>
#include <fstream>
#include <vector>
#include"/usr/local/astro/include/Astro/orbitalElementConversions.hpp"
#include"/usr/local/sml/include/SML/basicFunctions.hpp"
#include"/usr/local/pykep/src/core_functions/par2ic.h"
#include"/usr/local/pykep/src/core_functions/ic2par.h"
#include"/usr/local/atom/include/Atom/convertCartesianStateToTwoLineElements.hpp"

#include <iostream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

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

#include <Atom/printFunctions.hpp>

typedef double Real;
typedef std::vector< Real > Vector6;

int main()
{
std::cout.precision(10); //
double tmp[] = {7096137.00,0.0011219, sml::convertDegreesToRadians(92.0316), sml::convertDegreesToRadians(120.6878),sml::convertDegreesToRadians(296.1384), sml::convertDegreesToRadians(239.6546)};
std::vector<double> v( tmp, tmp+6 );

double c = 398600441000000; 
std::vector<double> f = astro::convertKeplerianToCartesianElements(v,c,0.00001);
for (int i=0; i<6; i++) std::cout << f[i] << std::endl; 

Vector6 cartesianState( 6 );
    cartesianState[ 0 ] = -7.1e3;
    cartesianState[ 1 ] = 2.7e3;
    cartesianState[ 2 ] = 1.3e3;
    cartesianState[ 3 ] = -2.5;
    cartesianState[ 4 ] = -5.5;
    cartesianState[ 5 ] = 5.5;

int i;
  gsl_vector * fdasfasd = gsl_vector_alloc (3);

Tle tle = Tle("UK-DMC 2                ",
    "1 35683U 09041C   12289.23158813  .00000484  00000-0  89219-4 0  5863",
    "2 35683  98.0221 185.3682 0001499 100.5295 259.6088 14.69819587172294");

Tle convertedTle = atom::convertCartesianStateToTwoLineElements< Real, Vector6 >(
 	cartesianState, DateTime( ) );

std::cout << convertedTle << std::endl; 


// Tle g = atom::convertCartesianStateToTwoLineElements(f, DateTime( ) );

}


