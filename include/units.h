#ifndef CUDA_UNITS_H
#define CUDA_UNITS_H

#include <stdarg.h>

// Hardwired constants:
#define AU_METERS 149597870700.0 
#define PC_METERS AU_METERS * (648000.0 / M_PI) //-- A parsec in m
#define MPC_METERS PC_METERS * 10E6
#define SIDEREAL_YEAR_SECONDS 365.256363004*24*3600

#define G_SI 6.67408E-11  //-- the gravitational constant in m^3/(kilogram*s^2)
#define MASS_SUN_KILOGRAMS   (4.0*M_PI*M_PI*AU_METERS*AU_METERS*AU_METERS) \
						   / (G_SI*SIDEREAL_YEAR_SECONDS*SIDEREAL_YEAR_SECONDS)

#define C_SI 299792458.0 //-- the speed of light in m/s
#define NUM_POLARIZATION_STATES 2

typedef double float64_t;

typedef struct {
	float64_t x; 
	float64_t y;
} float64_2_t;

// Mass functions:

typedef struct {
	float64_t msun;
	float64_t kilograms;
	float64_t seconds;
	float64_t meters;
} massUnit_t;

inline float64_t kilogramsToSeconds(const float64_t mass_kilogram)
{
	return mass_kilogram*G_SI / (C_SI*C_SI*C_SI);
}

inline float64_t kilogramsToMeters(const float64_t mass_kilogram)
{
	return mass_kilogram*G_SI / (C_SI*C_SI);
}

inline float64_t msunToKilograms(const float64_t mass_msun)
{
	return mass_msun*MASS_SUN_KILOGRAMS;
}

massUnit_t initMassSolarMass(
	const float64_t mass_msun
) {
	massUnit_t mass = {
		.msun      = mass_msun,
		.kilograms = msunToKilograms   (mass.msun),
		.seconds   = kilogramsToSeconds(msunToKilograms(mass_msun)),
		.meters    = kilogramsToMeters (msunToKilograms(mass_msun))
	};
	
	return mass;
}

massUnit_t scaleMass(
	const massUnit_t    mass, 
	const float64_t scalar
) {
	massUnit_t scaled = {
		.msun      = mass.msun      * scalar,
		.kilograms = mass.kilograms * scalar,
		.seconds   = mass.seconds   * scalar,
		.meters    = mass.meters    * scalar
	};
	
	return scaled;
}

massUnit_t addMasses(
	const massUnit_t mass_1, 
	const massUnit_t mass_2
) {
	massUnit_t sum = {
		.msun      = mass_1.msun      + mass_2.msun,
		.kilograms = mass_1.kilograms + mass_2.kilograms,
		.seconds   = mass_1.seconds   + mass_2.seconds,
		.meters    = mass_1.meters    + mass_2.meters
	};
	
	return sum;
}

massUnit_t subtractMasses(
	const massUnit_t mass_1, 
	const massUnit_t mass_2
) {
	massUnit_t difference = {	
		.msun      = mass_1.msun      - mass_2.msun,
		.kilograms = mass_1.kilograms - mass_2.kilograms,
		.seconds   = mass_1.seconds   - mass_2.seconds,
		.meters    = mass_1.meters    - mass_2.meters
	};
	
	return difference;
}

massUnit_t multiplyMasses(
	const massUnit_t mass_1, 
	const massUnit_t mass_2
) {
	massUnit_t product = {
		.msun      = mass_1.msun      * mass_2.msun,
		.kilograms = mass_1.kilograms * mass_2.kilograms,
		.seconds   = mass_1.seconds   * mass_2.seconds,
		.meters    = mass_1.meters    * mass_2.meters
	};
	
	return product;
}

massUnit_t divideMasses(
	const massUnit_t mass_1, 
	const massUnit_t mass_2
) {
	massUnit_t quotient = {
		.msun      = mass_1.msun      / mass_2.msun,
		.kilograms = mass_1.kilograms / mass_2.kilograms,
		.seconds   = mass_1.seconds   / mass_2.seconds,
		.meters    = mass_1.meters    / mass_2.meters
	};
	
	return quotient;
}

// Distance functions:
typedef struct {
	float64_t Mpc;
	float64_t meters;
} lengthUnit_t;

inline float64_t MpcToMeters(const float64_t length_mpc)
{
	return length_mpc*MPC_METERS;
}

lengthUnit_t initLengthMpc(
	const float64_t length_mpc
) {
	lengthUnit_t length = {
		 .Mpc    = length_mpc,
		 .meters = MpcToMeters(length.Mpc)
	 };
	
	return length;
}

lengthUnit_t scaleLength(
	const lengthUnit_t  length, 
	const float64_t scalar
) {
	lengthUnit_t scaled = {
		.Mpc    = length.Mpc    * scalar,
		.meters = length.meters * scalar
	};
	
	return scaled;
}

// Time Functions:

typedef struct {
	float64_t seconds;
} timeUnit_t;

timeUnit_t initTimeSeconds(
	float64_t time_seconds
	) {
	
	timeUnit_t time = {
		 .seconds = time_seconds
	 };
	
	return time;
}

timeUnit_t _addTimes(
	const timeUnit_t time_1,
	const timeUnit_t time_2
) {
	timeUnit_t sum = {
		.seconds = time_1.seconds + time_2.seconds
	};
	
	return sum;
}

timeUnit_t addTimes(
	const int32_t num_args,
	...
	) {

	va_list valist;
   	timeUnit_t sum = {
		.seconds = 0.0
	};

   va_start(valist, num_args);
   for (int32_t index = 0; index < num_args; index++) 
   {
      sum = _addTimes(sum, (timeUnit_t)va_arg(valist, timeUnit_t));
   }
   va_end(valist);

   return sum;
}

timeUnit_t multiplyTimes(
	const timeUnit_t time_1,
	const timeUnit_t time_2
) {
	timeUnit_t product = {
		.seconds = time_1.seconds * time_2.seconds
	};
	
	return product;
}

timeUnit_t scaleTime(
	const timeUnit_t time_1,
	const float64_t  scale
) {
	timeUnit_t scaled = {
		.seconds = time_1.seconds * scale
	};
	
	return scaled;
}

// Frequency Functions:
typedef struct {
	float64_t hertz;
} frequencyUnit_t;

frequencyUnit_t initFrequencyHertz(
	float64_t frequency_hertz
	) {
	
	frequencyUnit_t frequency = {
		 .hertz = frequency_hertz
	 };
	
	return frequency;
}

// Angle Functions:
typedef struct {
	float64_t radians;
} angularUnit_t;

angularUnit_t initAngleRadians(
	float64_t angle_radians
	) {
	
	angularUnit_t angle = {
		 .radians = angle_radians
	 };
	
	return angle;
}


#endif