#ifndef CUDA_UNITS_H
#define CUDA_UNITS_H

#include <stdarg.h>

// Hardwired constants.
#define PC_METERS 3.0856775807f*1e16f //-- a parsec in m
#define MASS_SUN_KILOGRAMS 1.988409902147041637325262574352366540e30 // Mass of the sun in kilogram
#define MASS_SUN_SECONDS 4.925491025543575903411922162094833998e-6 // Geometrized solar mass, s.
#define MASS_SUN_METERS 1.476625061404649406193430731479084713e3 // Geometrized solar mass, m.
#define G_SI 6.67430e-11  //-- the gravitational constant in m^3/(kilogram*s^2)
#define C_SI 299792458.0f //-- the speed of light in m/s
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
} mass_t;

inline float64_t kilogramsToSeconds(const float64_t mass_kilogram)
{
	return mass_kilogram*G_SI / pow(C_SI, 3.0);
}

inline float64_t kilogramsToMeters(const float64_t mass_kilogram)
{
	return mass_kilogram*G_SI / (C_SI*C_SI);
}

inline float64_t msunToKilograms(const float64_t mass_msun)
{
	return mass_msun*MASS_SUN_KILOGRAMS;
}

mass_t initMassSolarMass(
	const float64_t mass_msun
) {
	mass_t mass = {
		.msun      = mass_msun,
		.kilograms = msunToKilograms  (mass.msun),
		.seconds   = kilogramsToSeconds(mass.kilograms),
		.meters    = kilogramsToMeters (mass.meters)
	};
	
	return mass;
}

mass_t scaleMass(
	const mass_t    mass, 
	const float64_t scalar
) {
	mass_t scaled = {
		.msun      = mass.msun      * scalar,
		.kilograms = mass.kilograms * scalar,
		.seconds   = mass.seconds   * scalar,
		.meters    = mass.meters    * scalar
	};
	
	return scaled;
}

mass_t addMasses(
	const mass_t mass_1, 
	const mass_t mass_2
) {
	mass_t sum = {
		.msun      = mass_1.msun      + mass_2.msun,
		.kilograms = mass_1.kilograms + mass_2.kilograms,
		.seconds   = mass_1.seconds   + mass_2.seconds,
		.meters    = mass_1.meters    + mass_2.meters
	};
	
	return sum;
}

mass_t subtractMasses(
	const mass_t mass_1, 
	const mass_t mass_2
) {
	mass_t difference = {	
		.msun      = mass_1.msun      - mass_2.msun,
		.kilograms = mass_1.kilograms - mass_2.kilograms,
		.seconds   = mass_1.seconds   - mass_2.seconds,
		.meters    = mass_1.meters    - mass_2.meters
	};
	
	return difference;
}

mass_t multiplyMasses(
	const mass_t mass_1, 
	const mass_t mass_2
) {
	mass_t product = {
		.msun      = mass_1.msun      * mass_2.msun,
		.kilograms = mass_1.kilograms * mass_2.kilograms,
		.seconds   = mass_1.seconds   * mass_2.seconds,
		.meters    = mass_1.meters    * mass_2.meters
	};
	
	return product;
}

mass_t divideMasses(
	const mass_t mass_1, 
	const mass_t mass_2
) {
	mass_t quotient = {
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
} length_t;

inline float64_t MpcToMeters(const float64_t length_mpc)
{
	return length_mpc*PC_METERS*10E6f;
}

length_t initLengthMpc(
	const float64_t length_mpc
) {
	length_t length = {
		 .Mpc    = length_mpc,
		 .meters = MpcToMeters(length.Mpc)
	 };
	
	return length;
}

length_t scaleLength(
	const length_t  length, 
	const float64_t scalar
) {
	length_t scaled = {
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

#endif