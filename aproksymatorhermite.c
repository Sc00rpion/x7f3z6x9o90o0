#include "makespl.h"
#include "piv_ge_solver.h"
#include "aproksymatorhermite.h"

#include <stdio.h>
#include <stdlib.h>
#include <float.h>

/* Funkcja bazowa -wielomian Hermite'a */
double
hi(int stopien, double x)
{
	if ( stopien == 0)
		return 1;
	else if ( stopien == 1)
		return 2 * x;
	else
		return 2 * x * hi(stopien - 1, x) - 2 * (stopien - 1) * hi( stopien - 2, x);
}

/* Pierwsza pochodna hi */
double
dhi(int stopien, double x)
{
	if ( stopien == 0)
		return 0;
	else if ( stopien == 1)
		return 2;
	else
		return 2 * hi(stopien - 1, x) + 2 * x * dhi(stopien - 1, x) - 2 * (stopien - 1) * dhi(stopien - 2, x);
}

/* Druga pochodna hi */
double
d2hi(int stopien, double x)
{

}

/* Trzecia pochodna hi */
double
d3hi(int stopien, double x)
{

}

void
make_spl(points_t * pts, spline_t * spl)
{

}