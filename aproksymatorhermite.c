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
	if ( stopien == 0)
		return 0;
	else if ( stopien == 1)
		return 0;
	else
		return 4 * dhi(stopien - 1, x) + 2 * x * d2hi(stopien - 1, x) - 2 * (stopien - 1) * d2hi( stopien - 2, x);
}

/* Trzecia pochodna hi */
double
d3hi(int stopien, double x)
{
	if ( stopien == 0)
		return 0;
	else if ( stopien == 1)
		return 0;
	else
		return 6 * d2hi(stopien - 1, x) + 2 * x * d3hi(stopien - 1, x) - 2 * (stopien - 1) * d3hi( stopien - 2, x);
}

void
make_spl(points_t * pts, spline_t * spl)
{
	matrix_t       *eqs= NULL;
	double         *x = pts->x;
	double         *y = pts->y;
	double		a = x[0];
	double		b = x[pts->n - 1];
	int		i, j, k;
	int		nb = pts->n - 3 > 10 ? 10 : pts->n - 3;
  char *nbEnv= getenv( "APROX_BASE_SIZE" );

	if( nbEnv != NULL && atoi( nbEnv ) > 0 )
		nb = atoi( nbEnv );

	eqs = make_matrix(nb, nb + 1);


	for (j = 0; j < nb; j++) {
		for (i = 0; i < nb; i++)
			for (k = 0; k < pts->n; k++)
				add_to_entry_matrix(eqs, j, i, hi( i, x[k]) * hi(j, x[k]));

		for (k = 0; k < pts->n; k++)
			add_to_entry_matrix(eqs, j, nb, y[k] * hi(j, x[k]));
	}


	if (piv_ge_solver(eqs)) {
		spl->n = 0;
		return;
	}
	if (alloc_spl(spl, nb) == 0) {
		for (i = 0; i < spl->n; i++) {
			double xx = spl->x[i] = a + i*(b-a)/(spl->n-1);
			xx+= 10.0*DBL_EPSILON;  // zabezpieczenie przed ulokowaniem punktu w poprzednim przedziale
			spl->f[i] = 0;
			spl->f1[i] = 0;
			spl->f2[i] = 0;
			spl->f3[i] = 0;
			for (k = 0; k < nb; k++) {
				double		ck = get_entry_matrix(eqs, k, nb);
				spl->f[i]  += ck * hi  (k, xx);
				spl->f1[i] += ck * dhi (k, xx);
				spl->f2[i] += ck * d2hi(k, xx);
				spl->f3[i] += ck * d3hi(k, xx);
			}
		}
	}

}
