#include "simplex.hpp"

#include <stdio.h>
#include <math.h>

int num_eval = 0;

double func(double * p, int dim)
    {
    double x = p[0], y = p[1];
    
    num_eval++;

    return -(y - x -2 * x * x - 2 * x * y - y * y);
    }

int
main()
    {
    // initial guess.
    double p[] = {0.0, 2.0};
    double fmin = 0.0;
    // 2 = 2 dimensional problem.
    simplex optimizer2D (2);
    
    fmin = optimizer2D.optimize(p, func, 1.e-15);
    
    printf("min of f(x,y) = %12.6f\n", fmin);
    printf("at x = %12.6f\n", p[0]);
    printf("at y = %12.6f\n", p[1]);
    printf("# of function evaluation = %5d\n", num_eval);

    return 0;
    }
