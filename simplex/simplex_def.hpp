#pragma once

class simplex
    {
public:
    simplex(size_t dimension);
  
    double optimize(double *point, double(*func)(double *, int), double tol);
    double global_optimize(double *point, double(*func)(double *, int), double tol);

private:
    double * alloc_vector(int cols);
    void free_vector(double * vector, int cols);

    double ** alloc_matrix(int rows, int cols);
    void free_matrix(double ** matrix, int rows, int cols);

    double ** make_simplex(double * point, int dim, double * del);
    
    void evaluate_simplex(double ** simplex, int dim, double * fx, double(*func)(double *, int));
    void simplex_extremes(double *fx, int dim, int *ihi, int *ilo, int *inhi);
    void simplex_bearings(double ** simplex, int dim, double * midpoint, double * line, int ihi);

    int update_simplex(double * point, int dim, double * fmax, double * midpoint, 
                       double * line, double scale, double(*func)(double *, int));

    void contract_simplex(double ** simplex, int dim, double * fx, int ilo, double(*func)(double *, int));
    
    int check_tol(double fmax, double fmin, double ftol);

    double amoeba(double *point, int dim, double * del, double(*func)(double *, int), double tol);

    // temperature based update_simplex2 for simulated annealing.
    int update_simplex2(double * point, int dim, double * fmax, double * midpoint, 
                        double * line, double scale, double temperature,  double(*func)(double *, int));
 
    // uniform random distribution generator in range [0,1).
    // if you use C++11 it will use <random> and <chrono>
    // otherwise it use std::rand in <cstdlib> and std::time in <ctime>.
    double urand();
    // gaussian distribution.
    double grand();
    // global optimization using amoeba combined with simulated annealing.
    double amoeba2(double *point, int dim, double * del, double(*func)(double *, int), double tol);

    size_t m_dimension;
    };
