#pragma once

class simplex
    {
public:
    simplex(size_t dimension);
    double optimize(double* point, double(*func)(double*, int), double tol);
private:
    double* alloc_vector(int cols);
    void free_vector(double* vector, int cols);

    double** alloc_matrix(int rows, int cols);
    void free_matrix(double** matrix, int rows, int cols);

    double** make_simplex(double* point, int dim, double* del);

    void evaluate_simplex(double** simplex, int dim, double* fx, double(*func)(double* , int));
    void simplex_extremes(double* fx, int dim, int* ihi, int* ilo, int* inhi);
    void simplex_bearings(double** simplex, int dim, double* midpoint, double* line, int ihi);
    void contract_simplex(double** simplex, int dim, double* fx, int ilo, double(*func)(double* , int));

    int update_simplex(double* point, int dim, double* fmax, double* midpoint,
                       double* line, double scale, double(*func)(double* , int));
    int check_tol(double fmax, double fmin, double ftol);

    double amoeba(double* point, int dim, double* del, double(*func)(double* , int), double tol);

    size_t m_dimension;
    };
