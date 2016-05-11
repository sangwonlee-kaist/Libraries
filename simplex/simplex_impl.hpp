#pragma once

#ifdef SIMPLEX_DEBUG
    #include <iostream>
    #define DEBUG(x) std::cout << #x << " = " << (x) << std::endl;
    #define DEBUG2(x, y) std::cout << #x << " = " << (x) << ", "  #y << " = " << (y) << std::endl;
#else
    #define DEBUG(x)
    #define DEBUG2(x, y)
#endif

simplex::simplex(size_t dimension) : m_dimension (dimension)
    {
    }

double 
simplex::optimize(double* point, double(*func)(double*, int), double tol)
    {
    int i;
    int dim = m_dimension;

    double* del = alloc_vector(dim);
    // 10 % of initial position.
    for (i = 0; i < dim; i++)
        {
        del[i] = 0.1 * point[i];    
        DEBUG(del[i])
        }

    for (i = 0; i < dim; i++)
        {
        del[i] = (del[i] < 1.e-5 ? 1.e-5 : del[i]);
        DEBUG(del[i])
        }

    double fmin = amoeba(point, dim, del, func, tol);

    free_vector(del, dim);

    return fmin;
    }

double* 
simplex::alloc_vector(int cols)
    {
    return (double* )std::malloc(sizeof(double) * cols);
    }

void 
simplex::free_vector(double* vector, int cols)
    {
    std::free(vector);
    }

double** 
simplex::alloc_matrix(int rows, int cols)
    {
    int i;
    double** matrix = (double**)std::malloc(sizeof(double* ) * rows);
    for (i = 0; i < rows; i++)
        matrix[i] = alloc_vector(cols);
    return matrix;
    }

void 
simplex::free_matrix(double** matrix, int rows, int cols)
    {
    int i;
    for (i = 0; i < rows; i++)
        free_vector(matrix[i], cols);
    std::free(matrix);
    }

double** 
simplex::make_simplex(double* point, int dim, double* del)
    {
    int i, j;
    double** simplex = alloc_matrix(dim + 1, dim);
    for (i = 0; i < dim + 1; i++)
        for (j = 0; j < dim; j++)
            simplex[i][j] = point[j];
    for (i = 0; i < dim; i++)
        simplex[i][i] += del[i];
    return simplex;
    }   

void 
simplex::evaluate_simplex(double** simplex, int dim, double* fx, 
    double(*func)(double* , int))
    {
    // simplex::evaluate_simplex(): 
    // evaluate function value for each simpelx point.
    // calculated function values are stored in double array fx.
    for (int i = 0; i < dim + 1; i++)
        fx[i] = (*func)(simplex[i], dim);
    }

void 
simplex::simplex_extremes(double* fx, int dim, int* ihi, int* ilo,
    int* inhi)
    {
    // simlex::simplex_extreams(): 
    // find order of between function values of simplex points.
    // find maximum, next maximum, and minimum then save their indices.
    // ihi:  index of higher value.
    // inhi: index of next higher value.
    // ilo:  index of lower value.
    int i;
    if (fx[0] > fx[1])
        {
        *ihi = 0; *ilo = *inhi = 1;
        }
    else
        {
        *ihi = 1; *ilo = *inhi = 0;
        }
    for (i = 2; i < dim + 1; i++)
        if (fx[i] <= fx[*ilo])
            *ilo = i;
        else if (fx[i] > fx[*ihi])
            {
            *inhi = *ihi; *ihi = i;
            }
        else if (fx[i] > fx[*inhi])
            *inhi = i;
    }

void 
simplex::simplex_bearings(double** simplex, int dim,
    double* midpoint, double* line, int ihi)
    {
    // simplex::simplex_bearings(): 
    // make two vectors used in simplex trial moves.
    // calculate middle(center) point of simplex
    // and a lines through higher point and moddle point.
    int i, j;
    for (j = 0; j < dim; j++)
        midpoint[j] = 0.0;
    for (i = 0; i < dim + 1; i++)
        if (i != ihi)
            for (j = 0; j < dim; j++)
                midpoint[j] += simplex[i][j];
    for (j = 0; j < dim; j++)
        {
        midpoint[j] /= dim;
        line[j] = simplex[ihi][j] - midpoint[j];
        }
    }

int 
simplex::update_simplex(double* point, int dim, double* fmax,
    double* midpoint, double* line, double scale,
    double(*func)(double* , int))
    {
    // simplex::update_simplex():
    // make trial points and change to maximum point if
    // new point is lower than previous maximum point
    int i, update = 0; double* next = alloc_vector(dim), fx;
    for (i = 0; i < dim; i++)
        next[i] = midpoint[i] + scale * line[i];
    fx = (*func)(next, dim);
    if (fx < *fmax)
        {
        for (i = 0; i < dim; i++) point[i] = next[i];
        *fmax = fx; update = 1;
        }
    free_vector(next, dim);
    return update;
    }

void 
simplex::contract_simplex(double** simplex, int dim,
    double* fx, int ilo,
    double(*func)(double* , int))
    {
    int i, j;
    for (i = 0; i < dim + 1; i++)
        if (i != ilo)
            {
            for (j = 0; j < dim; j++)
                simplex[i][j] = (simplex[ilo][j] + simplex[i][j])*0.5;
            fx[i] = (*func)(simplex[i], dim);
            }
    }   

int 
simplex::check_tol(double fmax, double fmin, double ftol)
    {
    const double ZEPS = std::numeric_limits<double>::epsilon();
    double delta = std::fabs(fmax - fmin);
    double accuracy = (std::fabs(fmax) + fabs(fmin)) * ftol;
    return (delta < (accuracy + ZEPS));
    }

double 
simplex::amoeba(double* point, int dim, double* del,
    double(*func)(double* , int),
    double tol)
    {
    int ihi, ilo, inhi, j;
    double fmin;
    double* fx = alloc_vector(dim + 1);
    double* midpoint = alloc_vector(dim);
    double* line = alloc_vector(dim);
    double** simplex = make_simplex(point, dim, del);
    evaluate_simplex(simplex, dim, fx, func);

    size_t max_iter = 10000;
    size_t iter = 0;

    while (true)
        {   
        simplex_extremes(fx, dim, &ihi, &ilo, &inhi);
        simplex_bearings(simplex, dim, midpoint, line, ihi);
        if (check_tol(fx[ihi], fx[ilo], tol)) break;
        update_simplex(simplex[ihi], dim, &fx[ihi],
            midpoint, line, -1.0, func);
        if (fx[ihi] < fx[ilo])
            update_simplex(simplex[ihi], dim, &fx[ihi],
            midpoint, line, -2.0, func);
        else if (fx[ihi] >= fx[inhi])
            if (!update_simplex(simplex[ihi], dim, &fx[ihi],
                midpoint, line, 0.5, func))
                contract_simplex(simplex, dim, fx, ilo, func);
        if (++iter >= max_iter)
            throw std::runtime_error ("error: simplex exceed maximum iteration.");
        DEBUG2(fx[ihi], fx[ilo])
        }

    for (j = 0; j < dim; j++)
        point[j] = simplex[ilo][j];
    fmin = fx[ilo];
    free_vector(fx, dim);
    free_vector(midpoint, dim);
    free_vector(line, dim);
    free_matrix(simplex, dim + 1, dim);
    return fmin;
    }
