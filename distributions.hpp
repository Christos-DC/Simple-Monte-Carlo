//#include "mathfun.hpp"
#include <cmath>
#include <cstdlib>
#include <ctime>
#include "mathfun.hpp"
#include <random>

#define INTPART 1000

// binomial pmf function
double dbinom(int x, double p, int n, bool if_log = false){
    long double val = nCr(n,x) * pow(p, x) * pow(1-p, n - x);
    if (if_log){
        return log(val);
    }
    return val;
}

// cumulative binomial distribution
double pbinom(int x, double p, int n, bool lower_tail = true, bool if_log = false){
    double val = 0.0L;
    for (int i = 0; i <= x; i++){
        val = val + dbinom(i, p, n);
    }
    if (!lower_tail) val = 1 - val;
    if (if_log) val = log(val);

    return val;
}

// exponential pdf function.
double dexp(double x, double rate, bool if_log = false){
    if (x < 0) return 0;

    long double val = rate * exp(-rate * x);
    if (if_log) val = log(val);

    return val;
}

double pexp(double q, double rate, bool lower_tail = true, bool if_log = false){
    // Turn this integration into a lambda function.
    if (q < 0) return 0;

    auto f = [rate](double q){
        return dexp(q,rate);
    };
    double result = integration(0, q, INTPART, f);
    
    if (!lower_tail) result = 1 - result;
    if (if_log) result = log(result);

    return result;
}


// Gamma distribution
double dgamma(double x, double alpha, double beta, bool if_log = false){
    double coeff = pow(beta, alpha) / gamma(alpha);
    double val = coeff * pow(x, alpha - 1) * exp(-beta * x);
    if (if_log) val = log(val);
    return val;
}


double pgamma(double q, double alpha, double beta, bool lower_tail = true, bool if_log = false){
    if (q < 0) return 0;

    auto f = [alpha, beta](double q){
        return dgamma(q, alpha, beta);
    };
    double result = integration(0, q, INTPART, f);

    if (!lower_tail) result = 1 - result;
    if (if_log) result = log(result);

    return result;
}


// Poisson distribution
double dpois(int x, double rate, bool if_log = false){
    int fact;
    if (x <= 20){
        fact = factorial(x);
    } else {
        fact = gamma(x+1);
    }

    double val = pow(rate, x) * exp(-rate) / fact;
    if (if_log) val = log(val);
    return val;
}


double ppois(int x, double rate, bool lower_tail = true, bool if_log = false){
    double val = 0;
    for (int i = 0; i <= x; i++){
        val += dpois(i,rate);
    }
    if (!lower_tail) val = 1 - val;
    if (if_log) val = log(val);
    return val;
}
// This works but the computation speed is quite slow for even small values of x. It might have to be the combination of the gamma function and summation of poisson variables.
// May need to implement checks to not compute high values.
// Need to work on handling these numbers better.


// Beta distribution
double dbeta(double x, double alpha, double beta, bool if_log = false){
    double Beta = gamma(alpha) * gamma(beta) / gamma(alpha + beta);
    double val = pow(x, alpha - 1) * pow(1-x, beta -1) / Beta;
    if (if_log) val = log(val);
    return val;
}


double pbeta(double q, double alpha, double beta, bool lower_tail = true, bool if_log = false){
    if (q < 0) return 0;
    if (q >= 1) return 1;

    auto f = [alpha, beta](double q){
        return dbeta(q, alpha, beta);
    };
    double result = integration(0, q, INTPART, f);

    if (!lower_tail) result = 1 - result;
    if (if_log) result = log(result);

    return result;
}



// Random generating functions for the major distribution functions.
// Mersenne Twister algorithm in c++ to replace the rand() functions. they are in the <random> library. https://youtu.be/oW6iuFbwPDg?si=vZwNhTdfA8SbZuWN
int rbinom(double p, int n){
    random_device rd;
    int val = 0;
    double randprob;
    for (int i = 0; i < n; i++){
        uniform_int_distribution<int> dist(0, n);
        randprob = dist(rd) / (n * 1.0); // It does integer division without having the 1.0.
        if (randprob <= p) val++;
    }
    return val;
}


double rexp(double rate){
    random_device rd;
    uniform_real_distribution<double> dist(0.0,1.0);
    double uniform = dist(rd);

    if (uniform == 0.0) return 0;

    double result = -1/rate * log(uniform);
    return result;
}


int rpois(double rate){
    // Using a Poisson Process and Exponential waiting times to calculate the number of iterations.
    int N = 0;
    double S = 0, E;
    random_device rd;
    while (S < 1){
        uniform_real_distribution<double> dist(0.0,1.0);
        double uniform = dist(rd);
        E = - log(uniform) / rate;
        N++;
        S = S + E;
    }
    return N;
}


double rgamma(double alpha, double beta){
    int n = ceil(alpha);

    // Simulating the lower and upper bounds to account for non-integer values of alpha.
    double lowerbound = 0;
    for (int i = 1; i <= (n-1); i++){
        lowerbound = lowerbound + rexp(beta);
    }
    double upperbound = lowerbound + rexp(beta);

    // Setting up the weighted averages of simulated values.
    double lowdistance = alpha - (n - 1);
    double upperdistance = n - alpha;

    double simval = (1-lowdistance) * lowerbound + (1-upperdistance) * upperbound;

    return simval;
}



double rbeta(double alpha, double beta){
    double x1 = rgamma(alpha, 1);
    double x2 = rgamma(beta, 1);
    double simval = x1 / (x1 + x2);

    return simval;
}

