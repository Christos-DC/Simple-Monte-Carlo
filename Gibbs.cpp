//#include <iostream>
//#include "mathfun.hpp"
#include "distributions.hpp"
//using namespace std;

// Perform the following Bayesian Problem in C++
// X | p, N is Binomial(N,p)
// p||x,N is Beta(x+2, N-x+2)
// N - x | x,p is Poisson(16(1-p))

int main(){
    // Initial values and vectors.
    vector<int> Xvec = {2};
    vector<double> pvec = {0.5};
    vector<int> Nvec = {3};

    int iterations, burnin;
    cout << "Enter the number of iterations: ";
    cin >> iterations; 

    cout << "Enter the burn-in period: ";
    cin >> burnin; 

    if (burnin >= iterations){
        cout << "Retry again due to invalid input values!" << endl;
        return 0;
    }
    
    int Xval, Nval;
    double pval;
    for (int i = 1; i <= iterations; i++){
        Xval = rbinom(pvec.back(), Nvec.back());
        pval = rbeta(Xval + 2, Nvec.back() - Xval + 4);
        Nval = rpois(16*(1 - pval)) + Xval; // Poisson distribution function is very slow

        Xvec.push_back(Xval);
        pvec.push_back(pval);
        Nvec.push_back(Nval);

        if (i % 1000 == 0){
            cout << "The " << i << "th iteration is done." << endl;
        }
    }

    Xvec.erase(Xvec.begin(), Xvec.begin() + burnin);
    pvec.erase(pvec.begin(), pvec.begin() + burnin);
    Nvec.erase(Nvec.begin(), Nvec.begin() + burnin);

    cout << "Mean value of X is : " << mean(Xvec) << endl;
    cout << "Mean value of p is : " << mean(pvec) << endl;
    cout << "Mean value of N is : " << mean(Nvec) << endl;

    return 0;
}