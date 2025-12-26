#include <iostream>
#include <cmath>
#include <bits/stdc++.h>
#include <vector>
#include <stdexcept>
using namespace std;

#define MAX 100000
#define SPLIT 1000000

// Turn the int type to a double to gain greater range of double values.
double factorial(double n){
    double integer = 1;
    for (int i = 2; i <= n; i++){
        integer *= i;
    }
    return integer;
}


// Integration 
// https://stackoverflow.com/questions/60005533/composite-simpsons-rule-in-c
template <typename func_type>
double integration(double a, double b, int n, func_type f){
    // Here I am using Simpson's Method for integration.
    
    double h = (b-a)/n; // even split in the formula

    double sum_odds = 0.0;
    for (int i = 1; i < n; i += 2){
        sum_odds += f(a + i*h);
    }
    double sum_evens = 0.0;
    for (int i = 2; i < n; i+=2){
        sum_evens += f(a + i*h);
    };

    return (f(a) + f(b) + 2*sum_evens + 4*sum_odds) * h/3;
}


// n choose r
long long int nCr(int n, int r){
    if (r > n){
        throw std::runtime_error("r cannot be bigger than n!");
    }
    if (r < 0 || n < 0){
        throw std::runtime_error("r or n is negative!");
    }

    // Just coming directly from the binomial coefficient formula.
    double prod = 1;
    for (int i = 1; i <= r; i++){
        prod = prod * (n - r + i)/i;
    }   
    return (long long int)prod;
}

// Here is to check if a double value is actually an integer.
bool isInteger(double z){
    int x = z;
    double val = z - x;
    if (val > 0){
        return false;
    }
    return true;  
}


// Gamma function (Keep this until i figure out a better way to optimise this for value 170! Does around 85. It's to do with the pow function)
// For now, keep the gamma function to only take in real numbers.
bool isNegativeInteger(double z){
    if (z >= 0) return false;
    if (z == floor(z)){
        return true;
    }

    return false;
}

// Sacrificing a little bit of precision with proper functionality.
double gamma(double z){
    if (z == 0 || isNegativeInteger(z)) {
        throw std::invalid_argument("Invalid Input!");
    }

    // Here it's a positive integer number
    if (z ==  floor(z)){
        return factorial(z-1);
    }

    // Bigger non-integer numbers are going to be recursive.
    if (z > 2){
        return (z-1) * gamma(z-1);
    } // smaller numbers need to be raised up 
    else if (z < 1){
        return gamma(z+1) / z;
    }

    auto f = [z](double x){
        return pow(x, z-1) * exp(-x);
    };

    double fact = integration(0, MAX, SPLIT, f);
    return fact;
}


template <size_t N>
double sum(const double (&x)[N]){
    double val = 0.0;
    for (size_t i = 0; i < N; i++){
        val += x[i];
    }
    return val;
} 

// Integer power function
double mypow(double x, int n){
    if (n == 1) return x;
    if (n == 0) return 1;

    int half = n/2;
    if (n % 2 == 1){
        int half2 = half + 1;
        return mypow(x, half) * mypow(x, half2);
    }
    return mypow(x, half) * mypow(x, half);
}


// Mean Function
template <typename T>
double mean(const vector<T>& data){
    int n = data.size();
    double val = 0;
    for (int i = 0; i < n; i++){
        val += data[i];
    }
    double finalval = val / n;
    return finalval;
}










