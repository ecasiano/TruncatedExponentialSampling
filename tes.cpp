//
//  truncated_exponential.cpp
//  truncated exponential sampling
//
//  Created by Emanuel Casiano-Diaz on 8/22/20.
//  Copyright © 2020 Emanuel Casiano-Díaz. All rights reserved.
//

// #include "pimc.hpp"
// #include "cxxopts.hpp"
// #include <assert.h>
// #include "uuid.hpp"

// #ifndef pimc_hpp
// #define pimc_hpp

#include <stdio.h>
#include<iostream>
#include<vector>
#include<boost/random.hpp>
#include<cmath>
#include<chrono>
#include<iomanip>  // for std::setprecision
#include <fstream>
#include <cstdlib> // for exit function
// #include "RNG.h"
#include<sstream>
#include<string.h>

#include <boost/math/special_functions/lambert_w.hpp>
using boost::math::lambert_w0;
using boost::math::lambert_wm1;

using namespace std;
using namespace std::chrono;

double infer_tau(double x,double a, double b, double c){

    double y,Z,P,arg,atol,tau,A;

    // Compute normalization of truncated exponential dist.
    Z = (1/c) * (exp(c*(b-a)) - 1) - (b - a);
    
    //
    y = Z*x - (1/c)*exp(c*(b-a)) - a;

    //
    A = -(1/c)*exp(c*b);

    // Determine LambertW branch & compute tau
    arg = max(-1/exp(1), A*c*exp(c*y));
    if (c < 0){ // k = 0 branch
        tau = (1/c)*lambert_w0(arg)-y;
    }
    else {      // k = -1 branch
        tau = (1/c)*lambert_wm1(arg)-y;
    }

    // Check with specific x values
    // double F;
    // F = (1/Z) * ((1/c) * (exp(c*(b-a))-exp(c*(b-tau)))-(tau-a));
    // atol = 1.0e-10;
    // assert(abs(y-(A*exp(-c*tau)-tau)) <= atol);
    // assert(abs(F-x) <= atol);
    // assert(a-atol <= tau <= b+atol);

    return tau;
}

/*----------------------------------------------------------------------------*/

double infer_tau2(double x,double a, double b, double c){

    double tau,Z,atol;

    /* ---- */
    // Sample the new time of the worm end from truncated exponential dist.
    /*:::::::::::::::::::: Truncated Exponential RVS :::::::::::::::::::::::::*/
    Z = 1.0 - exp(-(-c)*(b-a));
    tau = a - log(1.0-Z*x)  / (-c);
    /*::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::*/

    return tau;
}

/*----------------------------------------------------------------------------*/

double shifted_infer_tau(double x,double a, double b, double c){

    double y,Z,P,arg,atol,tau,A,a_pre_shift,b_pre_shift,inverse_c;

    // Save original values of the upper bounds
    a_pre_shift = a;
    b_pre_shift = b;

    // Shift lower and upper bounds
    a -= b;
    b -= b;

    // Compute normalization of truncated exponential dist.
    inverse_c = 1/c;
    Z = inverse_c * (exp(c*(b-a)) - 1) - (b - a);
    
    //
    y = Z*x - inverse_c*exp(c*(b-a)) - a;

    //
    A = -inverse_c;

    // Determine LambertW branch & compute tau
    arg = max(-1/exp(1), A*c*exp(c*y));
    if (c < 0){ // k = 0 branch
        tau = inverse_c*lambert_w0(arg)-y;
    }
    else {      // k = -1 branch
        tau = inverse_c*lambert_wm1(arg)-y;
    }

    tau += b_pre_shift;
    // cout << tau << endl;

    // Check with specific x values
    // double F;
    // F = (1/Z) * ((1/c) * (exp(c*(b-a))-exp(c*(b-tau)))-(tau-a));
    // atol = 1.0e-10;
    // assert(abs(y-(A*exp(-c*tau)-tau)) <= atol);
    // assert(abs(F-x) <= atol);
    // assert(a-atol <= tau);
    // assert(tau <= b+atol);

    return tau;
}

/*----------------------------------------------------------------------------*/
// Main
int main(int argc, char** argv){

    volatile double a,a_new,b,b_new,c,x,tau1,tau2,arg1,arg2,arg3,arg0;
    auto constexpr num_samples = 5'000'000;

    int w0_ctr = 0;
    int wm1_ctr = 0;

    // set parameters
    a = 0.1;  // lower bound
    b = 1.3;  // upper bound
    c = -0.5; // exponential decay

    // To collect samples
    vector<double> samples(num_samples,0);
    vector<double> samples2(num_samples,0);

    // Initialize a Mersenne Twister RNG
    int seed = 1;
    boost::random::mt19937 rng(seed);

    // Create a uniform distribution with support: [0.0,1.0)
    boost::random::uniform_real_distribution<double> rnum(0.0,1.0);

    // Time execution of exp(x)
    arg1 = 1.0;
    auto start = high_resolution_clock::now();
    for (int i=0; i < num_samples; i++){
        exp(arg1);
    }
    auto end = high_resolution_clock::now();
    auto elapsed_time = duration_cast<nanoseconds>(end - start);
    double duration = elapsed_time.count() * 1e-9;

    cout << scientific;

    cout << endl;
    cout << "subroutine         | Time per call (seconds)" << endl;
    cout << "------------------------------------" << endl;
    cout << "exp()              | " << duration/num_samples <<endl;

    // Time execution of log()"
    arg0 = 1.2;
    start = high_resolution_clock::now();
    for (int i=0; i < num_samples; i++){
        x = rnum(rng);
        log(x);
    }
     end = high_resolution_clock::now();
     elapsed_time = duration_cast<nanoseconds>(end - start);
     duration = elapsed_time.count() * 1e-9;

    cout << "------------------------------------" << endl;
    cout << "log()              | " << duration/num_samples <<endl;

    // Time execution of max()"
    start = high_resolution_clock::now();
    for (int i=0; i < num_samples; i++){
        max(arg0,arg1);
    }
     end = high_resolution_clock::now();
     elapsed_time = duration_cast<nanoseconds>(end - start);
     duration = elapsed_time.count() * 1e-9;

    cout << "------------------------------------" << endl;
    cout << "max()              | " << duration/num_samples <<endl;

    // Time execution of rnum(rng) "i.e rand()"
    float randsum =0;
    start = high_resolution_clock::now();
    for (int i=0; i < num_samples; i++){
        x = rnum(rng);

        randsum += x;
    }
    // cout << "randsum: " << randsum << endl;

     end = high_resolution_clock::now();
     elapsed_time = duration_cast<nanoseconds>(end - start);
     duration = elapsed_time.count() * 1e-9;

    cout << "------------------------------------" << endl;
    cout << "rand()             | " << duration/num_samples <<endl;

    // Time execution of exp(x1-x2) "i.e rand()"
    randsum =0;
    double x1,x2;
    start = high_resolution_clock::now();
    for (int i=0; i < num_samples; i++){
        x1 = rnum(rng);
        x2 = rnum(rng);

        randsum += exp(x1-x2);
    }
    cout << "randsum: " << randsum << endl;

     end = high_resolution_clock::now();
     elapsed_time = duration_cast<nanoseconds>(end - start);
     duration = elapsed_time.count() * 1e-9;

    cout << "------------------------------------" << endl;
    cout << "exp(x1-x2)         | " << duration/num_samples <<endl;

    // Time execution of lambert_w0"
    arg2 = 10.0;
    start = high_resolution_clock::now();
    for (int i=0; i < num_samples; i++){
        lambert_w0(arg2);
    }
     end = high_resolution_clock::now();
     elapsed_time = duration_cast<nanoseconds>(end - start);
     duration = elapsed_time.count() * 1e-9;

    cout << "------------------------------------" << endl;
    cout << "lambert_w0         | " << duration/num_samples <<endl;

    // Time execution of lambert_wm1"
    arg3 = -0.2;
    float arg3sum = 0;
    start = high_resolution_clock::now();
    for (int i=0; i < num_samples; i++){
        lambert_wm1(arg3);

        arg3sum+=arg3;

    }
     end = high_resolution_clock::now();
     elapsed_time = duration_cast<nanoseconds>(end - start);
     duration = elapsed_time.count() * 1e-9;

    cout << "------------------------------------" << endl;
    cout << "lambert_wm1        | " << duration/num_samples <<endl;

    float sum;
    sum = 0;
    // Time execution of sampling tau1 and tau2"
    start = high_resolution_clock::now();
    for (int i=0;i<num_samples;i++){
        /* sample a random x value from U(0,1)*/
        x = rnum(rng);
        tau1 = infer_tau(x,a,b,c);
        // samples[i] = tau1;

        x = rnum(rng);
        a_new = tau1;
        tau2 = infer_tau2(x,a_new,b,c);

        // samples2[i] = tau2;
        sum+=tau1; 
    }
    cout << "sum: " << sum << endl;
    
     end = high_resolution_clock::now();
     elapsed_time = duration_cast<nanoseconds>(end - start);
     duration = elapsed_time.count() * 1e-9;

    cout << "------------------------------------" << endl;
    cout << "sampling tau1,tau2 | " << duration/num_samples <<endl;
    cout << endl;

    sum = 0;
    // Time execution of sampling shifted tau1 and tau2"
    start = high_resolution_clock::now();
    for (int i=0;i<num_samples;i++){
        /* sample a random x value from U(0,1)*/
        x = rnum(rng);
        tau1 = shifted_infer_tau(x,a,b,c);
        samples[i] = tau1;

        x = rnum(rng);
        a_new = tau1;
        tau2 = infer_tau2(x,a_new,b,c);

        samples2[i] = tau2;
        sum+=tau1; 
    }
    cout << "sum: " << sum << endl;
    
     end = high_resolution_clock::now();
     elapsed_time = duration_cast<nanoseconds>(end - start);
     duration = elapsed_time.count() * 1e-9;

    cout << "------------------------------------" << endl;
    cout << "shifted  tau1,tau2 | " << duration/num_samples <<endl;
    cout << endl;

    // Create file name
    string filename,filename2;
    filename=to_string(a)+"_"+to_string(b)+"_"+to_string(c)+
    "_samples.dat";

    filename2=to_string(a)+"_"+to_string(b)+"_"+to_string(c)+
    "_samples2.dat";

    // Open files
    ofstream samples_file,samples_file2;
    samples_file.open(filename);
    samples_file2.open(filename2);


    // Write sampled numbers to file
    for (int i=0; i<samples.size(); i++){
        samples_file<<fixed<<setprecision(17)<<samples[i]<<endl;
        samples_file2<<fixed<<setprecision(17)<<samples2[i]<<endl;
    }

    // Close file
    samples_file.close();
    samples_file2.close();

    return 0;
}