#include <iostream>
#include <random>
#include <string>
#include <vector>
#include <fstream>
#include <Eigen/Dense>
#include <math.h>

using Eigen::MatrixXd;

template<class T>
void output(std::vector<T> v, std::string filename){
    /*
    write a vector to a .txt file and store it in the Data subdirectory 
    */

    std::ofstream ofs;
    ofs.open("./Data/" + filename + ".txt");
    
    for(int i=0; i<v.size(); i++){
        ofs << v[i] << '\n';
    }
    
    ofs.close();
}

class RNG{

    private:
    std::mt19937 gen;

    public:

    RNG() : gen((std::random_device())()) {}
    
    double gen_norm(double mean, double variance){
        /*
        Sample a normal distribution with given mean and variance
        */
        std::normal_distribution<double> d(mean, sqrt(variance)); 
        return d(gen);
    }

    double gen_uniform(double min, double max){
        /*
        sample a uniform distribution
        */

        std::uniform_real_distribution<double> d(min, max);
        return d(gen);
    }

    std::vector<double> create_norm_vector(int n, double mean, double variance){
        /*
        Creates a vector of N(mean, variance) random variables
        */
        std::vector<double> data;

        for(int i=0; i<n; i++){
            data.push_back(gen_norm(mean, variance));
        }

        return data;
    }
};

double N(double x){
    /*
    Cumulative normal function
    */
    return 0.5*(1+erf(x/sqrt(2)));
}

double N_prime(double x){
    /*
    Derivative of N(x)
    */
    return exp(-0.5*x*x)/sqrt(2*M_PI);
}

class Option{

    public:
    double S, K, r, d, T, sigma;

    Option(double S_, double K_, double r_, double d_, double T_, double sigma_) : S(S_), K(K_), r(r_), d(d_), T(T_), sigma(sigma_){}

    virtual double payout(double x){
        /*
        A generic payout function that will later be replaced by a class for a specific option type
        */
        return 0;
    }

    double d1(double S_, double K_, double r_, double d_, double T_, double sigma_){
        return (log(S_/K_)+(r_-d_+0.5*sigma_*sigma_)*T_)/(sigma_*sqrt(T_));
        
    }

    double d2(double S_, double K_, double r_, double d_, double T_, double sigma_){
        return (log(S_/K_)+(r_-d_-0.5*sigma_*sigma_)*T_)/(sigma_*sqrt(T_));
    }

    double MC_SDE_underlying_price(double z){
        /*
        Return the end value of a lognormally distributed path, as determined by solving the SDE. 
        z should be a N(0,1) distributed random variable, as calculated by the RNG class
        */
        return S*exp((r-d)*T-0.5*sigma*sigma*T+sigma*z*sqrt(T));
    }

    double MC_price(std::vector<double> rand_vals){
        /*
        Return the price of a call using the SDE price as calculated above. 
        rand_vals should contain a number of N(0,1) random variables equal to the number of trials that we want to perform. 
        */
        double sum = 0, price = 0;
        int n = rand_vals.size();

        for(int i=0; i<n; i++){
            price = MC_SDE_underlying_price(rand_vals[i]);
            sum += payout(price);
        }

        return exp(-r*T)*sum/n;
    }

    double variance(std::vector<double> data){
        double mean = 0, var = 0;

        for(double x : data){
            mean += x;
        }

        mean = mean/data.size();

        for(double x : data){
            var += (x-mean)*(x-mean);
        }

        var = var/data.size();

        return var;
    }

    double std_error(double std_dev, double n){
        return std_dev/sqrt(n);
    }

};

class Call : public Option{
    public:

    using Option::Option;

    double payout(double s){
        return std::max(s-K, 0.0);
    }

    double BS_price(double S_, double K_, double r_, double d_, double T_, double sigma_){     
        double D1 = d1(S_, K_, r_, d_, T_, sigma_);
        double D2 = d2(S_, K_, r_, d_, T_, sigma_);
        return S_*exp(-d_*T_)*N(D1)-K*exp(-r_*T_)*N(D2);
    }

    double BS_Delta(double S_, double K_, double r_, double d_, double T_, double sigma_){
        double D1 = d1(S_, K_, r_, d_, T_, sigma_);
        return exp(-d*T)*N(D1);
    }

    double BS_Gamma(double S_, double K_, double r_, double d_, double T_, double sigma_){
        double D1 = d1(S_, K_, r_, d_, T_, sigma_);
        return exp(-d*T)*N_prime(D1)/(S*sigma*sqrt(T));
    }

    // Also known as Vega, if you don't like the Greek alphabet.
    double BS_Kappa(double S_, double K_, double r_, double d_, double T_, double sigma_){
        double D1 = d1(S_, K_, r_, d_, T_, sigma_);
        return exp(-d*T)*S*sqrt(T)*N_prime(D1);
    }

    double BS_Rho(double S_, double K_, double r_, double d_, double T_, double sigma_){
        double D2 = d2(S_, K_, r_, d_, T_, sigma_);
        return K*exp(-d*T)*T*N(D2);
    }

    double BS_Theta(double S_, double K_, double r_, double d_, double T_, double sigma_){
        double D1 = d1(S_, K_, r_, d_, T_, sigma_);
        double D2 = d2(S_, K_, r_, d_, T_, sigma_);
        return S*exp(-d*T)*d*N(D1) - K*exp(-r*T)*r*N(D2) - S*exp(-r*T)*sigma*N_prime(D1)/(2*sqrt(T));
    }

    double FD_Delta(double S_, double K_, double r_, double d_, double T_, double sigma_, double epsilon){
        return (BS_price(S+epsilon, K, r, d, T, sigma)-BS_price(S, K, r, d, T, sigma))/epsilon;
    }

    double FD_Gamma(double epsilon){
        return (BS_price(S+epsilon, K, r, d, T, sigma) - 2*BS_price(S, K, r, d, T, sigma) + BS_price(S-epsilon, K, r, d, T, sigma))/(epsilon*epsilon);
    }

    double FD_Kappa(double epsilon){
        return (BS_price(S, K, r, d, T, sigma+epsilon)-BS_price(S, K, r, d, T, sigma))/epsilon;
    }

    double FD_Rho(double epsilon){
        return (BS_price(S, K, r+epsilon, d, T, sigma)-BS_price(S, K, r, d, T, sigma))/epsilon;
    }

    double FD_Theta(double epsilon){
        return -(BS_price(S, K, r, d, T+epsilon, sigma)-BS_price(S, K, r, d, T, sigma))/epsilon;
    }
};

int main(){

}