#include <iostream>
#include <random>
#include <fstream>
#include <string>
#include <vector>
#include <algorithm>
#include <Eigen/Dense>

using Eigen::MatrixXd;

std::random_device rd;
std::mt19937 gen(rd());
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

class European_Option{
    private:

    double strike = 0.0, expiry_time = 0.0, underlying_value = 0.0, risk_free_rate = 0.0, drift = 0.0, variance = 0.0;
    std::string type;

    public:

    European_Option(std::string type_, double strike_, double expiry_time_, double underlying_value_, double risk_free_rate_, double drift_, double variance_){
        type = type_;
        strike = strike_;
        expiry_time = expiry_time_;
        underlying_value = underlying_value_;
        risk_free_rate = risk_free_rate_;
        drift = drift_;
        variance = variance_;
    }

    double payout(){

        if(type == "call"){
            return std::max(underlying_value - strike, 0.0);
        }

        else if(type == "put"){
            return std::max(strike - underlying_value, 0.0);
        }

        else{
            std::cout << "Invalid option type" << std::endl;
            return 0;
        }
    }

    // We calculate the expected payoff and discount the price by using the risk free rate, following the risk-neutral principal. This is computationally
    // simpler than implementing a recursive solution, although potentially less flexible. 
    double binomial_lattice_price(double s_0, int num_trials){

        double E = 0.0, x = 0.0;
        double dt = 1./365;
        
        // start by using the u=1/d convention for binomial factors; if this leads to unphysical results, use p=1/2 convention instead
        double A = (exp(-drift*dt)+exp((drift+variance*variance)*dt))/2.;
        double d = A - sqrt(A*A - 1.0);
        double u = A + sqrt(A*A - 1.0);
        double p = (exp(drift*dt)-d)/(u-d);

        if(p <= 0 || p >= 1){
            d = exp(drift*dt)*(1.0 - sqrt(exp(variance*variance*dt) - 1.0));
            u = exp(drift*dt)*(1.0 + sqrt(exp(variance*variance*dt) - 1.0));
            p = 0.5;
        }

        // calculate many paths through the lattice, adding the final payout to a running total
        for(int i=0; i<num_trials; i++){
            underlying_value = s_0;

            for(int j=0; j<expiry_time; j++){
                x = gen_uniform(0,1);

                if(x > p){
                    underlying_value = underlying_value*u;
                }

                else{
                    underlying_value = underlying_value*d;
                }
            }

            E += payout();
        }

        // return the average payout, discounted by the RFR to t=0
        E = E/num_trials;
        return exp(-expiry_time*risk_free_rate*dt)*E;
    }

};

class American_Option{
    private:

    double strike, s_0, r, mu, sigma;
    std::string type;
    int expiry_time;

    public:

    American_Option(std::string type_, double strike_, int expiry_time_, double s_0_, double r_, double mu_, double sigma_){
        type = type_;
        strike = strike_;
        expiry_time = expiry_time_;
        s_0 = s_0_;
        r = r_;
        mu = mu_;
        sigma = sigma_;
    }

    double payout(double s){

        if(type == "call"){
            return std::max(s - strike, 0.0);
        }

        else if(type == "put"){
            return std::max(strike - s, 0.0);
        }

        else{
            std::cout << "Invalid option type" << std::endl;
            return 0;
        }
    }

    double binomial_lattice_price(){

        double x = 0.0;
        double dt = 1./365;
        
        // start by using the u=1/d convention for binomial factors; if this leads to unphysical results, use p=1/2 convention instead
        double A = (exp(-mu*dt)+exp((mu+sigma*sigma)*dt))/2.;
        double d = A - sqrt(A*A - 1.0);
        double u = A + sqrt(A*A - 1.0);
        double p = (exp(mu*dt)-d)/(u-d);

        //if(p <= 0 || p >= 1){
            d = exp(mu*dt)*(1.0 - sqrt(exp(sigma*sigma*dt) - 1.0));
            u = exp(mu*dt)*(1.0 + sqrt(exp(sigma*sigma*dt) - 1.0));
            p = 0.5;
        //}
    
        Eigen::MatrixXd underlying_lattice(expiry_time, expiry_time);
        Eigen::MatrixXd option_lattice(expiry_time, expiry_time);

        // calculate underlying values
        for(int i=0; i<expiry_time; i++){
            for(int j=0; j<i+1; j++){
                underlying_lattice(i,j) = s_0*pow(u,j)*pow(d,i-j);
            }
        }

        // calculate all possible option payouts
        for(int j=0; j<expiry_time; j++){
            option_lattice(expiry_time-1, j) = payout(underlying_lattice(expiry_time-1, j));
        }

        // backpropagate option value to t=0
        for(int i=expiry_time-2; i>=0; i--){
            for(int j=0; j<i+1; j++){

                //calculate option price from child nodes and compare to value obtained by exercising the option
                x =  exp(-r*dt)*(p*option_lattice(i+1,j+1) + (1-p)*option_lattice(i+1,j));
                option_lattice(i,j) = std::max(x, payout(underlying_lattice(i,j)));
            }
        }

        return option_lattice(0,0);
    }

};

int main(){
  
    American_Option O = American_Option("put", 49, 5, 50, 0.26, 0.26, 0.4);

    std::cout << O.binomial_lattice_price() << std::endl;

    return 0;
}