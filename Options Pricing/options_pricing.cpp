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

struct GRW{
    /*
    geometric random walk with drift mu and volatility sigma
    T_max gives number of time units for which to calculate the path, dt is number of steps per unit time
    s_0 is initial security price at day zero
    */

    double mu, sigma, T_max, dt, s_0, s, z;
    int n;

    GRW(double s_0_, double mu_, double sigma_, double T_max_, double dt_){
        s_0 = s_0_;
        mu = mu_;
        sigma = sigma_;
        T_max = T_max_; 
        dt = dt_;
        n  = T_max/dt;
    }

    // calculate price as a function of time, measured in days
    std::vector<double> path(){
        std::vector<double> data;
        s = s_0;
        
        for(int i=0; i<n; i++){

            z = gen_norm(0.0, 1.0);
            s = s*(1.0 + mu*dt+sigma*z*sqrt(dt));

            data.push_back(s);
        }

        return data;
    }
};

class Option{
    public:
    double strike, s_0, r, mu, sigma, dt;
    int expiry_time;
    std::string type;

    Option(std::string type_, double strike_, int expiry_time_, double dt_, double s_0_, double r_, double mu_, double sigma_){
        type = type_;
        strike = strike_;
        expiry_time = expiry_time_;
        s_0 = s_0_;
        r = r_;
        mu = mu_;
        sigma = sigma_;
        dt = dt_;
    }

    double payout(double x){

        if(type == "call"){
            return std::max(x - strike, 0.0);
        }

        else if(type == "put"){
            return std::max(strike - x, 0.0);
        }

        else{
            std::cout << "Invalid option type" << std::endl;
            return 0;
        }
    }
};

class European_Option : public Option{
    
    public:
    using Option::Option;

    // We calculate the expected payoff and discount the price by using the risk free rate, following the risk-neutral principal. This is computationally
    // simpler than implementing a recursive solution, although potentially less flexible. 
    double binomial_lattice_price(int num_trials){

        double E = 0.0, x = 0.0;
        
        // // start by using the u=1/d convention for binomial factors; if this leads to unphysical results, use p=1/2 convention instead
        // double A = (exp(-mu*dt)+exp((mu+sigma*sigma)*dt))/2.;
        // double d = A - sqrt(A*A - 1.0);
        // double u = A + sqrt(A*A - 1.0);
        // double p = (exp(mu*dt)-d)/(u-d);

        //if(p <= 0 || p >= 1){
            double d = exp(mu*dt)*(1.0 - sqrt(exp(sigma*sigma*dt) - 1.0));
            double u = exp(mu*dt)*(1.0 + sqrt(exp(sigma*sigma*dt) - 1.0));
            double p = 0.5;
        //}

        // calculate many paths through the lattice, adding the final payout to a running total
        for(int i=0; i<num_trials; i++){
            double s = s_0;

            for(int j=0; j<expiry_time; j++){
                x = gen_uniform(0,1);

                if(x > p){
                    s = s*u;
                }

                else{
                    s = s*d;
                }
            }

            E += payout(s);
        }

        // return the average payout, discounted by the RFR to t=0
        E = E/num_trials;
        return exp(-expiry_time*r*dt)*E;
    }

    // prices should follow geometric brownian motion - use a GRW to simulate end price
    double GRW_price(int n){
        double x = 0.0;
        std::vector<double> p;

        for(int i=0; i<n; i++){
            GRW g = GRW(s_0, mu, sigma, expiry_time/365., dt);

            p = g.path();
            x += payout(p.back());
        }

        return exp(-r*expiry_time/365.)*x/n;
    }
};

class American_Option : public Option{
    
    public:
    using Option::Option;

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
    
        Eigen::MatrixXd underlying_lattice(expiry_time+1, expiry_time+1);
        Eigen::MatrixXd option_lattice(expiry_time+1, expiry_time+1);

        // calculate underlying values
        for(int i=0; i<expiry_time+1; i++){
            for(int j=0; j<i+1; j++){
                underlying_lattice(i,j) = s_0*pow(u,j)*pow(d,i-j);
            }
        }

        // calculate all possible option payouts
        for(int j=0; j<expiry_time+1; j++){
            option_lattice(expiry_time, j) = payout(underlying_lattice(expiry_time, j));
        }

        // backpropagate option value to t=0
        for(int i=expiry_time-1; i>=0; i--){
            for(int j=0; j<i+1; j++){

                //calculate option price from child nodes and compare to value obtained by exercising the option
                x =  exp(-r*dt)*(p*option_lattice(i+1,j+1) + (1-p)*option_lattice(i+1,j));
                option_lattice(i,j) = std::max(x, payout(underlying_lattice(i,j)));
            }
        }

        return option_lattice(0,0);
    }

};

class Asian_Option : public Option{

};

int main(){

    European_Option O = European_Option("put", 49.0, 4.0, 1/365., 50.0, 0.26, 0.26, 0.4);

    std::cout << "Binomial lattice price: " << O.binomial_lattice_price(10000) << std::endl << "GRW price: " << O.GRW_price(10000) << std::endl;

    American_Option A = American_Option("put", 49.0, 4.0, 1/365., 50.0, 0.26, 0.26, 0.4);

    std::cout << "American price: " << A.binomial_lattice_price() << std::endl;
    return 0;
}