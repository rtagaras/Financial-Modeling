#include <iostream>
#include <random>
#include <string>
#include <vector>
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
    /*
    Base class for a general option

    s_0 is initial underlying value, r is risk-free rate (measured in percent/year), mu is drift, sigma is volatility, expiry time is measured in days
    */
    public:
    double strike, s_0, r, mu, sigma, dt;
    int expiry_time;
    std::string type;

    Option(std::string type_, double strike_, int expiry_time_, double dt_, double s_0_, double r_, double mu_, double sigma_){
        type = type_;
        strike = strike_;
        expiry_time = expiry_time_;
        s_0 = s_0_;
        dt = dt_/365.;
        
        // 
        r = r_;
        mu = mu_;
        sigma = sigma_;
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
            GRW g = GRW(s_0, mu, sigma, expiry_time, dt);

            p = g.path();
            x += payout(p.back());
        }

        return exp(-r*expiry_time)*x/n;
    }
};

class American_Option : public Option{
    
    public:
    using Option::Option;

    double binomial_lattice_price(){

        double x = 0.0;
        int n = expiry_time;
        
        // start by using the u=1/d convention for binomial factors; if this leads to unphysical results, use p=1/2 convention instead
        double A = (exp(-mu*dt)+exp((mu+sigma*sigma)*dt))/2.;
        double d = A - sqrt(A*A - 1.0);
        double u = A + sqrt(A*A - 1.0);
        double p = (exp(mu*dt)-d)/(u-d);

        if(p <= 0 || p >= 1){
            d = exp(mu*dt)*(1.0 - sqrt(exp(sigma*sigma*dt) - 1.0));
            u = exp(mu*dt)*(1.0 + sqrt(exp(sigma*sigma*dt) - 1.0));
            p = 0.5;
        }
    
        Eigen::MatrixXd underlying_lattice(n+1, n+1);
        Eigen::MatrixXd option_lattice(n+1, n+1);

        // calculate underlying values
        for(int i=0; i<n+1; i++){
            for(int j=0; j<i+1; j++){
                underlying_lattice(i,j) = s_0*pow(u,j)*pow(d,i-j);
            }
        }

        // calculate all possible option payouts
        for(int j=0; j<n+1; j++){
            option_lattice(n, j) = payout(underlying_lattice(n, j));
        }

        // backpropagate option value to t=0
        for(int i=n-1; i>=0; i--){
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
    public:
    using Option::Option;

    double fixed_strike_GRW_price(int n){
        double x = 0.0, s = 0.0;
        std::vector<double> p;

        for(int i=0; i<n; i++){
            GRW g = GRW(s_0, mu, sigma, expiry_time, dt);

            s = 0.0;
            p = g.path();

            // get average price over entire random walk
            for(int j=0; j<p.size(); j++){
                s += p[j]/p.size();
            }

            x += payout(s);
        }

        return exp(-r*expiry_time)*x/n;
    }
};

class Barrier_Option : public Option{

    private:
    double E = 0.0, s = 0.0, b, d1 = 0.0, d2 = 0.0, z, u;
    int n;
    bool crossed = 0;

    public:
    Barrier_Option(std::string type_, double strike_, double barrier_, int expiry_time_, double dt_, double s_0_, double r_, double mu_, double sigma_) 
        : Option(type_, strike_, expiry_time_, dt_, s_0_, r_, mu_, sigma_){
        
        b = barrier_;
        n = expiry_time/dt;
    }

    double GRW_price(int num_trials){
        for(int i=0; i<num_trials; i++){

            crossed = 0;
            s = s_0;
            d1 = b-s_0;

            for(int j=0; j<n; j++){

                // normal GRW calculation
                z = gen_norm(0,1);
                s = s*(1.0 + mu*dt + sigma*z*sqrt(dt));
                d2 = b-s;

                // use Brownian bridge to decide if barier has been crossed, which makes option value zero
                u = gen_uniform(0,1);
                if(u < exp(-2.*d1*d2/(sigma*sigma*dt))){
                    crossed = 1;
                }

                d1 = d2;
            }

            if(crossed == 0){
                E += payout(s);
            }
        }

        return E*exp(-r*expiry_time)/num_trials;
    }
};

int main(){

    std::string type = "put";
    double s_0 = 100;
    double strike = 100;
    double r = 0.03;
    double sigma = 0.2;
    int days_to_expiry = 60;
    double dt = 0.01;
    double barrier = 105;

    European_Option E1 = European_Option(type, strike, days_to_expiry, dt, s_0, r, r, sigma);
    std::cout << "European GRW price: " << E1.GRW_price(10000) << std::endl;

    European_Option E2 = European_Option(type, strike-51., days_to_expiry-56., dt+0.99, s_0-50., r+0.23, r+0.23, sigma+0.2);
    std::cout << "European lattice price: " << E2.binomial_lattice_price(10000) << std::endl;

    American_Option A = American_Option(type, strike-51., days_to_expiry-56., dt+0.99, s_0-50., r+0.23, r+0.23, sigma+0.2);
    std::cout << "American lattice price: " << A.binomial_lattice_price() << std::endl;

    Asian_Option As = Asian_Option(type, strike, days_to_expiry, dt*10, s_0, r, r, sigma);
    std::cout << "Asian GRW price: " << As.fixed_strike_GRW_price(100) << std::endl;

    Barrier_Option B = Barrier_Option(type, strike, barrier, days_to_expiry, dt, s_0, r, r, sigma);    
    std::cout << "Knock-out option GRW price: " << B.GRW_price(100) << std::endl;

    return 0;
}