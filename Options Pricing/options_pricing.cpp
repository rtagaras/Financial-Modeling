#include <iostream>
#include <random>
#include <fstream>
#include <string>
#include <vector>
#include <algorithm>

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
    int option_type = 0;

    public:

    European_Option(std::string type, double strike_, double expiry_time_, double underlying_value_, double risk_free_rate_, double drift_, double variance_){
        if(type == "call"){
            option_type = 0;
        }

        else if(type == "put"){
            option_type = 1;
        }

        else
            std::cout << "Invalid option type" << std::endl;

        strike = strike_;
        expiry_time = expiry_time_;
        underlying_value = underlying_value_;
        risk_free_rate = risk_free_rate_;
        drift = drift_;
        variance = variance_;
    }

    double payout(){

        if(option_type == 0){
            return std::max(underlying_value - strike, 0.0);
        }

        else if(option_type == 1){
            return std::max(strike - underlying_value, 0.0);
        }

        else{
            std::cout << "Invalid option type" << std::endl;
            return 0;
        }
    }

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

        for(int i=0; i<num_trials; i++){
            underlying_value = s_0;

            for(int j=0; j<expiry_time; j++){
                x = gen_uniform(0,1);

                if(x < 1.0 - p){
                    underlying_value = underlying_value*u;
                }

                else{
                    underlying_value = underlying_value*d;
                }
            }

            E += payout();
        }

        E = E/num_trials;

        return exp(-expiry_time*risk_free_rate*dt)*E;
    }

};

int main(){

    European_Option O = European_Option("put", 49, 4, 50, 0.26, 0.26, 0.4);

    double val = O.binomial_lattice_price(50, 10000);

    std::cout << "Option price is: " << val << std::endl;

    return 0;
}