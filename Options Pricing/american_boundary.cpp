#include <cmath>
#include <iostream>
#include <random>
#include <string>
#include <vector>
#include <fstream>
#include <Eigen/Dense>

std::random_device rd;
std::mt19937 gen(rd());
double gen_norm(double mean, double variance){
    /*
    Sample a normal distribution with given mean and variance
    */
    std::normal_distribution<double> d(mean, sqrt(variance)); 
    return d(gen);
}

double boundary(double t, double a1, double a2, double b1, double b2, double c1, double c2){
    /*
    A particular choice of exercise boundary parameterization. t gives the time to expiry. The a,b,c constants will be determined by some optimization
    method. 
    */
    return a1*std::log(b1*std::pow(t,c1) + 1.) + a2*std::log(b2*std::pow(t,c2) + 1.);
}

class American_Option{
    
    private:
    double strike, s_0, r, mu, sigma, dt, expiry_time, lambda, dist_mu, dist_sigma, a1, a2, b1, b2, c1, c2;
    int n;
    std::string type;
    
    public:
    std::vector<double> samples;

    American_Option(std::string type_, double strike_, double expiry_time_, double dt_, double s_0_, double r_, double mu_, double sigma_, double lambda_ = 0, double dist_mu_ = 0, double dist_sigma_ = 0){
        type = type_;
        strike = strike_;
        expiry_time = expiry_time_;
        s_0 = s_0_;
        dt = dt_;
        r = r_;
        mu = mu_;
        sigma = sigma_;
        lambda = lambda_;
        dist_mu = dist_mu_;
        dist_sigma = dist_sigma_;
        n = expiry_time/dt;
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

    // Finds the optimal values of {a1, a2, b1, b2, c1, c2} so that the option price is maximized
    std::pair<std::vector<double>, double> optimize_exercise_boundary(double a, double b, std::vector<double> x, double epsilon, int m, double T_max){
        
        std::vector<double> y;
        double E = 1./GRW_price(x, 1);
        double E_trial = 0;
        double delta_E = 0;
        double z = 0;
        double T = T_max;

        for(int i=0; i<m; i++){
            // create a random 6d vector with length less than epsilon. The vector components give the parameters for the current trial.
            y = generate_vector(epsilon);
            E_trial = 1./GRW_price(y,1);
            delta_E = E_trial - E;

            // with probability e^{-Delta_E/T}, update the choice of vector
            z = gen_uniform(0,1);
            if(z <= exp(-delta_E/T)){
                x = y;
                E = E_trial;
            }

            T = a*std::pow(b,i);
        }

        return std::make_pair(x,1./GRW_price(x,1));
    }

    double GRW_price(int m){
        double sum = 0;

        for(int i=0; i<m; i++){
            double z;
            double s = s_0;
            for(int i=0; i<n; i++){
                z = gen_norm(0,1);
                s = s*(1. + mu*dt + sigma*z*sqrt(dt));

                // if current price (as a percentage) is in the money, add payoff discounted by risk free rate using current time
                if(payout(s)/strike >= boundary(expiry_time - n*dt, a1, a2, b1, b2, c1, c2)){
                    sum += exp(-r*n*dt)*payout(s);
                }
            }
            sum += exp(-r*expiry_time)*payout(s);
        }
        return sum/m;
    }
    
};

int main(){


}