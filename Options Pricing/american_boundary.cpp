#include <cmath>
#include <iostream>
#include <random>
#include <string>
#include <vector>
#include <fstream>
#include <Eigen/Dense>
#include <math.h>
#include <cfloat>

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

double boundary(double t, double a1, double a2, double b1, double b2, double c1, double c2){
    /*
    A particular choice of exercise boundary parameterization. t gives the time to expiry. The a,b,c constants will be determined by some optimization
    method. 
    */
    return a1*std::log(b1*std::pow(t,c1) + 1.) + a2*std::log(b2*std::pow(t,c2) + 1.);
}

std::vector<double> generate_vector(std::vector<double> x, double epsilon){
    /*
    Creates a random state "close" to x. This is done by sampling the 6-ball with radius epsilon centered at x.
    The boundary parameters are supposed to be positive, so I took the absolute value. I'm not sure if this is better or worse than returning zero in the
    case of a negative value. 
    */

    double r = gen_uniform(0, epsilon) + sqrt(x[0]*x[0] + x[1]*x[1] + x[2]*x[2] + x[3]*x[3] + x[4]*x[4] + x[5]*x[5]);
    double theta_1 = gen_uniform(0, M_PI);
    double theta_2 = gen_uniform(0, M_PI);
    double theta_3 = gen_uniform(0, M_PI);
    double theta_4 = gen_uniform(0, M_PI);
    double theta_5 = gen_uniform(0, 2.*M_PI);

    // convert from spherical coordinates to cartesian
    double x1 = abs(r*cos(theta_1));
    double x2 = abs(r*sin(theta_1)*cos(theta_2));
    double x3 = abs(r*sin(theta_1)*sin(theta_2)*cos(theta_3));
    double x4 = abs(r*sin(theta_1)*sin(theta_2)*sin(theta_3)*cos(theta_4));
    double x5 = abs(r*sin(theta_1)*sin(theta_2)*sin(theta_3)*sin(theta_4)*cos(theta_5));
    double x6 = abs(r*sin(theta_1)*sin(theta_2)*sin(theta_3)*sin(theta_4)*sin(theta_5));

    std::vector<double> v = {x1,x2,x3,x4,x5,x6};
    return v;
}

class American_Option{
    
    private:
    double strike, s_0, r, mu, sigma, dt, expiry_time;
    int n;
    std::string type;
    
    public:
    std::vector<double> samples;

    American_Option(std::string type_, double strike_, double expiry_time_, double dt_, double s_0_, double r_, double mu_, double sigma_){
        type = type_;
        strike = strike_;
        expiry_time = expiry_time_;
        s_0 = s_0_;
        dt = dt_;
        r = r_;
        mu = mu_;
        sigma = sigma_;
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
    std::pair<std::vector<double>, double> optimize_exercise_boundary(double epsilon, double T_max, double dt){
        
        std::vector<double> x = {gen_uniform(0,10),gen_uniform(0,10),gen_uniform(0,10),gen_uniform(0,10),gen_uniform(0,10),gen_uniform(0,10)};
        
        std::vector<double> y;
        double E = 1./GRW_price(x,100);
        double E_trial = 0;
        double delta_E = 0;
        double z = 0;
        double T = T_max;
        int m = T_max/dt;

        for(int i=0; i<m; i++){

            // create a random 6d vector with length less than epsilon. The vector components give the parameters for the current trial.
            y = generate_vector(x, epsilon);
            E_trial = 1./GRW_price(y,100);
            delta_E = E_trial - E;

            // with probability e^{-Delta_E/T}, update the choice of vector
            z = gen_uniform(0,1);
            if(z <= exp(-delta_E/T)){
                x = y;
                E = E_trial;
            }

            T -= dt;
        }

        return std::make_pair(x, 1/E);
    }

    double GRW_price(std::vector<double> parameters, int m){
        double sum = 0;
        double a1 = parameters[0], a2 = parameters[1], b1 = parameters[2], b2 = parameters[3], c1 = parameters[4], c2 = parameters[5];
        double z;
        
        for(int i=0; i<m; i++){
            
            double s = s_0;
            for(int i=0; i<n; i++){
                z = gen_norm(0,1);
                s = s*(1. + mu*dt + sigma*z*sqrt(dt));

                // if current price (as a percentage) is in the money, add payoff discounted by risk free rate using current time
                if(payout(s)/strike >= boundary(expiry_time - n*dt, a1, a2, b1, b2, c1, c2)){
                    sum += exp(-r*n*dt)*payout(s);
                }
            }
            sum += payout(s);
        }
        return exp(-r*expiry_time)*sum/m;
    }
    
};

int main(){

    std::string type = "put";
    double strike = 100;
    double expiry_time = 90/365.;
    double dt = 1/365.;
    double s_0 = 100;
    double r = 0.1;
    double mu = 0.1;
    double sigma = 0.3;
    double epsilon = 0.1;
    double T_max = 50;

    American_Option O = American_Option(type, strike, expiry_time, dt, s_0, r, mu, sigma);
    std::pair<std::vector<double>, double> data = O.optimize_exercise_boundary(epsilon, T_max, 0.01);
    std::cout << data.second << std::endl;
}