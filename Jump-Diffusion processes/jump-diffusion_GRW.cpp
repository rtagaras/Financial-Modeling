#include <iostream>
#include <random>
#include <string>
#include <utility>
#include <vector>
#include <fstream>
#include <algorithm>
#include <Eigen/Dense>

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

std::vector<double> arrival_times(double t, double lambda){
    /*
    Sample exponential distribution to give arrival times for Poisson jumps that occur before time t
    */

    double s = 0, sum = 0, dt, u;
    std::vector<double> times;

    while(s < t){
        u = gen_uniform(0,1);
        dt = -log(1-u)/lambda;
        s += dt;

        // it's possible that t_{n-1} < t but t_n > t, so we need to make sure that we don't add t_n to the vector if this is the case
        if(s < t){
            times.push_back(s);
        }
    }

    return times;
}

struct J_D_GRW{
    /*
    Simulate a jump-diffusion process with drift mu, volatility sigma, and intensity lambda
    T gives number of time units for which to calculate the path, dt is number of steps per unit time, s_0 is initial security price at day zero, dist_mu
    and dist_sigma are the parameters for the (normal) distribution from which jump sizes are sampled
    */

    double mu, sigma, lambda, T, dt, s_0, s, dist_mu, dist_sigma;
    int n;

    J_D_GRW(double s_0_, double mu_, double sigma_, double T_, double dt_, double lambda_, double dist_mu_, double dist_sigma_){
        s_0 = s_0_;
        mu = mu_;
        sigma = sigma_;
        lambda = lambda_;
        T = T_; 
        dt = dt_;
        dist_mu = dist_mu_;
        dist_sigma = dist_sigma_;
        n  = T/dt;
    }

    // calculate price as a function of time, measured in years
    // returns a vector of prices and a vector of time steps
    std::pair<std::vector<double>, std::vector<double>> path(){
        
        std::pair<std::vector<double>, std::vector<double>> ret;
        std::vector<double> data, time;
        std::vector<double> jump_times = arrival_times(T, lambda);
        int jump_count = 0;
        s = s_0;
        double dt1, j, z, next_jump_time, t = 0;

        // create a vector with all the time steps at which we want to calculate a value
        // Most of these will be of the form m*dt, where m is a positive integer. The steps at which a jump occurs will be located off of this grid.
        for(int k=0; k<n; k++){
            time.push_back(k*dt);
        }

        for(int k=0; k<jump_times.size(); k++){
            time.push_back(jump_times[k]);
        }

        std::sort(time.begin(), time.end());

        for(int i=0; i<time.size()-1; i++){

            if(jump_count < jump_times.size()){
                next_jump_time = jump_times[jump_count];
            }

            dt1 = time[i+1]-time[i];            

            // If a jump isn't coming next, do a standard geometric Brownian motion calculation.
            if(time[i+1] < next_jump_time || jump_count > jump_times.size()){
                z = gen_norm(0,1);
                s = s*(1.0 + mu*dt1+sigma*z*sqrt(dt1));
                data.push_back(s);
            }

            // Otherwise, calculate the jump.
            else{
                
                // First, interpolate to the jump point using GBM
                z = gen_norm(0,1);
                s = s*(1.0 + mu*dt1+sigma*z*sqrt(dt1));
                
                // Calculate the jump at the jump point
                j = gen_norm(dist_mu, dist_sigma*dist_sigma);
                s = s*j;
                data.push_back(s);
                jump_count++;
            }
        }

        // Delete the last time step, since we don't use it
        time.pop_back();

        ret = std::make_pair(data, time);
        return ret;
    }
};

int main(){

    J_D_GRW g = J_D_GRW(100, 0.03, 0.4, 100/365., 0.1/365., 0.1*365, 1, 0.06);
    std::pair<std::vector<double>, std::vector<double>> values;
    std::vector<double> p,t;
    std::vector<double> endpoints;

    for(int i=0; i<100000; i++){
        if(i % 100 == 0){
            std::cout << i << std::endl;
            output(p, "data_" + std::to_string(i));
            output(t, "times_" + std::to_string(i));
        }
        values = g.path();
        p = values.first;
        t = values.second;
      
        

        endpoints.push_back(p.back());
    }

    output(endpoints, "endpoints");

    return 0;
}