#include <iostream>
#include <random>
#include <fstream>
#include <string>
#include <vector>
#include <sstream>

//output a vector
template <class T>
void output(std::vector<T> v, std::string filename){
	std::stringstream ss;
	std::string s;
		
	std::string name = "./Data/" + filename + ".txt";
	std::ofstream ofs(name);
	for(auto p =v.begin(); p!=v.end(); ++p){
        ss << *p;
        ss >> s;
        ofs << s;
        ss.clear();
			
		ofs << std::endl;
	}
	ofs.close();
}

std::random_device rd;
std::mt19937 gen(rd());
double gen_norm(double mean, double variance){
    std::normal_distribution<double> d(mean, sqrt(variance)); 
    return d(gen);
}

// geometric random walk with drift mu and volatility sigma
// T_max gives number of days to calculate path for
// s_0 is initial stock price at day zero
struct GRW{

    double mu, sigma, T_max, dt, s_0, s, z, RFR, yield, div = 0.0;
    int n, days_between_divs, days_since_div, steps_per_day = 0;

    GRW(double s_0_, double mu_, double sigma_, double T_max_, double dt_, int days_between_divs_ = 0, double RFR_ = 0.0, double yield_ = 0.0){
        s_0 = s_0_;
        mu = mu_;
        sigma = sigma_;
        T_max = T_max_; 
        dt = dt_;
        n  = T_max/dt;
        days_between_divs = days_between_divs_;
        RFR = RFR_;
        yield = yield_;
        steps_per_day = 1/dt;
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

    // calculate price with dividends, assuming that dividends are reinvested at risk-free rate
    std::vector<double> path_with_dividends(){
        std::vector<double> data;

        s = s_0;
        div = 0;
        days_since_div = 0;

        for(int i=0; i<n; i++){
            z = gen_norm(0.0, 1.0);
            s = s*(1.0 + mu*dt+sigma*z*sqrt(dt));
            
            // This should count the number of time steps in a single day and only increment days_since_div if a complete day has passed. 
            if(i % steps_per_day == 0){
                days_since_div++;
            }

            if(days_since_div == days_between_divs){
                
                // grow accumulated dividend at risk-free rate and add dividend for current period
                div += div*RFR + s*yield;
                
                // reduce current price by dividend amount
                s -= s*yield;

                days_since_div = 0;
            }

            data.push_back(s);
        }

        // After reaching the end of the calculation period, we need to make sure that the total dividend growth after the last dividend has been 
        // accounted for.
        div += div*RFR*days_since_div/days_between_divs;
        data.push_back(s+div);

        return data;
    }
};

int main(){

    std::string filename = "";
    std::vector<double> end_values;
    std::vector<double> p;
    double end_val = 0.0;
    int profit_count = 0;

    // Parameters
    int num_trials = 10000;
    double starting_price = 100.0;
    double mu = 0.1;
    double sigma = 0.04;
    double T_max = 365.0;
    double dt = 0.1;
    int days_between_divs = 91;
    double RFR = 0.02;
    double yield = 0.08;

    // Create 100000 instances of geometric brownian motion with a dividend. 
    // Here, we consider the paths and also store the ending values.
    GRW g = GRW(starting_price, mu, sigma, T_max, dt, days_between_divs, RFR, yield);

    for(int i=0; i<num_trials; i++){
        p = g.path_with_dividends();
        end_val = p.back();

        if (end_val > 100){
            profit_count++;
        }

        filename = "data"+ std::to_string(i);
        output(p, filename);

        end_values.push_back(end_val);
    }

    output(end_values, "end_values");
    std::cout << "Probability of profit: " << (double) profit_count/num_trials << std::endl;
 
    return 0;
}