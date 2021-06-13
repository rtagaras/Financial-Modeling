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

    double mu, sigma, T_max, dt, s_0, end_value, s, z, RFR, yield;
    int n, days_between_divs;

    std::vector<double> data;

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
    }

    // calculate price as a function of time, measured in days
    void path(){
        s = s_0;
        
        for(int i=0; i<n; i++){

            z = gen_norm(0.0, 1.0);
            s = s*(1.0 + mu*dt+sigma*z*sqrt(dt));

            data.push_back(s);
        }

        end_value = data.back();
    }

    // calculate price with dividends, assuming that dividends are reinvested at risk-free rate
    void path_with_dividends(){

        s = s_0;
        double div = 0;
        double days_since_div = 0;

        for(int i=0; i<n; i++){
            z = gen_norm(0.0, 1.0);
            s = s*(1.0 + mu*dt+sigma*z*sqrt(dt));
            days_since_div++;

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
        end_value = s + div;

        //return data;
    }
};

int main(){

    //Create 1000 instances of geometric brownian motion. Here, we consider the paths and also store the ending values.
    std::string filename = "";
    std::vector<double> end_values;

    for(int i=0; i<100000; i++){
        GRW g = GRW(100.0, 0.0, 0.04, 60.0, 0.1);
        //std::vector<double> p = g.path();
        g.path();

        filename = "data"+ std::to_string(i);
        output(g.data, filename);

        end_values.push_back(g.end_value);
    }

    output(end_values, "end_values");
 
    return 0;
}