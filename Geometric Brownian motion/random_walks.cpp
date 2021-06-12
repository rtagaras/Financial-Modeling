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

// arithmetic random walk with drift mu and Weiner variance parameter sigma
std::vector<double> arw(double x, int T, double dt, double mu, double sigma){

    int n = T/dt;
    double z = 0;
    double dx = 0;

    std::vector<double> data;

    for(int i=0; i<n; i++){
        z = gen_norm(0,1);
        dx = mu*dt+sigma*z*sqrt(dt);        
        x += dx;

        data.push_back(x);
    }

    return data;
}

// geometric random walk with drift mu and volatility sigma
struct GRW{

    double mu = 0.0;
    double sigma = 0.0;
    double T_max = 0.0;
    double dt = 0.0;
    double s_0 = 0.0;
    double end_value = 0.0;

    GRW(double s_in, double m, double s, double Tm, double step){
        s_0 = s_in;
        mu = m;
        sigma = s;
        T_max = Tm;
        dt = step;
    }

    std::vector<double> path(){

        int n = T_max/dt;
        std::vector<double> data;

        double z = 0;
        double s = s_0;
        
        for(int i=0; i<n; i++){

            z = gen_norm(0.0, 1.0);
            s = s*(1.0 + mu*dt+sigma*z*sqrt(dt));

            data.push_back(s);
        }

        end_value = data.back();
        return data;
    }
};

int main(){

    //Create 1000 instances of geometric brownian motion. Here, we consider the paths and also store the ending values.
    std::string filename = "";
    std::vector<double> end_values;

    for(int i=0; i<100000; i++){
        GRW g = GRW(100.0, 0.0, 0.04, 60.0, 0.1);
        std::vector<double> p = g.path();

        filename = "data"+ std::to_string(i);
        output(p, filename);

        end_values.push_back(p.back());
    }

    output(end_values, "end_values");
 
    return 0;
}