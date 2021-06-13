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