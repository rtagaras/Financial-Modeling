#include <iostream>
#include <vector>
#include <random>
#include <fstream>
#include <sstream>

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


std::vector<std::pair<std::string, std::vector<double>>> read_csv(std::string filename){
    // Reads a CSV file into a vector of <string, vector<double>> pairs where
    // each pair represents <column name, column values>

    // Create a vector of <string, int vector> pairs to store the result
    std::vector<std::pair<std::string, std::vector<double>>> result;

    // Create an input filestream
    std::ifstream myFile(filename);

    // Make sure the file is open
    if(!myFile.is_open()) throw std::runtime_error("Could not open file");

    // Helper vars
    std::string line, colname;
    double val;

    // Read the column names
    if(myFile.good()){
        // Extract the first line in the file
        std::getline(myFile, line);

        // Create a stringstream from line
        std::stringstream ss(line);

        // Extract each column name
        while(std::getline(ss, colname, ',')){
            
            // Initialize and add <colname, int vector> pairs to result
            result.push_back({colname, std::vector<double> {}});
        }
    }

    // Read data, line by line
    while(std::getline(myFile, line)){
        // Create a stringstream of the current line
        std::stringstream ss(line);
        
        // Keep track of the current column index
        int colIdx = 0;
        
        // Extract each integer
        while(ss >> val){
            
            // Add the current integer to the 'colIdx' column's values vector
            result.at(colIdx).second.push_back(val);
            
            // If the next token is a comma, ignore it and move on
            if(ss.peek() == ',') ss.ignore();
            
            // Increment the column index
            colIdx++;
        }
    }

    // Close file
    myFile.close();

    return result;
}

class RNG{

    private:
    std::mt19937 gen;

    public:

    RNG() : gen((std::random_device())()) {}
    
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

    std::vector<double> create_norm_vector(int n, double mean, double variance){
        /*
        Creates a vector of N(mean, variance) random variables
        */
        std::vector<double> data;

        for(int i=0; i<n; i++){
            data.push_back(gen_norm(mean, variance));
        }

        return data;
    }
};

class Ho_Lee_Model : RNG{
    /*
    Model bond prices using the Ho-Lee model. sigma gives the volatility, f and f_last are values taken from the instantaneous forward rate curve.
    T_max gives number of time units for which to calculate the path, dt is number of steps per unit time
    r_0 is initial interest rate at day zero
    */

    public:
    double f, f_last, sigma, T_max, dt, r_0, z, r;
    int n;
    std::vector<std::pair<std::string, std::vector<double>>> data;

    Ho_Lee_Model(double r_0_, double sigma_, double T_max_, double dt_) : r_0(r_0_), sigma(sigma_), T_max(T_max_), dt(dt_){
        n  = T_max/dt;
        data = read_csv("FED-SVENF-2.csv");
        f_last = data[0].second[0];
        //data = read_csv("test.csv");
    }

    // get the mth value from the 1-year instantaenous forward rate
    double get_rate_curve_value(int m){    
        return data[0].second[m];
    }

    // calculate price as a function of time, measured in years
    std::vector<double> path(){
        std::vector<double> vals;
        r = r_0;
        
        for(int i=1; i<n; i++){

            f = data[0].second[i];

            z = gen_norm(0,1);
            r += (f-f_last) + 0.5*sigma*sigma*(2*i-1)*dt*dt + sigma*z*sqrt(dt);
            f_last = f;
            vals.push_back(r);
        }

        return vals;
    }
};

int main(){

    Ho_Lee_Model H = Ho_Lee_Model(1.0,0.2,60/365.,1/365.);
    std::vector<double> v = H.path();
    output(v, "HL_data");

    // Vasicek_model V = Vasicek_model(1.0,0.2,60/365.,1/365.);
    // v = V.path();
    // output(v, "V_data");
}