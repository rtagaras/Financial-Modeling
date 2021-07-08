#include <iostream>
#include <random>
#include <fstream>
#include <string>
#include <vector>
#include <Eigen/Dense>
#include <optional>
 
using Eigen::MatrixXd;

//write a vector to a .txt file and store it in the Data subdirectory
template<class T>
void output(std::vector<T> v, std::string filename){
    std::ofstream ofs;
    ofs.open("./Data/" + filename + ".txt");
    
    for(int i=0; i<v.size(); i++){
        ofs << v[i] << '\n';
    }
    
    ofs.close();
}

// Sample a normal distribution with given mean and variance
std::random_device rd;
std::mt19937 gen(rd());
double gen_norm(double mean, double variance){
    std::normal_distribution<double> d(mean, sqrt(variance)); 
    return d(gen);
}

// geometric random walk with drift mu and volatility sigma
// T_max gives number of days for which to calculate the path
// s_0 is initial stock price at day zero
struct GRW{

    double mu = 0.0, sigma = 0.0, T_max = 0.0, dt = 0.0, s_0 = 0.0, s = 0.0, z = 0.0;
    int n = 0, steps_per_day = 0;

    GRW(double s_0_, double mu_, double sigma_, double T_max_, double dt_){
        s_0 = s_0_;
        mu = mu_;
        sigma = sigma_;
        T_max = T_max_; 
        dt = dt_;
        n  = T_max/dt;
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
};

// Given the correlations between pairs of variables and the individual standard deviations of each random variable, return the correlation matrix. 
// "correlations" should be a strictly upper triangular nxn matrix where the (i,j)th entry is the correlation coefficient for variables i and j.
// The ith element of "std_deviations" should be the standard deviation for the ith random variable.
//
// A Hermitian, positive definite matrix M can be decomposed as M = LL^T, where L is lower triangular.
// This class also includes the matrix L, which is used in generating correlated samples.
struct Correlation_Matrix{
    Eigen::MatrixXd C;
    Eigen::MatrixXd L;

    Correlation_Matrix(Eigen::VectorXd std_deviations, Eigen::MatrixXd correlations){
        int x = correlations.cols();
        int y = correlations.rows();

        Eigen::MatrixXd C_(x,y);

        for(int i=0; i<x; i++){
            C_(i,i) = std::pow(std_deviations(i), 2);
        }

        for(int i=0; i<y; i++){
            for(int j=i+1; j<x; j++){
                C_(i,j) = correlations(i,j) * std_deviations(i) * std_deviations(j);
                C_(j,i) = C_(i,j);
            }
        }

        C = C_;

        // Perform the Cholesky decomposition of correlation matrix
        // I would have included this as a member function instead of being part of the constructor, but checking that the correlation matrix is valid
        // already requires us to perform the decomposition. 
        Eigen::LLT<MatrixXd> lltOfA(C);
        L = lltOfA.matrixL();

        // Check to make sure that we have a valid correlation matrix
        if(lltOfA.info() == Eigen::NumericalIssue){
            std::cout << "Correlation matrix is likely not positive semidefinite." << std::endl;
        }   
    }
};

// Given a the matrix L obtained from the Cholesky decomposition of a correlation matrix C, generate a vector of correlated random variables. 
// L is generated by the "Correlation_matrix".
Eigen::VectorXd correlated_samples(Eigen::MatrixXd L){

    // Vector of uncorrelated N(0,1) random variables
    Eigen::VectorXd Z(L.cols());

    // Vector of random variables with correlations and standard deviations defined by choice of C. 
    Eigen::VectorXd X;

    for(int i=0; i<L.cols(); i++){
        Z(i) = gen_norm(0,1);
    }

    X = L*Z;

    return X;
}

// We may want to calculate prices of securities that are correlated with some predetermined market conditions.
// We consider random fluctuations around a piecewise function constructed by linearly interpolating between a set of specified (time, price) points. 
//
// mu and sigma give drift and volatility, respectively.
// fixed_times holds the days at which we would like to specify market prices. fixed_prices holds the prices we specify. 
// It is important to ensure that units of fixed_times are consistent with the units of time used by the random walk. 
//
// adjusted_prices contains the values of the GRW at each of the t_points, because we want to interpolate between y_i-g(t), not y_i alone. Otherwise, we
// wouldn't actually have the desired market conditions.
struct market_scenario{

    double mu = 0.0, sigma = 0.0, T_max = 0.0, dt = 0.0, steps_per_day = 0.0;
    std::vector<double> fixed_times, fixed_prices, adjusted_prices, market_prices;

    // fluctuations around the fixed prices
    GRW g;
    std::vector<double> p;

    market_scenario(double mu_, double sigma_, double T_max_, double dt_, std::vector<double> fixed_times_, std::vector<double> fixed_prices_) : g(fixed_prices_[0], mu_, sigma_, T_max_, dt_){
        mu = mu_;
        sigma = sigma_;
        T_max = T_max_;
        dt = dt_;
        steps_per_day = 1.0/dt;
        fixed_times = fixed_times_;
        fixed_prices = fixed_prices_;

        p = g.path();

        // Shift fixed prices by the value of the GRW at the corresponding time step. We will interpolate between shifted points and later add the GRW
        // back, so the final result is a GRW with specified values at chosen points. 
        for(int i=0; i<fixed_times.size(); i++){
            adjusted_prices.push_back(fixed_prices[i] - p[fixed_times[i]]);
        }
    }

    // This function takes values of the market at specific times and interpolates (linearly) between them, giving a base market value at times other 
    // than those specified. 
    std::optional<double> piecewise_linear_interpolation(double t){

        double a = fixed_times[0];
        double A = adjusted_prices[0];
        double b = 0.0, B = 0.0;

        for(int i=1; i<fixed_times.size(); i++){
            b = fixed_times[i];
            B = adjusted_prices[i];

            if(a <= t && t <= b){
                return (B*(t-a)-A*(t-b))/(b-a);
            }

            a = b;
            A = B;
        }

        return std::nullopt;
    }

    // returns a vector that contains the prices for the market, subject to the constraints given in the market_scenario constructor
    std::vector<double> scenario(){
        std::vector<double> mp;

        for(int i=0; i< T_max/dt; i++){
            mp.push_back(p[i] + piecewise_linear_interpolation(i*dt).value());
        }

        market_prices = mp;
        return mp;
    }

    // given a scenario, return the price increments 
    std::vector<double> price_increments(){

        std::vector<double> m = scenario();
        std::vector<double> PI;

        for(int i=1; i<m.size(); i++){
            PI.push_back((m[i]-m[i-1])/(m[i-1]*sigma*sqrt(dt)));
        }

        return PI;
    }
};

// A vector-valued geometric random walk where components have correlated price increments. 
//
// mu is a vector that holds the drift parameter for each component. Similarly, sigma is the volatility vector. 
// T_max gives the number of days for which to calculate the path.
// s_0 is the initial price vector at day zero.
// correlations is a strictly upper triangular matrix containing the correlation coefficients rho_{ij}
struct correlated_GRW{

    // final correlation matrix 
    Eigen::MatrixXd Corr;

    // L matrix given by Choelsky decomposition
    Eigen::MatrixXd l;

    Eigen::VectorXd s_0, s, mu, sigma, z;
    double T_max = 0.0, dt = 0.0;
    int n = 0, steps_per_day = 0;

    correlated_GRW(Eigen::VectorXd s_0_, Eigen::VectorXd mu_, Eigen::VectorXd sigma_, Eigen::MatrixXd correlations_, double T_max_, double dt_){
        s_0 = s_0_;
        mu = mu_;
        sigma = sigma_;
        T_max = T_max_; 
        dt = dt_;
        n  = T_max/dt;
        steps_per_day = 1/dt;

        Correlation_Matrix CM = Correlation_Matrix(sigma_, correlations_);
        Corr = CM.C;
        l = CM.L;
    }

    // calculate price as a function of time, measured in days
    std::vector<Eigen::VectorXd> path(){
        
        std::vector<Eigen::VectorXd> data;
        s = s_0;
        
        for(int i=0; i<n; i++){

            z = correlated_samples(l);
            s = s*(Eigen::VectorXd::Ones(s_0.size()) + mu*dt+sigma*z*sqrt(dt));
            data.push_back(s);
        }

        return data;
    }
};

// price of a stock correlated with a given market scenario
//
// s gives initial stock price, mu_s gives stock drift, sigma_s gives stock variance
struct market_GRW{

    double s = 0.0, mu = 0.0, rho = 0.0, sigma = 0.0, T_max = 0.0, dt = 0.0, n = 0.0;
    std::vector<double> increments, market_prices;

    market_GRW(double s_, double rho_, double mu_, double sigma_, double T_max_, double dt_, market_scenario m){
        
        s = s_;
        mu = mu_;
        rho = rho_;
        sigma = sigma_;
        T_max = T_max_;
        dt = dt_;

        n  = T_max/dt;
        increments = m.price_increments();
        market_prices = m.scenario();
    }

    std::vector<double> path(){

        std::vector<double> stock_path;
        double y = 0.0;

        for(int i=0; i<n; i++){

            y = rho*increments[i] + gen_norm(0,1)*sqrt(1.0-rho*rho);
            s = s*(1.0 + mu*dt + sigma*y*sqrt(dt));
            stock_path.push_back(s);
        }

        return stock_path;
    }
};
 
int main(){

    std::vector<double> fixed_times = {0, 91, 182, 273, 364};
    std::vector<double> fixed_prices = {100, 90, 110, 120, 100};
    market_scenario m = market_scenario(0.0, 0.04, 365, 1, fixed_times, fixed_prices);
    market_GRW mg = market_GRW(100, 0.8, 0.0, 0.05, 365, 1, m);

    output(mg.path(), "data");
    output(mg.market_prices, "market");

    return 0;
}