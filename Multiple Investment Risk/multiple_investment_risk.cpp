#include <iostream>
#include <random>
#include <fstream>
#include <string>
#include <vector>
#include <Eigen/Dense>
#include <optional>
 
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

struct GRW{
    /*
    geometric random walk with drift mu and volatility sigma
    T_max gives number of days for which to calculate the path
    s_0 is initial security price at day zero
    */

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


struct Correlation_Matrix{

    /*
    Given the correlations between pairs of variables and the individual standard deviations of each random variable, return the correlation matrix. 
    "correlations" should be a symmetric nxn matrix with ones on the diagonal where the (i,j)th entry is the correlation coefficient for variables i and j.
    The ith element of "std_deviations" should be the standard deviation for the ith random variable.

    A Hermitian, positive definite matrix M can be decomposed as M = LL^T, where L is lower triangular.
    This class also includes the matrix L, which is used in generating correlated samples.
    */

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

        std::cout << "C: " << C << std::endl;

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

Eigen::VectorXd correlated_samples(Eigen::MatrixXd L){
    /*
    Given a the matrix L obtained from the Cholesky decomposition of a correlation matrix C, generate a vector of correlated random variables. 
    L is generated by the "Correlation_matrix".
    */

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

struct market_scenario{
    /*
    We may want to calculate prices of securities that are correlated with some predetermined market conditions.
    We consider random fluctuations around a piecewise function constructed by linearly interpolating between a set of specified (time, price) points. 

    mu and sigma give drift and volatility, respectively.
    fixed_times holds the days at which we would like to specify market prices. fixed_prices holds the prices we specify. 
    It is important to ensure that units of fixed_times are consistent with the units of time used by the random walk. 

    adjusted_prices contains the values of the GRW at each of the t_points, because we want to interpolate between y_i-g(t), not y_i alone. Otherwise, we
    wouldn't actually have the desired market conditions.
    */

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
        
        // start with a zero price increment so that there are as many increments as there are time steps
        std::vector<double> PI = {0};

        for(int i=1; i<m.size(); i++){
            PI.push_back((m[i]-m[i-1])/(m[i-1]*sigma*sqrt(dt)));
        }

        return PI;
    }
};

//std::vector<Eigen::VectorXd> market_correlated_samples(Eigen::MatrixXd M, market_scenario m){
std::vector<Eigen::VectorXd> market_correlated_samples(Eigen::MatrixXd M, std::vector<double> increments){

    /*
    Creates an array of vector-valued samples that can be used to generate GRWs that are correlated with a given market scenario.
    Individual vector components correspond to the different securities that we consider. There is one vector in the array for each time step.

    L is the symmetric matrix with ones on the diagonal of correlation coefficients for the securities and the market scenario.
    Here, we take the market scenario to be the zeroth component of any vectors, so when defining the correlation matrix for the system, we must put the 
    correlations with the market scenario in the first row. 
    */

    //std::vector<double> increments = m.price_increments();
    std::vector<Eigen::VectorXd> Y;

    // Obtain the Cholesky decomposition of M, and along the way, check to make sure that we have a valid correlation matrix
    Eigen::LLT<MatrixXd> lltOfA(M);
    Eigen::MatrixXd L = lltOfA.matrixL();
    if(lltOfA.info() == Eigen::NumericalIssue){
        std::cout << "Correlation matrix is likely not positive semidefinite." << std::endl;
    }   

    for(int i=0; i<increments.size(); i++){

        // This vector holds the market scenario price increments as its zeroth value, and one N(0,1) independent random variable for each security.
        Eigen::VectorXd V(M.cols());

        // The first element of L*V is a random variable that is directly proportional to the market price increments. We don't actually need this
        // to construct the GRW for a security, so this vector will hold only the needed components, making it one element smaller than V. 
        Eigen::VectorXd Z(M.cols()-1);

        V(0) = increments[i];

        for(int j=1; j<M.cols(); j++){
            V(j) = gen_norm(0,1);
        }

        V = L*V;

        for(int j=1; j<V.size(); j++){
            Z(j-1) = V(j);
        }

        Y.push_back(Z);
    }

    return Y;
}

struct correlated_GRW{
    /*
    A vector-valued geometric random walk where components have correlated price increments. 

    mu is a vector that holds the drift parameter for each component. Similarly, sigma is the volatility vector. 
    T_max gives the number of days for which to calculate the path.
    s_0 is the initial price vector at day zero.
    correlations is a strictly upper triangular matrix containing the correlation coefficients rho_{ij}
    */

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

struct market_GRW{
    /*
    price of a security correlated with a given market scenario

    s gives initial security price, mu_s gives security drift, sigma_s gives security variance
    */

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


struct market_correlated_GRW{
    /*
    Create a portfolio of possibly correlated securites that are also correlated with a market scenario.

    s_0 contains the initial security values, mu contains the security drifts, and sigma contains the security variances. All of these quantities are 
    vector-valued, with a component for each security.

    correlations is a symmetric matrix with ones on the diagonal that contains the correlation coefficients rho_{ij}. We take rho_{0j} (using 
    zero-indexed matrix notation)to be the correlations of security j with the market scenario and rho_{ij} with i,j != 0 to be the correlation of 
    security i with security j (using one-indexed vector notation). To give an example, if we have two securities and the following correlation matrix,

        1.0, 0.1, 0.03,
        0.1, 1.0, 0.3,
        0.03, 0.3, 1.0;

    then security 1 has correlation 0.1 with the market and security 2 has correlation 0.03. The securities are correlated with each other with value 0.3.

    weights contains the relative proportions of each security in the portfolio. 
    */

    Eigen::VectorXd s_0, s, mu, sigma, z;
    double T_max = 0.0, dt = 0.0;
    int n = 0, steps_per_day = 0;

    std::vector<double> market_path, increments;
    std::vector<Eigen::VectorXd> samples;
    Eigen::MatrixXd correlations;
    Eigen::VectorXd weights;

    // I need a way to get the market scenario out of the constructor so that I can use it when I need to calculate the samples later on
    market_correlated_GRW(Eigen::VectorXd s_0_, Eigen::VectorXd mu_, Eigen::VectorXd sigma_, Eigen::MatrixXd correlations_, double T_max_, double dt_, market_scenario m, Eigen::VectorXd weights_){
        s_0 = s_0_;
        mu = mu_;
        sigma = sigma_;
        T_max = T_max_; 
        dt = dt_;
        n  = T_max/dt;
        steps_per_day = 1/dt;
        market_path = m.scenario();
        increments = m.price_increments();
        correlations = correlations_;
        weights = weights_;

        // check to make sure that the fractions of the overall portfolio made up by each security add to 1.
        double sum = 0;
        for(int i=0; i<weights.size(); i++){
            sum += weights(i);
        }

        // this may cause some issues with floating point precision, but as long as we don't have a portfolio where particular securities make up very
        // very small percentages of the overall portfolio, we should be okay. 
        if(sum != 1.0){
            std::cout << "Portfolio weights do not add to 1." << std::endl;
        }
    }

    // calculate prices as a function of time, measured in days
    // the first component of the vector at each time step should be the market value at that time
    std::vector<Eigen::VectorXd> path(){

        // calculate a new set of correlated samples
        samples = market_correlated_samples(correlations, increments);
        
        // create vector of VectorXds with size one greater than the number of stocks we want to consider. The first element should be the market 
        // scenario, as described above.
        std::vector<Eigen::VectorXd> data;
        
        // the market scenario and the correlated GRWs at a given time step
        Eigen::VectorXd values(s_0.size()+1);
        s = s_0;
        
        for(int i=0; i<n; i++){

            // samples needs to be calculated here, not referenced from somewhere else, that way each time we call path(), we get a new simulation
            z = samples[i];

            // for(int k=0; k<samples[i].size(); k++){
            //     std::cout << z(i) << std::endl;
            // }

            // calculate each vector component of the prices
            for(int j=0; j<s.size(); j++){
                s(j) = s(j)*(1+mu(j)*dt+sigma(j)*z(j)*sqrt(dt));
            }

            values(0) = market_path[i];

            for(int j=0; j<s.size(); j++){
                values(j+1) = s(j);
            }

            data.push_back(values);
        }

        return data;
    }

    // returns the likelihood that the portfolio value ends in a given range
    double risk(int num_trials, double min_val, double max_val){

        double counter = 0.0;
        double portfolio_value = 0.0;
        Eigen::VectorXd security_end_values(s_0.size());
        Eigen::VectorXd p(s_0.size());

        for(int i=0; i<num_trials; i++){

            p = path().back();

            // copy path end values (except market scenario value) into new vector of just the security values
            for(int j=1; j<p.size(); j++){
                security_end_values(j-1) = p(j);
            }            

            portfolio_value = weights.dot(security_end_values);    

            if(portfolio_value <= max_val && portfolio_value >= min_val){
                counter++;
            }
        }

        return counter/num_trials;
    }

    // create files for each component of the GRW and the market scenario in the Data subdirectory
    void output_mGRW(std::vector<Eigen::VectorXd> v){

        std::ofstream ofs;
        ofs.open("./Data/market_scenario.txt");
        
        for(int i=0; i<v.size(); i++){
            ofs << v[i](0) << '\n';
        }
        
        ofs.close();

        for(int i=1; i<v[0].size(); i++){

            //std::ofstream ofs;
            ofs.open("./Data/security_" + std::to_string(i) + ".txt");
            
            for(int j=0; j<v.size(); j++){
                ofs << v[j](i) << '\n';
            }
            
            ofs.close();
        }
    }
};
 
int main(){

    // number of securities in portfolio
    int num = 2;

    // parameters for market scenario
    std::vector<double> fixed_times = {0, 91, 182, 273, 364};
    std::vector<double> fixed_prices = {100, 90, 110, 120, 100};
    double market_drift = 0.0;
    double market_variance = 0.04;

    market_scenario m = market_scenario(market_drift, market_variance, 365, 1, fixed_times, fixed_prices);

    // first security has correlation 0.7 with the market, second has correlation 0.1. The securities are also correlated with each other with value 0.4.
    Eigen::MatrixXd correlations(num+1, num+1);
    correlations << 1.0, 0.8, 0.9,
                    0.8, 1.0, 0.6,
                    0.9, 0.6, 1.0; 

    // initial conditions for both securities
    Eigen::VectorXd s(num);
    s << 100.0, 50.0;

    // drifts for each security
    Eigen::VectorXd mu(num);
    mu << 0.0, 0.0;

    // Variances for the market and each security. We need to include the market variance because it is used when calculating the correlation matrix.
    Eigen::VectorXd sigma(num+1);
    sigma << market_variance, 0.04, 0.06;

    // portfolio consists of 60% security 1 and 40% security 2
    Eigen::VectorXd weights(num);
    weights << 0.6, 0.4;

    market_correlated_GRW mGRW = market_correlated_GRW(s, mu, sigma, correlations, 365, 1, m, weights);

    std::vector<Eigen::VectorXd> p = mGRW.path();
    mGRW.output_mGRW(p);

    // total amount of money that we need to spend to buy the portfolio at t=0
    double starting_val = weights.dot(s);

    // probability of making money
    double r = mGRW.risk(10000, starting_val, 1e300);

    std::cout << r << std::endl;

    return 0;
}