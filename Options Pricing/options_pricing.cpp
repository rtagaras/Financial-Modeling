#include <iostream>
#include <random>
#include <string>
#include <vector>
#include <Eigen/Dense>

using Eigen::MatrixXd;

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

struct GRW{
    /*
    geometric random walk with drift mu and volatility sigma
    T_max gives number of time units for which to calculate the path, dt is number of steps per unit time
    s_0 is initial security price at day zero
    */

    double mu, sigma, T_max, dt, s_0, s, z;
    int n;

    GRW(double s_0_, double mu_, double sigma_, double T_max_, double dt_){
        s_0 = s_0_;
        mu = mu_;
        sigma = sigma_;
        T_max = T_max_; 
        dt = dt_;
        n  = T_max/dt;
    }

    // calculate price as a function of time, measured in years
    std::vector<double> path(){
        std::vector<double> data;
        s = s_0;
        
        for(int i=0; i<n; i++){

            z = gen_norm(0,1);
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

struct correlated_GRW{
    /*
    A vector-valued geometric random walk where components have correlated price increments. 

    mu is a vector that holds the drift parameter for each component. Similarly, sigma is the volatility vector. 
    T_max gives the number of days for which to calculate the path.
    s_0 is the initial price vector at day zero.
    correlations is a symmetric diagonal matrix containing the correlation coefficients rho_{ij}
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

        Correlation_Matrix CM = Correlation_Matrix(Eigen::VectorXd::Ones(sigma.size()), correlations_);

        Corr = CM.C;
        l = CM.L;
    }

    // calculate price as a function of time, measured in days
    std::vector<Eigen::VectorXd> path(){
        
        std::vector<Eigen::VectorXd> data;
        s = s_0;

        for(int i=0; i<n; i++){

            z = correlated_samples(l);
            s = s.cwiseProduct((Eigen::VectorXd::Ones(s_0.size()) + mu*dt + sigma.cwiseProduct(z)*sqrt(dt)));
            data.push_back(s);
        }

        return data;
    }
};

class Option{
    /*
    Base class for a general option

    s_0 is initial underlying value, r is risk-free rate (measured in percent/year), mu is drift, sigma is volatility, expiry time is measured in days
    */
    protected:
    double strike, s_0, r, mu, sigma, dt, expiry_time;
    std::string type;

    public:
    Option(std::string type_, double strike_, double expiry_time_, double dt_, double s_0_, double r_, double mu_, double sigma_){
        type = type_;
        strike = strike_;
        expiry_time = expiry_time_;
        s_0 = s_0_;
        dt = dt_;
        r = r_;
        mu = mu_;
        sigma = sigma_;
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
};

class European_Option : protected Option{
    
    public:
    using Option::Option;

    // We calculate the expected payoff and discount the price by using the risk free rate, following the risk-neutral principal. This is computationally
    // simpler than implementing a recursive solution, although potentially less flexible. 
    double binomial_lattice_price(int num_trials){

        double E = 0.0, x = 0.0;
        
        // start by using the u=1/d convention for binomial factors; if this leads to unphysical results, use p=1/2 convention instead
        double A = (exp(-mu*dt)+exp((mu+sigma*sigma)*dt))/2.;
        double d = A - sqrt(A*A - 1.0);
        double u = A + sqrt(A*A - 1.0);
        double p = (exp(mu*dt)-d)/(u-d);

        if(p <= 0 || p >= 1){
            d = exp(mu*dt)*(1.0 - sqrt(exp(sigma*sigma*dt) - 1.0));
            u = exp(mu*dt)*(1.0 + sqrt(exp(sigma*sigma*dt) - 1.0));
            p = 0.5;
        }

        // calculate many paths through the lattice, adding the final payout to a running total
        for(int i=0; i<num_trials; i++){
            double s = s_0;

            for(int j=0; j<expiry_time/dt; j++){
                x = gen_uniform(0,1);

                if(x > p){
                    s = s*u;
                }

                else{
                    s = s*d;
                }
            }

            E += payout(s);
        }

        // return the average payout, discounted by the RFR to t=0
        E = E/num_trials;
        return exp(-expiry_time*r*dt)*E;
    }

    // prices should follow geometric brownian motion - use a GRW to simulate end price
    double GRW_price(int n){
        
        double x = 0.0;
        std::vector<double> p;
        GRW g = GRW(s_0, mu, sigma, expiry_time, dt);
        
        for(int i=0; i<n; i++){
            
            p = g.path();
            x += payout(p.back());
        }

        return exp(-r*expiry_time)*x/n;
    }
};

class American_Option : protected Option{
    
    public:
    using Option::Option;

    double binomial_lattice_price(){

        double x = 0.0;
        int n = expiry_time/dt;
        
        // start by using the u=1/d convention for binomial factors; if this leads to unphysical results, use p=1/2 convention instead
        double A = (exp(-mu*dt)+exp((mu+sigma*sigma)*dt))/2.;
        double d = A - sqrt(A*A - 1.0);
        double u = A + sqrt(A*A - 1.0);
        double p = (exp(mu*dt)-d)/(u-d);

        //if(p <= 0 || p >= 1){
            d = exp(mu*dt)*(1.0 - sqrt(exp(sigma*sigma*dt) - 1.0));
            u = exp(mu*dt)*(1.0 + sqrt(exp(sigma*sigma*dt) - 1.0));
            p = 0.5;
        //}
    
        Eigen::MatrixXd underlying_lattice(n+1, n+1);
        Eigen::MatrixXd option_lattice(n+1, n+1);

        // calculate underlying values
        for(int i=0; i<n+1; i++){
            for(int j=0; j<i+1; j++){
                underlying_lattice(i,j) = s_0*pow(u,j)*pow(d,i-j);
            }
        }

        // calculate all possible option payouts
        for(int j=0; j<n+1; j++){
            option_lattice(n, j) = payout(underlying_lattice(n, j));
        }

        // backpropagate option value to t=0
        for(int i=n-1; i>=0; i--){
            for(int j=0; j<i+1; j++){

                //calculate option price from child nodes and compare to value obtained by exercising the option
                x =  exp(-r*dt)*(p*option_lattice(i+1,j+1) + (1-p)*option_lattice(i+1,j));
                option_lattice(i,j) = std::max(x, payout(underlying_lattice(i,j)));
            }
        }

        return option_lattice(0,0);
    }

};

class Asian_Option : protected Option{
    public:
    using Option::Option;

    double fixed_strike_GRW_price(int n){
        double x = 0.0, s = 0.0;
        std::vector<double> p;

        for(int i=0; i<n; i++){
            GRW g = GRW(s_0, mu, sigma, expiry_time, dt);

            s = 0.0;
            p = g.path();

            // get average price over entire random walk
            for(int j=0; j<p.size(); j++){
                s += p[j]/p.size();
            }

            x += payout(s);
        }

        return exp(-r*expiry_time)*x/n;
    }
};

class Barrier_Option : protected Option{

    private:
    double E = 0.0, s = 0.0, b, d1 = 0.0, d2 = 0.0, z, u;
    int n;
    bool crossed = 0;

    public:
    Barrier_Option(std::string type_, double strike_, double barrier_, double expiry_time_, double dt_, double s_0_, double r_, double mu_, double sigma_) 
        : Option(type_, strike_, expiry_time_, dt_, s_0_, r_, mu_, sigma_){
        
        b = barrier_;
        n = expiry_time/dt;
    }

    double GRW_price(int num_trials){
        for(int i=0; i<num_trials; i++){

            crossed = 0;
            s = s_0;
            d1 = b-s_0;

            for(int j=0; j<n; j++){

                // normal GRW calculation
                z = gen_norm(0,1);
                s = s*(1.0 + mu*dt + sigma*z*sqrt(dt));
                d2 = b-s;

                // use Brownian bridge to decide if barier has been crossed, which makes option value zero
                u = gen_uniform(0,1);
                if(u < exp(-2.*d1*d2/(sigma*sigma*dt))){
                    crossed = 1;
                }

                d1 = d2;
            }

            if(crossed == 0){
                E += payout(s);
            }
        }

        return E*exp(-r*expiry_time)/num_trials;
    }
};

class Basket_Option{
    private:

    std::string type;
    double strike, dt, r, expiry_time;
    Eigen::MatrixXd correlations;
    Eigen::VectorXd weights, s_0, mu, sigma;
    std::vector<Eigen::VectorXd> p;

    public:
    Basket_Option(std::string type_, double strike_, Eigen::VectorXd weights_, Eigen::MatrixXd correlations_, double expiry_time_, double dt_, Eigen::VectorXd s_0_, double r_, Eigen::VectorXd mu_, Eigen::VectorXd sigma_){

        type = type_;
        strike = strike_;
        weights = weights_;
        correlations = correlations_;
        expiry_time = expiry_time_;
        dt = dt_;
        s_0 = s_0_;
        r = r_;
        mu = mu_;
        sigma = sigma_;
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

    double GRW_value(int n){

        double E = 0.0, x = 0.0;
        correlated_GRW g = correlated_GRW(s_0, mu, sigma, correlations, expiry_time, dt);
        
        for(int i=0; i<n; i++){
                
            p = g.path();   
            x = weights.dot(p.back());
            E += payout(x);
        }

        return E*exp(-r*expiry_time)/n;
    }
};

class Exchange_Option {
    private:
    double A_0, B_0, correlation, mu1, mu2, sigma1, sigma2, r, dt, expiry_time;
    std::vector<double> p1, p2;
    int n;

    public:
    Exchange_Option(double A_0_, double B_0_, double correlation_, double mu1_, double mu2_, double sigma1_, double sigma2_, double r_, double dt_, double expiry_time_){
        correlation = correlation_;
        expiry_time = expiry_time_;
        dt = dt_;
        A_0 = A_0_;
        B_0 = B_0_;
        r = r_;
        mu1 = mu1_;
        mu2 = mu2_;
        sigma1 = sigma1_;
        sigma2 = sigma2_;
        n = expiry_time/dt;
    }

    double payout(double x){
        return std::max(x, 0.0);
    }

    double GRW_value(int num_trials){
        
        double E = 0.0, A, B, z1, z2, x;
        for(int i=0; i<num_trials; i++){
            A = A_0;
            B = B_0;

            for(int j=0; j<n; j++){
                z1 = gen_norm(0,1);
                z2 = gen_norm(0,1);
                x = correlation*z1 + sqrt(1-correlation*correlation)*z2;

                A += A*(mu1*dt + sigma1*z1*sqrt(dt));
                B += B*(mu2*dt + sigma2*x*sqrt(dt)); 
            }
            
            E += payout(A-B);
        }
        return E*exp(-r*expiry_time)/num_trials;
    }
    
};

int main(){

    // general parameters
    std::string type = "put";
    double s_0 = 100;
    double strike = 100;
    double r = 0.03;
    double sigma = 0.2;
    double expiry_time = 60/365.;
    double dt = 0.01/365.;
    double barrier = 105;

    // Parameters for basket option
    Eigen::VectorXd weights(3);
    weights << 0.33, 0.33, 0.33;

    Eigen::MatrixXd correlations(3,3);
    correlations << 1.0, 0.7, 0.3,
                    0.7, 1.0, -0.1,
                    0.3, -0.1, 1.0;

    Eigen::VectorXd s_0_B(3);
    s_0_B << 100, 100, 100;

    Eigen::VectorXd mu(3);
    mu << 0.03, 0.03, 0.03;

    Eigen::VectorXd sigma_B(3);
    sigma_B << 0.2, 0.2, 0.2;



    European_Option E1 = European_Option(type, strike, expiry_time, dt, s_0, r, r, sigma);
    std::cout << "European GRW price: " << E1.GRW_price(1000) << std::endl;

    European_Option E2 = European_Option(type, 49, 4/365., dt, 50, 0.26, 0.26, 0.4);
    std::cout << "European lattice price: " << E2.binomial_lattice_price(10000) << std::endl;

    American_Option A = American_Option(type, 49, 4/365., dt, 50, 0.26, 0.26, 0.4);
    std::cout << "American lattice price: " << A.binomial_lattice_price() << std::endl;

    Asian_Option As = Asian_Option(type, strike, expiry_time, dt, s_0, r, r, sigma);
    std::cout << "Asian GRW price: " << As.fixed_strike_GRW_price(100) << std::endl;

    Barrier_Option B = Barrier_Option(type, strike, barrier, expiry_time, dt, s_0, r, r, sigma);    
    std::cout << "Knock-out option GRW price: " << B.GRW_price(100) << std::endl;

    Basket_Option Ba = Basket_Option(type, strike, weights, correlations, expiry_time, dt, s_0_B, r, mu, sigma_B);
    std::cout << "Basket GRW price: " << Ba.GRW_value(100) << std::endl;

    Exchange_Option E = Exchange_Option(100, 100, 0.6, -0.05, -0.03, 0.2, 0.2, 0.03, 1/365., 90/365.);
    std::cout << "Exchange GRW price: " << E.GRW_value(10000) << std::endl;


    // KNOWN VALUES
    //
    //      OPTION      | VALUE | CORRECT? 
    //=====================================
    // European GRW     | 2.99  |     Y
    // European Lattice | 0.421 |     Y
    // American Lattice | 0.425 |     Y
    // Asian GRW        | 1.75  |     Y
    // Barrier GRW      | 2.57  |     Y
    // Basket GRW       | 2.62  |     Y
    // Exchange GRW     | 3.24  |     Y

    return 0;
}