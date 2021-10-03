#include <iostream>
#include <random>
#include <string>
#include <vector>
#include <fstream>
#include <Eigen/Dense>
#include <math.h>

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

double N(double x){
    /*
    Cumulative normal function
    */
    return 0.5*(1+erf(x/sqrt(2)));
}

double N_prime(double x){
    /*
    Derivative of N(x)
    */
    return exp(-0.5*x*x)/sqrt(2*M_PI);
}

class Option{

    public:
    double S, K, r, d, T, sigma;

    Option(double S_, double K_, double r_, double d_, double T_, double sigma_) : S(S_), K(K_), r(r_), d(d_), T(T_), sigma(sigma_){}

    virtual double payout(double x){
        /*
        A generic payout function that will later be replaced by a class for a specific option type
        */
        return 0;
    }

    double d1(double S_, double K_, double r_, double d_, double T_, double sigma_){
        return (log(S_/K_)+(r_-d_+0.5*sigma_*sigma_)*T_)/(sigma_*sqrt(T_));
        
    }

    double d2(double S_, double K_, double r_, double d_, double T_, double sigma_){
        return (log(S_/K_)+(r_-d_-0.5*sigma_*sigma_)*T_)/(sigma_*sqrt(T_));
    }

    double MC_SDE_underlying_price(double z){
        /*
        Return the end value of a lognormally distributed path, as determined by solving the SDE. 
        z should be a N(0,1) distributed random variable, as calculated by the RNG class
        */
        return S*exp((r-d)*T-0.5*sigma*sigma*T+sigma*z*sqrt(T));
    }

    double MC_price(std::vector<double> rand_vals){
        /*
        Return the price of a call using the SDE price as calculated above. 
        rand_vals should contain a number of N(0,1) random variables equal to the number of trials that we want to perform. 
        */
        double sum = 0, price = 0;
        int n = rand_vals.size();

        for(int i=0; i<n; i++){
            price = MC_SDE_underlying_price(rand_vals[i]);
            sum += payout(price);
        }

        return exp(-r*T)*sum/n;
    }

    double variance(std::vector<double> data){
        double mean = 0, var = 0;

        for(double x : data){
            mean += x;
        }

        mean = mean/data.size();

        for(double x : data){
            var += (x-mean)*(x-mean);
        }

        var = var/data.size();

        return var;
    }

    double std_error(double std_dev, double n){
        return std_dev/sqrt(n);
    }
};

class Call : public Option{
    public:

    using Option::Option;

    double payout(double s){
        return std::max(s-K, 0.0);
    }

    double d_payout_ds(double s){
        if(S >= K){
            return 1;
        }

        else{
            return 0;
        }
    }

    //==================================================================================================================================================
    //          Properties derived from the Black-Scholes equation
    //==================================================================================================================================================


    double BS_price(double S_, double K_, double r_, double d_, double T_, double sigma_){     
        double D1 = d1(S_, K_, r_, d_, T_, sigma_);
        double D2 = d2(S_, K_, r_, d_, T_, sigma_);
        return S_*exp(-d_*T_)*N(D1)-K*exp(-r_*T_)*N(D2);
    }

    double BS_Delta(double S_, double K_, double r_, double d_, double T_, double sigma_){
        double D1 = d1(S_, K_, r_, d_, T_, sigma_);
        return exp(-d*T)*N(D1);
    }

    double BS_Gamma(double S_, double K_, double r_, double d_, double T_, double sigma_){
        double D1 = d1(S_, K_, r_, d_, T_, sigma_);
        return exp(-d*T)*N_prime(D1)/(S*sigma*sqrt(T));
    }

    // Also known as Vega, if you don't like the Greek alphabet.
    double BS_Kappa(double S_, double K_, double r_, double d_, double T_, double sigma_){
        double D1 = d1(S_, K_, r_, d_, T_, sigma_);
        return exp(-d*T)*S*sqrt(T)*N_prime(D1);
    }

    double BS_Rho(double S_, double K_, double r_, double d_, double T_, double sigma_){
        double D2 = d2(S_, K_, r_, d_, T_, sigma_);
        return K*exp(-d*T)*T*N(D2);
    }

    double BS_Theta(double S_, double K_, double r_, double d_, double T_, double sigma_){
        double D1 = d1(S_, K_, r_, d_, T_, sigma_);
        double D2 = d2(S_, K_, r_, d_, T_, sigma_);
        return S*exp(-d*T)*d*N(D1) - K*exp(-r*T)*r*N(D2) - S*exp(-r*T)*sigma*N_prime(D1)/(2*sqrt(T));
    }

    //==================================================================================================================================================
    //          Properties obtained by using finite differences
    //==================================================================================================================================================


    double FD_Delta(double S_, double K_, double r_, double d_, double T_, double sigma_, double epsilon){
        return (BS_price(S+epsilon, K, r, d, T, sigma)-BS_price(S, K, r, d, T, sigma))/epsilon;
    }

    double FD_Gamma(double S_, double K_, double r_, double d_, double T_, double sigma_, double epsilon){
        return (BS_price(S_+epsilon, K_, r_, d_, T_, sigma_) - 2*BS_price(S_, K_, r_, d_, T_, sigma) + BS_price(S_-epsilon, K_, r_, d_, T_, sigma_))/(epsilon*epsilon);
    }

    double FD_Kappa(double S_, double K_, double r_, double d_, double T_, double sigma_, double epsilon){
        return (BS_price(S_, K_, r_, d_, T_, sigma_+epsilon)-BS_price(S_, K_, r_, d_, T_, sigma_))/epsilon;
    }

    double FD_Rho(double S_, double K_, double r_, double d_, double T_, double sigma_, double epsilon){
        return (BS_price(S_, K_, r_+epsilon, d_, T_, sigma_)-BS_price(S_, K_, r_, d_, T_, sigma_))/epsilon;
    }

    double FD_Theta(double S_, double K_, double r_, double d_, double T_, double sigma_, double epsilon){
        return -(BS_price(S_, K_, r_, d_, T_+epsilon, sigma_)-BS_price(S_, K_, r_, d_, T_, sigma_))/epsilon;
    }

    //==================================================================================================================================================
    //          Properties obtained by using Monte Carlo
    //==================================================================================================================================================

    double MC_price(double S_, double K_, double r_, double d_, double T_, double sigma_, std::vector<double> rand_vals){
        /*
        Return the price of a call using the SDE price as calculated above. 
        rand_vals should contain a number of N(0,1) random variables equal to the number of trials that we want to perform. 
        */
        double sum = 0, price = 0;
        int n = rand_vals.size();
        double z;

        for(int i=0; i<n; i++){
            z = rand_vals[i];
            price = S_*exp((r_-d_)*T_-0.5*sigma_*sigma_*T_+sigma_*z*sqrt(T_));
            sum += payout(price);
        }

        return exp(-r_*T_)*sum/n;
    }

    double MC_Delta(double S_, double K_, double r_, double d_, double T_, double sigma_, double epsilon, std::vector<double> rand_vals_1, std::vector<double> rand_vals_2){
        double price_1 = MC_price(S_ + epsilon, K_, r_, d_, T_, sigma_, rand_vals_1);
        double price_2 = MC_price(S_, K_, r_, d_, T_, sigma_, rand_vals_2);

        return (price_1-price_2)/epsilon;
    }

    double MC_Delta_2(double S_, double K_, double r_, double d_, double T_, double sigma_, double epsilon, std::vector<double> rand_vals){
        double price_1 = MC_price(S_ + epsilon, K_, r_, d_, T_, sigma_, rand_vals);
        double price_2 = MC_price(S_, K_, r_, d_, T_, sigma_, rand_vals);

        return (price_1-price_2)/epsilon;
    }

    double dlog_phi(double x, double S_0_, double r_, double d_, double T_, double sigma_){
        /*
        Derivative of the logarithm of the Black-Scholes distribution with respect to S.
        This is used in the likelihood and pathwise methods.
        */
        
        return -(1.0 + (log(x/S_0_)-T_*(r_ - d_ - 0.5*sigma_*sigma_))/(T_*sigma_*sigma_))/x;
    }
};

int main(){

    double S = 100, K = 90, r = 0.03, d = 0.01, T = 60/365., sigma = 0.3, epsilon=0.01;
    Call C = Call(S, K, r, d, T, sigma);

    // Check to see that the Black-Scholes Greeks and the finite difference Greeks agree
    std::cout << "Delta: " << C.BS_Delta(S, K, r, d, T, sigma) - C.FD_Delta(S, K, r, d, T, sigma, epsilon) << std::endl <<
                 "Gamma: " << C.BS_Gamma(S, K, r, d, T, sigma) - C.FD_Gamma(S, K, r, d, T, sigma, epsilon) << std::endl <<
                 "Kappa: " << C.BS_Kappa(S, K, r, d, T, sigma) - C.FD_Kappa(S, K, r, d, T, sigma, epsilon) << std::endl <<
                 "Rho: " << C.BS_Rho(S, K, r, d, T, sigma) - C.FD_Rho(S, K, r, d, T, sigma, epsilon) << std::endl << 
                 "Theta: " << C.BS_Theta(S, K, r, d, T, sigma) - C.FD_Theta(S, K, r, d, T, sigma, epsilon) << std::endl;

    //==================================================================================================================================================
    //          Plotting
    //==================================================================================================================================================

    // Plot as a function of spot price
    int n = 1000;
    double S_max = 100, dS=S_max/n, Delta, Gamma, Kappa;
    std::vector<double> Delta_vals_spot, Gamma_vals_spot, Kappa_vals_spot, S_vals;
    S = 0;

    for(int i=0; i<n; i++){
        Delta = C.BS_Delta(S, K, r, d, T, sigma);
        Gamma = C.BS_Gamma(S, K, r, d, T, sigma);
        Kappa = C.BS_Kappa(S, K, r, d, T, sigma);
        
        Delta_vals_spot.push_back(Delta);
        Gamma_vals_spot.push_back(Gamma);
        Kappa_vals_spot.push_back(Kappa);

        S_vals.push_back(S);

        S += dS;
    }

    output(Delta_vals_spot, "Delta_spot");
    output(Gamma_vals_spot, "Gamma_spot");
    output(Kappa_vals_spot, "Kappa_spot");
    output(S_vals, "S");

    // Plot as a function of time
    double t=0, Expiry_time=T, dt=Expiry_time/n;
    std::vector<double> Delta_OTM_vals, Delta_ATM_vals, Delta_ITM_vals, Kappa_vals_time, t_vals;

    for(int i=0; i<n; i++){
        
        // ITM Delta
        S = 100;
        Delta = C.BS_Delta(S, K, r, d, T, sigma);
        Delta_ITM_vals.push_back(Delta);

        // ATM Delta
        S = 90;
        Delta = C.BS_Delta(S, K, r, d, T, sigma);
        Delta_ATM_vals.push_back(Delta);

        // OTM Delta
        S = 80;
        Delta = C.BS_Delta(S, K, r, d, T, sigma);
        Delta_OTM_vals.push_back(Delta);

        // Kappa
        Kappa = C.BS_Kappa(S, K, r, d, T, sigma);
        Kappa_vals_time.push_back(Kappa);

        t_vals.push_back(t*365);

        t += dt;
        T = Expiry_time - t;
    }

    output(Delta_ITM_vals, "Delta_ITM");
    output(Delta_ATM_vals, "Delta_ATM");
    output(Delta_OTM_vals, "Delta_OTM");
    output(Kappa_vals_time, "Kappa_time");
    output(t_vals, "t");

    // Plot as a function of volatility
    double sigma_max=1, d_sig=sigma_max/n;
    std::vector<double> Kappa_vals_volatility, sigma_vals;
    T=60/365., S=100, sigma=0;

    for(int i=0; i<n; i++){
        Kappa = C.BS_Kappa(S, K, r, d, T, sigma);
        Kappa_vals_volatility.push_back(Kappa);

        sigma_vals.push_back(sigma);
        sigma += d_sig;
    }

    output(Kappa_vals_volatility, "Kappa_vol");
    output(sigma_vals, "sigma");

    //==================================================================================================================================================
    //          Monte Carlo Greeks
    //==================================================================================================================================================

    // First, calculate Delta using Monte Carlo for finite differences
    RNG R = RNG();
    int m = 20;
    sigma = 0.3;

    std::vector<double> rand_vals_1, rand_vals_2, Delta_1_vals, Delta_2_vals;
    double Delta_1, Delta_2;

    for(int i=0; i<m; i++){
        n = std::pow(2,i);
        
        for(int j=0; j<n; j++){
            rand_vals_1.push_back(R.gen_norm(0,1));
            rand_vals_2.push_back(R.gen_norm(0,1));
        }

        // Different RNG
        Delta_1 = C.MC_Delta(S, K, r, d, T, sigma, epsilon, rand_vals_1, rand_vals_2);
        
        // Same RNG
        Delta_2 = C.MC_Delta_2(S, K, r, d, T, sigma, epsilon, rand_vals_1);

        Delta_1_vals.push_back(Delta_1);
        Delta_2_vals.push_back(Delta_2);
    }

    output(Delta_1_vals, "Delta_1_vals");
    output(Delta_2_vals, "Delta_2_vals");

    // Next, the likelihood ratio and pathwise methods
    double  u, sum1=0, sum2=0, z, S_T;

    std::vector<double> likelihood_vals, pathwise_vals;

    for(int i=0; i<m; i++){
        n = std::pow(2,i);
        sum1 = 0;
        sum2 = 0;

        std::cout << "i: " << i << std::endl;

        for(int j=0; j<n; j++){
            z = R.gen_norm(0,1);
            S_T = C.MC_SDE_underlying_price(z);
            sum1 += exp(-r*T)*C.payout(C.MC_SDE_underlying_price(z))*z/(S*sigma*sqrt(T))/n;
            sum2 += exp(-r*T)*S_T*C.d_payout_ds(S_T)/(S*n);
        }

        likelihood_vals.push_back(sum1);
        pathwise_vals.push_back(sum2);
    }

    output(likelihood_vals, "likelihood_vals");
    output(pathwise_vals, "pathwise_vals");

    std::cout<< "Black-Scholes value is: " << C.BS_Delta(S, K, r, d, T, sigma) << std::endl;
}