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

class Option{

    private:
    double S_0, r, d, T, sigma;

    public:
    double K;

    Option(double S_0_, double K_, double r_, double d_, double T_, double sigma_) : S_0(S_0_), K(K_), r(r_), d(d_), T(T_), sigma(sigma_) {}

    virtual double payout(double x){
        /*
        A generic payout function that will later be replaced by a class for a specific option type
        */
        return 0;
    }

    double MC_SDE_price(double z){
        /*
        Return the end value of a lognormally distributed path, as determined by solving the SDE. 
        z should be a N(0,1) distributed random variable, as calculated by the RNG class
        */
        return S_0*exp((r-d)*T-0.5*sigma*sigma*T+sigma*z*sqrt(T));
    }

    double MC_price(std::vector<double> rand_vals){
        /*
        Return the price of a call using the SDE price as calculated above. 
        rand_vals should contain a number of N(0,1) random variables equal to the number of trials that we want to perform. 
        */
        double sum = 0, price = 0;
        int n = rand_vals.size();

        for(int i=0; i<n; i++){
            price = MC_SDE_price(rand_vals[i]);
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

    double BS_price(double S, double K, double sigma, double r, double T, double d){
        /*
        Return Black-Scholes value for a call
        */
        double d1 = (log(S/K)+(r-d+0.5*sigma*sigma)*T)/(sigma*sqrt(T));
        double d2 = (log(S/K)+(r-d-0.5*sigma*sigma)*T)/(sigma*sqrt(T));

        return S*exp(-d*T)*N(d1)-K*exp(-r*T)*N(d2);
    }
};

class Put : public Option{
    public:

    using Option::Option;

    double payout(double s){
        return std::max(K-s, 0.0);
    }

    double BS_price(double S, double K, double sigma, double r, double T, double d){
        /*
        Return Black-Scholes value for a put
        */
        double d1 = (log(S/K)+(r-d+0.5*sigma*sigma)*T)/(sigma*sqrt(T));
        double d2 = (log(S/K)+(r-d-0.5*sigma*sigma)*T)/(sigma*sqrt(T));

        return -S*exp(-d*T)*N(-d1)+K*exp(-r*T)*N(-d2);
    }

};

class Digital_Call : public Option{
    public:

    using Option::Option;

    double payout(double s){
        return s>K ? 1 : 0;
    }

    double BS_price(double S, double K, double sigma, double r, double T, double d){
        /*
        Return Black-Scholes value for a digital call
        */
        double d2 = (log(S/K)+(r-d-0.5*sigma*sigma)*T)/(sigma*sqrt(T));

        return exp(-r*T)*N(d2);
    }
};

class Digital_Put : public Option{
    public:

    using Option::Option;

    double payout(double s){
        return s<K ? 1 : 0;
    }

    double BS_price(double S, double K, double sigma, double r, double T, double d){
        /*
        Return Black-Scholes value for a digital put
        */
        double d2 = (log(S/K)+(r-d-0.5*sigma*sigma)*T)/(sigma*sqrt(T));

        return exp(-r-T)*N(-d2);
    }
};

class Forward{

    private:
    double S_0, K, r, d, T;
    
    public:

    Forward(double S_0_, double K_, double r_, double d_, double T_) : S_0(S_0_), K(K_), r(r_), d(d_), T(T_) {}

    double BS_price(double T, double r, double d, double S_0, double K){
        return exp(-r*T)*(exp((r-d)*T)*S_0-K);
    }
};

class Zero_Coupon_Bond{
   
    private:
    double S_0, r, d, T;

    public:

    Zero_Coupon_Bond(double S_0_, double r_, double d_, double T_) : S_0(S_0_), r(r_), d(d_), T(T_) {}

    double BS_price(double S_0, double r, double d, double T){
        /*
        Return value of a zero-coupon bond
        */
        return S_0*exp((r-d)*T);
    }
};

int main(){

  /*
  For consistency, we should check the following:

    1. BS_call - BS_put = forward_contract_price at all times
    2. BS_call should be monotonically descreasing as a function of K
    3. BS_call should return a value between S and S-K*exp(-rT)
    4. BS_call should be monotonically increasing as a function of volatility
    5. If d=0, BS_call should be an increasing function of time to expiry
    6. BS_call should be a convex function of K
    7. The price of a call spread should approximate BS_digital_call
    8. BS_digital_call + BS_digital_put = zero_coupon_bond
  */

    // parameters used in multiple checks. These were chosen essentially at random. 
    int n = 1000;
    double Expiry_time=100, r=0.03, sigma=0.2, K=100, S=100, d=0.01;


    // //==================================================================================================================================================
    // //      check 1 - BS_call - BS_put = forward_contract_price at all times
    // //==================================================================================================================================================

    double T=0, dT=Expiry_time/n, Time_to_expiry=Expiry_time-T;
    double diff;
    bool fail = 0;

    // Specify the tolerance for failure. In theory, diff is always exactly zero, but floating point errors may not allow for this. 
    double epsilon = 0.0001;

    Call C = Call(S, K, r, d, Time_to_expiry, sigma);
    Put P = Put(S, K, r, d, Time_to_expiry, sigma);
    Forward F = Forward(S,K,r,d,Time_to_expiry);

    for(int i=0; i<n; i++){

        // this should be zero for all times
        diff = C.BS_price(S,K,sigma,r,Time_to_expiry,d) - P.BS_price(S,K,sigma,r,Time_to_expiry,d) - F.BS_price(Time_to_expiry, r, d, S, K);

        if(diff > epsilon){
            fail = 1;
        }

        T += dT;
        Time_to_expiry = Expiry_time - T;
    }

    if(!fail){
        std::cout << "All cases passed for check 1." << std::endl;
    }

    // reset time to expiry for next test
    Time_to_expiry = 100;

    // //==================================================================================================================================================
    // //      check 2 - BS_call should be monotonically descreasing as a function of K
    // //==================================================================================================================================================

    double K_max=10000, K_min=0, K_step = (K_max-K_min)/n;
    K = K_min;

    double prev_call_val = C.BS_price(S,K,sigma,r,Time_to_expiry,d);
    double call_val=0;
    fail = 0;

    for(int i=0; i<n; i++){
        K += K_step;
        call_val = C.BS_price(S,K,sigma,r,Time_to_expiry,d);

        if(call_val - prev_call_val > 0){
            fail = 1;
        }

        prev_call_val = call_val;

    }

    if(!fail){
        std::cout << "All cases passed for check 2." << std::endl;
    }

    // reset K for next test
    K = 100;

    // //==================================================================================================================================================
    // //      check 3 - BS_call should return a value between S*exp(-d*T) and S*exp(-d*T)-K*exp(-rT)
    // //==================================================================================================================================================

    fail = 0;
    call_val = C.BS_price(S, K, sigma, r, Time_to_expiry, d);

    if(call_val>S*exp(-d*Time_to_expiry) || call_val<S*exp(-d*Time_to_expiry)-K*exp(-r*Time_to_expiry)){
        fail = 1;
    }

    if(!fail){
        std::cout << "All cases passed for check 3." << std::endl;
    }

    // //==================================================================================================================================================
    // //      check 4 - BS_call should be monotonically increasing as a function of volatility
    // //==================================================================================================================================================

    double sigma_max = 1.0, sigma_step = sigma_max/n;
    sigma = 0;
    fail = 0;
    prev_call_val = C.BS_price(S,K,sigma,r,Time_to_expiry,d);

    for(int i=0; i<n; i++){
        sigma += sigma_step;
        call_val = C.BS_price(S,K,sigma,r,Time_to_expiry,d);

        if(call_val - prev_call_val < 0){
            fail = 1;
        }

        prev_call_val = call_val;
    }

    if(!fail){
        std::cout << "All cases passed for check 4." << std::endl;
    }

    // reset volatility for next test
    sigma = 0.3;

    // //==================================================================================================================================================
    // //      check 5 - If d=0, BS_call should be an increasing function of time to expiry
    // //==================================================================================================================================================

    T = Time_to_expiry;
    dT = Expiry_time/n;
    Time_to_expiry = Expiry_time-T;
    fail = 0;
    d=0;

    prev_call_val = C.BS_price(S,K,sigma,r,Time_to_expiry,d);
    fail = 0;

    for(int i=0; i<n; i++){
        T += dT;
        Time_to_expiry = Expiry_time-T;
        call_val = C.BS_price(S,K,sigma,r,Time_to_expiry,d);

        if(call_val <= prev_call_val){
            fail = 1;
        }

        prev_call_val = call_val;

    }

    if(!fail){
        std::cout << "All cases passed for check 5." << std::endl;
    }

    // reset time to expiry for next test
    Time_to_expiry = 100;

    // //==================================================================================================================================================
    // //      check 6 - BS_call should be a convex function of K
    // //==================================================================================================================================================

    epsilon=0.001, K_max=140, K_min=60, K_step=(K_max-K_min)/n, K=K_min;
    fail = 0;

    //2nd derivative should be nonnegative everywhere
    double second_K_derivative;

    for(int i=0; i<n; i++){
        second_K_derivative = (C.BS_price(S, K+epsilon, sigma, r, Time_to_expiry, d)-2.0*C.BS_price(S, K, sigma, r, Time_to_expiry, d) + C.BS_price(S, K-epsilon, sigma, r, Time_to_expiry, d))/(epsilon*epsilon);

        if(second_K_derivative < 0){
            fail = 1;
        }

        K += K_step;
    }

    if(!fail){
        std::cout << "All cases passed for check 6." << std::endl;
    }

    // reset K for next test
    K = 100;

    // //==================================================================================================================================================
    // //      check 7 - The price of a call spread should approximate BS_digital_call
    // //==================================================================================================================================================

    // To approximate the digital call, we should buy 1/m calls struck at K-m and sell 1/m calls struck at K. As m->0, the approximation will get better.

    // acceptable difference between our portfolio and the digitial call that we are approximating
    epsilon = 0.01;

    Digital_Call DC = Digital_Call(S, K, r, d, Time_to_expiry, sigma);
    double m = 0.01;
    fail = 0;


    call_val = (C.BS_price(S,K-m, sigma, r, T, d) - C.BS_price(S, K, sigma, r, T, d))/m - DC.BS_price(S, K, sigma, r, T, d);

    if(call_val >= epsilon){
        fail = 1;
    }

    if(!fail){
        std::cout << "All cases passed for check 7." << std::endl;
    }

    // //==================================================================================================================================================
    // //      check 8 - BS_digital_call + BS_digital_put = zero_coupon_bond
    // //==================================================================================================================================================

    // acceptable difference between digitalcall + digital put and zero coupon bond.
    epsilon = 0.01;

    Digital_Put DP = Digital_Put(S, K, r, d, Time_to_expiry, sigma);
    Zero_Coupon_Bond Z = Zero_Coupon_Bond(S, r, d, Time_to_expiry);
    fail = 0;

    call_val = DC.BS_price(S, K, sigma, r, T, d) + DP.BS_price(S, K, sigma, r, T, d) - Z.BS_price(S, r, d, T);

    if(call_val >= epsilon){
        fail = 1;
    }

    if(!fail){
        std::cout << "All cases passed for check 8." << std::endl;
    }


    //==================================================================================================================================================
    //      Monte Carlo Pricing
    //==================================================================================================================================================

    /*
    We want to compute the prices for each option type using Monte Carlo. We calculate the prices for 2^n samples, with n=0,...,20.
    */

    RNG R = RNG();
    S=100, K=100, sigma=0.3, r=0.03, Expiry_time=100/365., d=0.0, m = 20;

    C = Call(S, K, r, d, Expiry_time, sigma);
    P = Put(S, K, r, d, Expiry_time, sigma);
    DC = Digital_Call(S, K, r, d, Expiry_time, sigma);
    DP = Digital_Put(S, K, r, d, Expiry_time, sigma);
  
    double put_val, DC_val, DP_val;
    std::vector<double> rands, call_prices, put_prices, DC_prices, DP_prices;

    for(int i=0; i<m; i++){
        rands = R.create_norm_vector(std::pow(2,i), 0, 1);
        call_val = C.MC_price(rands);
        call_prices.push_back(call_val);

        put_val = P.MC_price(rands);
        put_prices.push_back(put_val);

        DC_val = DC.MC_price(rands);
        DC_prices.push_back(DC_val); 

        DP_val = DP.MC_price(rands);
        DP_prices.push_back(DP_val);
    }

    output(call_prices, "call_vals");
    output(put_prices, "put_vals");
    output(DC_prices, "DC_vals");
    output(DP_prices, "DP_vals");

    //==================================================================================================================================================
    //      Investigations
    //==================================================================================================================================================

    /*
    1. How does the BS price of a call option vary as a function of volatility? What happens when volatility is zero, or volatility is very large?
    2. What about a digital call option?
    3. For various at-the-money call options, how does the price vary with volatility? Plot the ratio of price to volatility.
    4. For various put options plot the price and intrinsic value on the same graph. Find at least one example where the two graphs cross.
    */

    // Investigations 1+2
    S = 50, K=100, sigma = 0, n = 1000, sigma_max = 1, sigma_step = sigma_max/n;
    call_prices = {}, DC_prices = {};

    for(int i=0; i<n; i++){
        call_val = C.BS_price(S, K, sigma, r, Expiry_time, d);
        call_prices.push_back(call_val);

        DC_val = DC.BS_price(S, K, sigma, r, Expiry_time, d);
        DC_prices.push_back(DC_val);

        sigma += sigma_step;
    }

    output(call_prices, "call_volatility");
    output(DC_prices, "DC_volatility");


    // There's an issue here. My results are really big and jump around a lot. I'll need to look into it.

    // Investigation 3
    sigma=0.01;
    double price_vol_ratio = 0;
    std::vector<int> prices = {30, 50, 80, 100, 120, 150};
    std::vector<double> ratios;
    std::vector<double> sigma_values;

    for(int i=0; i<n; i++){
        sigma_values.push_back(sigma);
        sigma+=sigma_step;
    }

    for(int x : prices){
        call_prices = {}, ratios = {};
        sigma = 0.01;

        for(int i=0; i<n; i++){
            call_val = C.BS_price(x, x, sigma, r, Expiry_time, d);
            call_prices.push_back(call_val);
            ratios.push_back(call_val/sigma);

            sigma += sigma_step;
        }
        output(call_prices, "ATM_call_prices_" + std::to_string(x));
        output(ratios, "ATM_ratios_" + std::to_string(x));
    }
    output(sigma_values, "sigma");


    // Investigation 4
    S=0, sigma=0.7, n=20, K=10, r=0.1;
    double S_max = 20, S_step = S_max/n, IV;
    put_prices = {};
    std::vector<double> intrinsic_values;
    P = Put(S, K, r, d, Expiry_time, sigma);

    for(int i=0; i<n; i++){
        put_val = P.BS_price(S, K, sigma, r, Expiry_time, d);
        put_prices.push_back(put_val);

        IV = P.payout(S);
        intrinsic_values.push_back(IV);

        S += S_step;
    }

    output(put_prices, "put_prices");
    output(intrinsic_values, "intrinsic_values");
}