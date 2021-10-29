#include <random>
#include <vector>
#include "database.h"

class RNG{

    private:
    std::mt19937 gen;

    public:

    RNG();
    double gen_norm(double mean, double variance);
    double gen_uniform(double min, double max);
    std::vector<double> create_norm_vector(int n, double mean, double variance);
};

class Security : public RNG{
    /*
    Class for a general security. A general class for each of options, stocks, futures, etc. will all inherit from this class. 

    "name" is the symbol/name/etc. of the given security. "value" gives the value of the security at the current time. "t" is the current time.
    */

    public:
    std::string name;
    double value, t;

    Security(std::string n, double initial_value);
    virtual void Buy(Database D, int amount) = 0;
};

class Stock : public Security{
    /*
    Class for a generic stock. 

    S_0 is the initial value at t=0, mu is the drift rate, sigma is the volatility, and d is the dividend rate.
    */

    public:
    double S_0, mu, sigma, d;

    Stock(std::string name, double S_0_, double mu_0, double sigma_0, double d_0);
    double SDE_price(double t);

    void Buy(Database D, int amount);
};

class Option : public Stock{
    /*
    Class for a generic option. Specific option types will inherit from this.

    K is the strike and T is the time until expiry, which occurs at Expiry_time.
    */

    public:
    double K, T, Expiry_time;

    Option(std::string underlying_name, double S_0, double K_, double mu_0, double sigma_0, double d_0, double exp);

    virtual double payout(double x);
    double d1(double S);
    double d2(double S);
    double N(double x);
    double N_prime(double x);
};

class Call : public Option{
    public:

    using Option::Option;

    double payout(double s);
    double BS_price(double S);
    double BS_Delta(double S);
    double BS_Gamma(double S);

    void Buy(Database D, int amount);
};
