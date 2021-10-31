#include "securities.h"
#include <optional>

//======================================================================================================================================================
//  Implementation of RNG class
//======================================================================================================================================================

/*
Class that holds all random number generation functions. Anything that needs RNG will inherit from this. 
*/
RNG::RNG() : gen((std::random_device())()) {}

double RNG::gen_norm(double mean, double variance){
    /*
    Sample a normal distribution with given mean and variance
    */
    std::normal_distribution<double> d(mean, sqrt(variance)); 
    return d(gen);
}

double RNG::gen_uniform(double min, double max){
    /*
    sample a uniform distribution
    */

    std::uniform_real_distribution<double> d(min, max);
    return d(gen);
}

std::vector<double> RNG::create_norm_vector(int n, double mean, double variance){
    /*
    Creates a vector of N(mean, variance) random variables
    */
    std::vector<double> data;

    for(int i=0; i<n; i++){
        data.push_back(gen_norm(mean, variance));
    }

    return data;
}

//======================================================================================================================================================
//  Implementation of Security class
//======================================================================================================================================================

/*
Class for a general security. A general class for each of options, stocks, futures, etc. will all inherit from this class. 

"name" is the symbol/name/etc. of the given security. "value" gives the value of the security at the current time. "t" is the current time.
*/

Security::Security(std::string n, double initial_value): name(n), value(initial_value){}

//======================================================================================================================================================
//  Implementation of Stock class
//======================================================================================================================================================

/*
Class for a generic stock. 

S_0 is the initial value at t=0, mu is the drift rate, sigma is the volatility, and d is the dividend rate.
*/

Stock::Stock(std::string name, double S_0_, double mu_0, double sigma_0, double d_0): Security(name, S_0_), mu(mu_0), sigma(sigma_0), S_0(S_0_), d(d_0){}

double Stock::SDE_price(double t){
    /*
    Return a lognormally distributed value, as determined by solving the SDE. 
    z should be a N(0,1) distributed random variable, as calculated by the RNG class
    */
    double z = gen_norm(0,1);
    return S_0*exp((mu-d)*t-0.5*sigma*sigma*t+sigma*z*sqrt(t));
}

void Stock::Buy(Database D, int amount){
    /*
    This function is used to buy or sell a stock. "D" specifies the database that we want to update when conducting the transaction and "amount" is the
    amount we want to buy or sell (positive for buying, negative for selling).

    Closing positions is automatically handled. 
    */

    // Try and insert a new row for the stock we are buying. If a row already exists, update the amount we own instead. 
    std::string query = "INSERT INTO Stocks VALUES('" + name + "'," + std::to_string(amount) + ", NULL, NULL) "
                        "ON CONFLICT(Symbol) DO UPDATE SET Number_owned = Number_owned+" + std::to_string(amount) + ";";

    sqlite3* db;
    sqlite3_stmt* stmt;
    sqlite3_open("holdings.db", &db);
    int rc = sqlite3_open("holdings.db", &db);

    D.CheckDBErrors();

    // Prepare the query
    sqlite3_prepare(db, query.c_str(), query.length(), &stmt, NULL);

    // Run it
    rc = sqlite3_step(stmt);

    // Finialize the usage
    sqlite3_finalize(stmt);

    // Check how much we own after buying/selling. This is needed to update the total value of the security and to see if buying/selling reduced our 
    //holdings to zero.
    double current_amount = D.GetStockParameter(name, "Number_owned");

    // Update the total value of whatever we just bought
    D.UpdateData("Stocks", "Total_value", value*current_amount, name);

    // If we don't own anything, erase the row from the database.
    if(current_amount == 0){
        D.DeleteRow("Stocks", name, std::nullopt, std::nullopt, std::nullopt);
    }
}

//======================================================================================================================================================
//  Implementation of Option class
//======================================================================================================================================================

/*
Class for a generic option. Specific option types will inherit from this.

K is the strike, and T is the time until expiry.
*/

Option::Option(std::string underlying_name, double S_0, double K_, double mu_0, double sigma_0, double d_0, double exp) : K(K_), Expiry_time(exp), Stock(underlying_name, S_0, mu_0, sigma_0, d_0){
    T = Expiry_time - t;
}

double Option::payout(double x){
    /*
    This is a generic payout function that will later be replaced by the payout from a child class for a specific option type.
    */
    return 0;
}

double Option::d1(double S){
    // Calculate d1, as used in the Black-Scholes formula

    return (log(S/K)+(mu-d+0.5*sigma*sigma)*T)/(sigma*sqrt(T));   
}

double Option::d2(double S){
    // Calculate d2, as used in the Black-Scholes formula

    return (log(S/K)+(mu-d-0.5*sigma*sigma)*T)/(sigma*sqrt(T));
}

double Option::N(double x){
    /*
    Cumulative normal function
    */
    return 0.5*(1+erf(x/sqrt(2)));
}

double Option::N_prime(double x){
    /*
    Derivative of N(x) - the PDF for the standard normal distribution
    */
    return exp(-0.5*x*x)/sqrt(2*M_PI);
}

//======================================================================================================================================================
//  Implementation of Call class
//======================================================================================================================================================

double Call::payout(double s){
    return std::max(s-K, 0.0);
}

double Call::BS_price(double S){     
    return S*exp(-d*T)*N(d1(S))-K*exp(-mu*T)*N(d2(S));
}

double Call::BS_Delta(double S){
    return exp(-d*T)*N(d1(S));
}

double Call::BS_Gamma(double S){
    return exp(-d*T)*N_prime(d1(S))/(S*sigma*sqrt(T));
}

void Call::Buy(Database D, int amount){
    /*
    This function is used to buy or sell an option. "D" specifies the database that we want to update when conducting the transaction and "amount" is the
    amount we want to buy or sell (positive for buying, negative for selling).

    Closing positions is automatically handled. 
    */

    // Try and insert a new row for the option we are buying. If a row already exists, update the amount we own instead. 
    std::string query = "INSERT INTO Options VALUES  ('" + name + "', 'CALL', " + std::to_string(amount) + ", " + std::to_string(K) + ", " + std::to_string(Expiry_time) + ", NULL, NULL)"
                        "ON CONFLICT(Symbol, Strike, Expiry_date, Type) DO UPDATE SET Number_owned = Number_owned+" + std::to_string(amount) + ";";

    sqlite3* db;
    sqlite3_stmt* stmt;
    sqlite3_open("holdings.db", &db);
    int rc = sqlite3_open("holdings.db", &db);

    D.CheckDBErrors();

    // Prepare the query
    sqlite3_prepare(db, query.c_str(), query.length(), &stmt, NULL);

    // Run it
    rc = sqlite3_step(stmt);

    // Finialize the usage
    sqlite3_finalize(stmt);

    // Check how much we own after buying/selling. This is needed to update the total value of the security and to see if buying/selling reduced our 
    //holdings to zero.
    double current_amount = D.GetOptionParameter(name, K, Expiry_time, "CALL", "Number_owned");

    // Update the total value of whatever we just bought
    D.UpdateData("Options", "Total_value", value*current_amount, name, K, Expiry_time, "CALL");

    // If we don't own anything, erase the row from the database.
    if(current_amount == 0){
        D.DeleteRow("Options", name, K, Expiry_time, "CALL");
    }
}