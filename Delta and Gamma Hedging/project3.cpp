#include <iostream>
#include <random>
#include <string>
#include <unordered_map>
#include <vector>
#include <fstream>
#include <Eigen/Dense>
#include <math.h>
#include <sqlite3.h>
#include <string.h>

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

class Database{
    /*
    Class for the SQLite database that stores any persistent data used by the system. Currently, it contains tables for stock and option prices and 
    parameters.

    Some of this was adapted from https://videlais.com/2018/12/14/c-with-sqlite3-part-5-encapsulating-database-objects/
    */

    private:

    // Save any error messages
    char* zErrMsg;

    // Save the result of opening the file
    int rc;

    // Compiled SQLite Statement
    sqlite3_stmt* stmt;

    // This can be called every time a query is executed. Right now, it just prints a row in the database to the console. 
    static int Callback(void *NotUsed, int argc, char **argv, char **azColName) {

        // int argc: holds the number of results
        // (array) azColName: holds each column returned
        // (array) argv: holds each value

        for(int i=0; i<argc; i++) {
            
            // Show column name and value
            std::cout << azColName[i] << ": " << argv[i] << std::endl;
        }

        std::cout << std::endl;

        return 0;
    }

    void CheckDBErrors(){
        
        if(rc){
            // Show an error message
            std::cout << "DB Error: " << sqlite3_errmsg(db) << std::endl;
            CloseDB(); 
        }
    }

    public:
    // Pointer to SQLite connection
    sqlite3* db;

    Database(){
        // Save the result of opening the file
        rc = sqlite3_open("database.db", &db);

        CheckDBErrors();
    }

    void CreateTable(){

        // Save SQL to create a table
        std::string sql = "CREATE TABLE IF NOT EXISTS Options(                         "              
                          "    Symbol              TEXT     NOT NULL,                  "
                          "    Type                TEXT     NOT NULL,                  "
                          "    Number_owned        REAL     NOT NULL,                  "
                          "    Strike              REAL     NOT NULL,                  "
                          "    Expiry_date         INT      NOT NULL,                  "
                          "    Value               REAL,                               "
                          "    Delta               REAL,                               "
                          "    Gamma               REAL,                               "
                          "                                                            "
                          "    PRIMARY KEY(Symbol, Strike, Expiry_date, Type)          "
                          ");                                                          "
                          "                                                            "
                          "CREATE TABLE IF NOT EXISTS Stocks(                          "
                          "    Symbol           TEXT    PRIMARY KEY     NOT NULL,      "
                          "    Number_owned     REAL                    NOT NULL,      "
                          "    Value            REAL                                   "
                          ");                                                          "
                          "                                                            ";
        
        // Run the SQL
        rc = sqlite3_exec(db, sql.c_str(), NULL, 0, &zErrMsg);
    }

    // Insert a symbol and value pair into the chosen table
    void InsertData(std::string table, std::string sym, double value) {

        // The query that we will execute
        std::string query = "INSERT INTO" + table + "('Symbol', 'Value') VALUES ('" + sym + "','" + std::to_string(value) + "');";

        // Prepare the query
        sqlite3_prepare(db, query.c_str(), query.length(), &stmt, NULL);

        // Execute it
        rc = sqlite3_step(stmt);

        // Finialize the usage
        sqlite3_finalize(stmt);    
    }

    void ShowTable(std::string table) {

        std::string query = "SELECT * FROM '" + table + "';";

        // Run the SQL
        rc = sqlite3_exec(db, query.c_str(), Callback, 0, &zErrMsg);  
    }

    void DeleteRow(std::string table, std::string sym) {

        std::string query = "DELETE FROM '" + table + "' WHERE symbol = '" + sym + "';";

        // Prepare the query
        sqlite3_prepare(db, query.c_str(), query.length(), &stmt, NULL);

        // Run it
        rc = sqlite3_step(stmt);

        // Finialize the usage
        sqlite3_finalize(stmt);
    }

    // In the chosen table, change the specified parameter of security "sym" to val
    void UpdateData(std::string table, std::string sym, std::string parameter, double val){

        std::string query = "UPDATE " + table + "SET " + parameter + " = " + std::to_string(val) + "WHERE Symbol = '" + sym + "';";    
        
        // Prepare the query
        sqlite3_prepare(db, query.c_str(), query.length(), &stmt, NULL);

        // Run it
        rc = sqlite3_step(stmt);

        // Finialize the usage
        sqlite3_finalize(stmt);
    }

    double GetValue(std::string table, std::string key, std::string parameter){
        /*
        This function takes a key and finds the corresponding parameter in the given table.

        For now, this will only work for stocks, since the query includes "WHERE Symbol = x". Options will come later. I'm not exactly sure how I should 
        specify the key yet. For a stock, a row in the table is uniqely identified by the stock's symbol, but for an option, we also need the expiry date
        and strike price. I need a good way of incorporating all of this into one string, and then I need to figure out exactly how SQLite processes the
        information to turn that into a primary key.  

        Adapted from https://stackoverflow.com/questions/14437433/proper-use-of-callback-function-of-sqlite3-in-c
        */

       double value = 0;

        try{
            // sqlite3* d = DB.db;
            // sqlite3_stmt *stmt;
            
            // Create an SQL statement in a form that SQLite can understand
            int rc = sqlite3_prepare_v2(db, ("SELECT " + parameter + " FROM " + table + " WHERE Symbol = '" + key + '\'').c_str(), -1, &stmt, NULL);

            if (rc != SQLITE_OK){
                throw std::string(sqlite3_errmsg(db));
            }

            // Excecute the statement
            rc = sqlite3_step(stmt);

            if (rc != SQLITE_ROW && rc != SQLITE_DONE) {
                std::string errmsg(sqlite3_errmsg(db));
                sqlite3_finalize(stmt);
                throw errmsg;
            }

            if (rc == SQLITE_DONE) {
                sqlite3_finalize(stmt);
                throw std::string("Stock not found");
            }

            // Finally, store the value
            //this->value = sqlite3_column_double(stmt, 0);
            value = sqlite3_column_double(stmt, 0);

            // Behind-the-scenes stuff to prepare for another statement
            sqlite3_finalize(stmt);
        }

        catch(const std::string& ex){
            std::cout << ex << std::endl;
        }

        return value;
    }

    void CloseDB() {

        // Close the SQL connection
        sqlite3_close(db);
    }
};

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

class Security : public RNG{
    /*
    Class for a general security. A general class for each of options, stocks, futures, etc. will all inherit from this class. 

    "name" is the symbol/name/etc. of the given security. "value" gives the value of the security at the current time. "t" is the current time.
    */

    public:
    std::string name;
    double value;

    Security(std::string n, double initial_value): name(n), value(initial_value){}

    virtual void Buy(Database D, int amount) = 0;
    virtual void Sell(Database D, int amount) = 0;
};

class Stock : public Security{
    /*
    Class for a generic stock. 

    S_0 is the initial value at t=0, mu is the drift rate, sigma is the volatility, and d is the dividend rate.
    */

    public:
    double S_0, mu, sigma, d;

    Stock(std::string name, double S_0_, double mu_0, double sigma_0, double d_0): Security(name, S_0_), mu(mu_0), sigma(sigma_0), S_0(S_0_), d(d_0){}

    double SDE_price(double t){
        /*
        Return a lognormally distributed value, as determined by solving the SDE. 
        z should be a N(0,1) distributed random variable, as calculated by the RNG class
        */
        double z = gen_norm(0,1);
        return S_0*exp((mu-d)*t-0.5*sigma*sigma*t+sigma*z*sqrt(t));
    }

    void Buy(Database D, int amount){

        // Try and insert a new row for the stock we are buying. If a row already exists, update the amount we own instead. 
        std::string query = "INSERT INTO Stocks VALUES('" + name + "'," + std::to_string(amount) + ", NULL) "
                            "ON CONFLICT(Symbol) DO UPDATE SET Number_owned = Number_owned+" + std::to_string(amount) + ";";
    
        sqlite3* db;
        sqlite3_stmt* stmt;
        sqlite3_open("database.db", &db);
        int rc = sqlite3_open("database.db", &db);

        // CheckDBErrors();

        // Prepare the query
        sqlite3_prepare(db, query.c_str(), query.length(), &stmt, NULL);

        // Run it
        rc = sqlite3_step(stmt);

        // Finialize the usage
        sqlite3_finalize(stmt);

        // Now check to see if buying reduced our holdings to zero
        double current_amount = D.GetValue("Stocks", name, "Number_owned");

        // If so, erase the row from the database.
        if(current_amount == 0){
            D.DeleteRow("Stocks", name);
        } 
    }

    void Sell(Database D, int amount){

        // Try and insert a new row for the stock we are buying. If a row already exists, update the amount we own instead. 
        std::string query = "INSERT INTO Stocks VALUES('" + name + "'," + std::to_string(-amount) + ", NULL) "
                            "ON CONFLICT(Symbol) DO UPDATE SET Number_owned = Number_owned-" + std::to_string(amount) + ";";
    
        sqlite3* db;
        sqlite3_stmt* stmt;
        sqlite3_open("database.db", &db);
        int rc = sqlite3_open("database.db", &db);

        // CheckDBErrors();

        // Prepare the query
        sqlite3_prepare(db, query.c_str(), query.length(), &stmt, NULL);

        // Run it
        rc = sqlite3_step(stmt);

        // Finialize the usage
        sqlite3_finalize(stmt);

        // Now check to see if buying reduced our holdings to zero
        double current_amount = D.GetValue("Stocks", name, "Number_owned");

        // If so, erase the row from the database.
        if(current_amount == 0){
            D.DeleteRow("Stocks", name);
        } 
    }
};

class Option : public Stock{
    /*
    Class for a generic option. Specific option types will inherit from this.

    K is the strike, and T is the time until expiry.
    */

    public:
    double K, T;

    Option(std::string underlying_name, double S_0, double K_, double mu_0, double sigma_0, double d_0, double T_) 
    : K(K_), T(T_), Stock(underlying_name, S_0, mu_0, sigma_0, d_0){}

    // Equality operator used to determine hash collisions
    bool operator == (const Option& s) const{
        return name == s.name;
    }

    virtual double payout(double x){
        /*
        This is a generic payout function that will later be replaced by the payout from a child class for a specific option type.
        */
        return 0;
    }

    double d1(double S){
        // Calculate d1, as used in the Black-Scholes formula

        return (log(S/K)+(mu-d+0.5*sigma*sigma)*T)/(sigma*sqrt(T));   
    }

    double d2(double S){
        // Calculate d2, as used in the Black-Scholes formula

        return (log(S/K)+(mu-d-0.5*sigma*sigma)*T)/(sigma*sqrt(T));
    }

    double N(double x){
        /*
        Cumulative normal function
        */
        return 0.5*(1+erf(x/sqrt(2)));
    }

    double N_prime(double x){
        /*
        Derivative of N(x) - the PDF for the standard normal distribution
        */
        return exp(-0.5*x*x)/sqrt(2*M_PI);
    }
};

template <>
struct std::hash<Option>{
    size_t operator()(const Option& k) const{
        
        // Compute hash for an Option
        // http://stackoverflow.com/a/1646913/126995
        size_t res = 17;
        res = res * 31 + hash<std::string>()(k.name);
        return res;
    }
};

class Call : public Option{
    public:

    using Option::Option;

    double payout(double s){
        return std::max(s-K, 0.0);
    }

    double BS_price(double S){     
        return S*exp(-d*T)*N(d1(S))-K*exp(-mu*T)*N(d2(S));
    }

    double BS_Delta(double S){
        return exp(-d*T)*N(d1(S));
    }

    double BS_Gamma(double S){
        return exp(-d*T)*N_prime(d1(S))/(S*sigma*sqrt(T));
    }
};

template <class T>
class Strategy{
    /*
    A strategy consists of four elements: 
        1. A list of securities
        2. The amount of each security to buy/sell
        3. Whether to buy or sell each security (as a bool, this corresponds to 1 or 0)
        4. A condition for each security that needs to be met for the transaction to occur

    At each time step, we will evaluate the conditions of the strategy, and pass the buy/sell orders to the Portfolio class.
    */

    public:
    std::vector<T> securities;
    std::vector<double> amounts;
    std::vector<bool> action, conditions;

    Strategy(std::vector<T> s, std::vector<double> am, std::vector<bool> ac, std::vector<bool> c) : securities(s), amounts(am), action(ac), conditions(c){
        
        if(s.size != am.size() || s.size() != ac.size() || s.size() != c.size()){
            std::cout << "Missing one or more strategy parameters." << std::endl;
        }        
    }
};

// class Portfolio{

//     private:
//     double total_value = 0;
//     Database* db;

//     public:
//     // Keys are invidivual securities, values are the amounts of each. 
    
//     Portfolio(Database* d, std::vector<Strategy<Option>> option_strategies, std::vector<Strategy<Stock>> stock_strategies) : db(d){

//         // Set up the initial portfolio as specified by the input asset list. 
        
//     }

    // // Return the total value of the portfolio
    // double portfolio_value(){
    //     total_value = 0;

        
    // }

    // // Add an amount of a given asset to the portfolio.
    // void Buy(Security* s, double amount){
        
    // }

    // // Remove an amount of a given asset from the portfolio.
    // void Sell(Security* s, double amount){
        
    // }

// };

class Simulator{
    // Here is where we actually evolve in time. At each time step, we will check the conditions on buying and selling, then execute as needed. 

    double t, t_min, t_max, dt;
    Simulator(){}
};

int main(){

    double S = 100, K = 90, r = 0.03, d = 0.01, T = 60/365., sigma = 0.3, epsilon=0.01;
    Stock S1 = Stock("SPY", S, r, sigma, d);
    
    Database D = Database();
    D.CreateTable();
    
    S1.Buy(D, 5);

    // Something about the ShowTable command stops the program from running any lines that occur after it. WHY??
    //D.ShowTable("Stocks");

    // We can also use the get_value function to print the number of stocks owned, as well as to get other parameters if we need them
    std::cout << "Amount of SPY owned: " << D.GetValue("Stocks", "SPY", "Number_owned") << std::endl;
}