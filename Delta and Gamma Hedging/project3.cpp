#include <iostream>
#include <fstream>
#include <Eigen/Dense>
#include "securities.h"

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

// template <class T>
// class Strategy{
//     /*
//     A strategy consists of four elements: 
//         1. A list of securities
//         2. The amount of each security to buy/sell
//         3. Whether to buy or sell each security (as a bool, this corresponds to 1 or 0)
//         4. A condition for each security that needs to be met for the transaction to occur

//     At each time step, we will evaluate the conditions of the strategy, and pass the buy/sell orders to the Portfolio class.
//     */

//     public:
//     std::vector<T> securities;
//     std::vector<double> amounts;
//     std::vector<bool> action, conditions;

//     Strategy(std::vector<T> s, std::vector<double> am, std::vector<bool> ac, std::vector<bool> c) : securities(s), amounts(am), action(ac), conditions(c){
        
//         if(s.size != am.size() || s.size() != ac.size() || s.size() != c.size()){
//             std::cout << "Missing one or more strategy parameters." << std::endl;
//         }        
//     }
// };

// class Simulator{
//     // Here is where we actually evolve in time. At each time step, we will check the conditions on buying and selling, then execute as needed. 

//     double t, t_min, t_max, dt;
//     int n;
//     std::vector<Stock> stocks;
//     std::vector<Option> options;
//     Eigen::MatrixXd stock_prices, option_prices;
    
//     Simulator(std::vector<Stock> stocks_, std::vector<Option> options_, double t_min_, double t_max_, double dt_): stocks(stocks_), options(options_), t_min(t_min_), t_max(t_max_), dt(dt_), t(t_min_){
//         n = (t_max-t_min)/dt;
//         stock_prices = MatrixXd::Zero(n,stocks.size());
//         option_prices = MatrixXd::Zero(n, options.size());
//     }


//     // At each time step, do whatever calculations we need to do for each stock and option
//     // We also need to check at each time step whether the strategies trigger or not, then execute them if they do.
//     void Simulate(){
//         for(int i=0; i<n; i++){
//             for(Stock s : stocks){
                
//             }
//             t += dt;
//         }
//     }
// };

int main(){

    // Security parameters
    double S_0 = 100, K = 90, r = 0.03, d = 0.01, T = 60/365., sigma = 0.3, epsilon=0.01, Expiry_time=0;
    Call C = Call("SPY", S_0, K, r, sigma, d, Expiry_time);
    Stock S = Stock("AAPL", S_0, r, sigma, d);

    // Manually set these just to check that GetTotalValue works
    C.value = 2;
    S.value = 3;

    // Create a test portfolio
    Database D = Database();
    D.CreateTable();
    C.Buy(D, 5);
    S.Buy(D,1);

    // We can use the get_value function to print the number of stocks owned, as well as to get other parameters if we need them
    std::cout << "Amount of SPY Call owned: " << D.GetOptionParameter("SPY", K, Expiry_time, "CALL", "Number_owned") << std::endl;
    std::cout << "Amount of AAPL owned: " << D.GetStockParameter("AAPL", "Number_owned") << std::endl;
    std::cout << "Total portfolio value is: " << D.GetTotalValue() << std::endl;
}