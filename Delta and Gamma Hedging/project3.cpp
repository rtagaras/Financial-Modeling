#include <iostream>
#include <fstream>
#include "securities.h"

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