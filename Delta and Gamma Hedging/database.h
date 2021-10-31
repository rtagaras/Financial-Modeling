#include <string>
#include <sqlite3.h>
#include <vector>
#include <optional>

class Database{
    /*
    Class for the SQLite database that stores any persistent data used by the system. Currently, it contains tables for stock and option prices and 
    parameters.
    */

    private:

    // Save any error messages
    char* zErrMsg;

    // Save the result of opening the file
    int rc;

    // Compiled SQLite Statement
    sqlite3_stmt* stmt;

    // This can be called every time a query is executed.
    static int Callback(void *NotUsed, int argc, char **argv, char **azColName);

    public:
    // Pointer to SQLite connection
    sqlite3* db;

    void CheckDBErrors();

    Database();

    void CreateTable();

    // // Insert a symbol and value pair into the chosen table
    // void InsertData(std::string table, std::string sym, double value);

    // Print table to console
    void ShowTable(std::string table);

    // Delete the row corresponding to sym
    void DeleteRow(std::string table, std::string sym, std::optional<double> strike, std::optional<double> expiry_date, std::optional<std::string> type);

    // In the chosen table, change the specified parameter of security "sym" to val
    void UpdateData(std::string table, std::string parameter, double val, std::string sym, std::optional<double> strike=std::nullopt, std::optional<double> expiry_date=std::nullopt, std::optional<std::string> type=std::nullopt);

    // Given a symbol for a stock, return the value in the "parameter" column
    double GetStockParameter(std::string key, std::string parameter);

    // Given a symbol, strike, expiry date, and type for an option, return the value in the "parameter" column
    double GetOptionParameter(std::string symbol, double strike, double expiry_date, std::string type, std::string parameter);

    // Return the names of each table in the database
    std::vector<std::string> GetTableNames();

    // Return the total value of all the securities in the database
    double GetTotalValue();

    void CloseDB(); 
};
