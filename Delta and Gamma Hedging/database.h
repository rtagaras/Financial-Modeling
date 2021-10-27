#include<string>
#include<sqlite3.h>

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

    // Insert a symbol and value pair into the chosen table
    void InsertData(std::string table, std::string sym, double value);

    // Print table to console
    void ShowTable(std::string table);

    // Delete the row corresponding to sym
    void DeleteRow(std::string table, std::string sym);

    // In the chosen table, change the specified parameter of security "sym" to val
    void UpdateData(std::string table, std::string sym, std::string parameter, double val);

    // Given a table and a primary key, return the value in the "parameter" column
    double GetValue(std::string table, std::string key, std::string parameter);

    void CloseDB(); 
};
