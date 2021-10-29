#include "database.h"
#include <iostream>

/*
Class for the SQLite database that stores any persistent data used by the system. Currently, it contains tables for stock and option prices and 
parameters.

Some of this was adapted from https://videlais.com/2018/12/14/c-with-sqlite3-part-5-encapsulating-database-objects/
*/


int Database::Callback(void *NotUsed, int argc, char **argv, char **azColName) {
    /*
    This can be called every time a query is executed. Right now, it just prints a row in the database to the console. 

    int argc: holds the number of results
    (array) azColName: holds each column returned
    (array) argv: holds each value
    */

    for(int i=0; i<argc; i++) {
        
        // Show column name and value
        std::cout << azColName[i] << ": " << argv[i] << std::endl;
    }

    std::cout << std::endl;

    return 0;
}

void Database::CheckDBErrors(){
    
    if(rc){
        // Show an error message
        std::cout << "DB Error: " << sqlite3_errmsg(db) << std::endl;
        CloseDB(); 
    }
}

Database::Database(){
    /*
    Save the result of opening the file
    */

    rc = sqlite3_open("holdings.db", &db);

    CheckDBErrors();
}

void Database::CreateTable(){
    /*
    This function creates the base tables that the database will use. Right now, we have one table for stocks and one for options. 

    As it is currently configured, a stock is uniquiely specified by its name. An option is uniquely specified by the combination of the name of the 
    underlying, the strike, the expiry data, and the type.
    */

    // SQL to create the tables
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

// void Database::InsertData(std::string table, std::string sym, double value){
//     /*
//     Insert a symbol and value pair into the chosen table
//     */

//     // The query that we will execute
//     std::string query = "INSERT INTO" + table + "('Symbol', 'Value') VALUES ('" + sym + "','" + std::to_string(value) + "');";

//     // Prepare the query
//     sqlite3_prepare(db, query.c_str(), query.length(), &stmt, NULL);

//     // Execute it
//     rc = sqlite3_step(stmt);

//     // Finialize the usage
//     sqlite3_finalize(stmt);    
// }

void Database::ShowTable(std::string table){
    /*
    Print the chosen table in the console.

    Known bug: when this is called in main(), nothing that comes after will be exceuted. Why?
    */

    std::string query = "SELECT * FROM '" + table + "';";

    // Run the SQL
    rc = sqlite3_exec(db, query.c_str(), Callback, 0, &zErrMsg);  
}

void Database::DeleteRow(std::string table, std::string sym) {
    /*
    Delete a row from the chosen table.

    TODO: Extend this to handle options. The symbol isn't enough to specify an option, so this won't work right. 
    */

    std::string query = "DELETE FROM '" + table + "' WHERE symbol = '" + sym + "';";

    // Prepare the query
    sqlite3_prepare(db, query.c_str(), query.length(), &stmt, NULL);

    // Run it
    rc = sqlite3_step(stmt);

    // Finialize the usage
    sqlite3_finalize(stmt);
}

void Database::UpdateData(std::string table, std::string sym, std::string parameter, double val){
    /*
    In the chosen table, change the specified parameter of security "sym" to val.

    TODO: Like DeleteRow, this needs to be updated to handle options. 
    */

    std::string query = "UPDATE " + table + "SET " + parameter + " = " + std::to_string(val) + "WHERE Symbol = '" + sym + "';";    
    
    // Prepare the query
    sqlite3_prepare(db, query.c_str(), query.length(), &stmt, NULL);

    // Run it
    rc = sqlite3_step(stmt);

    // Finialize the usage
    sqlite3_finalize(stmt);
}

double Database::GetStockParameter(std::string key, std::string parameter){
    /*
    This function takes a stock symbol as a key and returns the corresponding parameter in the Stocks table.

    Adapted from https://stackoverflow.com/questions/14437433/proper-use-of-callback-function-of-sqlite3-in-c
    */

    double value = 0;

    try{
        // sqlite3* d = DB.db;
        // sqlite3_stmt *stmt;
        
        // Create an SQL statement in a form that SQLite can understand
        int rc = sqlite3_prepare_v2(db, ("SELECT " + parameter + " FROM Stocks WHERE Symbol = '" + key + '\'').c_str(), -1, &stmt, NULL);

        if (rc != SQLITE_OK){
            throw std::string(sqlite3_errmsg(db));
        }

        // Excecute the statement
        rc = sqlite3_step(stmt);

        // Check for errors
        if (rc != SQLITE_ROW && rc != SQLITE_DONE) {
            std::string errmsg(sqlite3_errmsg(db));
            sqlite3_finalize(stmt);
            throw errmsg;
        }

        // Check whether the stock is actually in the table
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

double Database::GetOptionParameter(std::string symbol, double strike, double expiry_date, std::string type, std::string parameter){
    /*
    This function takes he name of the underlying, the strike, the expiry data, and the type as a key and returns the corresponding parameter in the 
    Options table.

    Adapted from https://stackoverflow.com/questions/14437433/proper-use-of-callback-function-of-sqlite3-in-c
    */

    double value = 0;

    try{
        // sqlite3* d = DB.db;
        // sqlite3_stmt *stmt;

        std::string query = "SELECT " + parameter + " FROM Options WHERE Symbol='" + symbol + "' AND Strike=" + std::to_string(strike) + " AND Expiry_date=" + std::to_string(expiry_date) + " AND Type='" + type + "';";
        
        // Create an SQL statement in a form that SQLite can understand
        int rc = sqlite3_prepare_v2(db, query.c_str(), -1, &stmt, NULL);

        if (rc != SQLITE_OK){
            throw std::string(sqlite3_errmsg(db));
        }

        // Excecute the statement
        rc = sqlite3_step(stmt);

        // Check for errors
        if (rc != SQLITE_ROW && rc != SQLITE_DONE) {
            std::string errmsg(sqlite3_errmsg(db));
            sqlite3_finalize(stmt);
            throw errmsg;
        }

        // Check whether the option is actually in the table
        if (rc == SQLITE_DONE) {
            sqlite3_finalize(stmt);
            throw std::string("Option not found");
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

void Database::CloseDB() {

    // Close the SQL connection
    sqlite3_close(db);
}

