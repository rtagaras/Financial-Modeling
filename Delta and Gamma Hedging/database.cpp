#include "database.h"
#include <iostream>

/*
Class for the SQLite database that stores any persistent data used by the system. Currently, it contains tables for stock and option prices and 
parameters.

Some of this was adapted from https://videlais.com/2018/12/14/c-with-sqlite3-part-5-encapsulating-database-objects/
*/

// This can be called every time a query is executed. Right now, it just prints a row in the database to the console. 
int Database::Callback(void *NotUsed, int argc, char **argv, char **azColName) {

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

void Database::CheckDBErrors(){
    
    if(rc){
        // Show an error message
        std::cout << "DB Error: " << sqlite3_errmsg(db) << std::endl;
        CloseDB(); 
    }
}

Database::Database(){
    // Save the result of opening the file
    rc = sqlite3_open("holdings.db", &db);

    CheckDBErrors();
}

void Database::CreateTable(){

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
void Database::InsertData(std::string table, std::string sym, double value) {

    // The query that we will execute
    std::string query = "INSERT INTO" + table + "('Symbol', 'Value') VALUES ('" + sym + "','" + std::to_string(value) + "');";

    // Prepare the query
    sqlite3_prepare(db, query.c_str(), query.length(), &stmt, NULL);

    // Execute it
    rc = sqlite3_step(stmt);

    // Finialize the usage
    sqlite3_finalize(stmt);    
}

void Database::ShowTable(std::string table) {

    std::string query = "SELECT * FROM '" + table + "';";

    // Run the SQL
    rc = sqlite3_exec(db, query.c_str(), Callback, 0, &zErrMsg);  
}

void Database::DeleteRow(std::string table, std::string sym) {

    std::string query = "DELETE FROM '" + table + "' WHERE symbol = '" + sym + "';";

    // Prepare the query
    sqlite3_prepare(db, query.c_str(), query.length(), &stmt, NULL);

    // Run it
    rc = sqlite3_step(stmt);

    // Finialize the usage
    sqlite3_finalize(stmt);
}

// In the chosen table, change the specified parameter of security "sym" to val
void Database::UpdateData(std::string table, std::string sym, std::string parameter, double val){

    std::string query = "UPDATE " + table + "SET " + parameter + " = " + std::to_string(val) + "WHERE Symbol = '" + sym + "';";    
    
    // Prepare the query
    sqlite3_prepare(db, query.c_str(), query.length(), &stmt, NULL);

    // Run it
    rc = sqlite3_step(stmt);

    // Finialize the usage
    sqlite3_finalize(stmt);
}

double Database::GetValue(std::string table, std::string key, std::string parameter){
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

void Database::CloseDB() {

    // Close the SQL connection
    sqlite3_close(db);
}

