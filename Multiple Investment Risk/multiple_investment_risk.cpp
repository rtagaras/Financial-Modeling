#include <iostream>
#include <random>
#include <fstream>
#include <string>
#include <vector>
#include <Eigen/Dense>
 
using Eigen::MatrixXd;

//output a vector
template<class T>
void output(std::vector<T> v, std::string filename){
    std::ofstream ofs;
    ofs.open(filename + ".txt");
    
    for(int i=0; i<v.size(); i++){
        ofs << v[i] << '\n';
    }
    
    ofs.close();
}

std::random_device rd;
std::mt19937 gen(rd());
double gen_norm(double mean, double variance){
    std::normal_distribution<double> d(mean, sqrt(variance)); 
    return d(gen);
}

// geometric random walk with drift mu and volatility sigma
// T_max gives number of days to calculate path for
// s_0 is initial stock price at day zero
struct GRW{

    double mu, sigma, T_max, dt, s_0, s, z = 0.0;
    int n, steps_per_day = 0;

    GRW(double s_0_, double mu_, double sigma_, double T_max_, double dt_){
        s_0 = s_0_;
        mu = mu_;
        sigma = sigma_;
        T_max = T_max_; 
        dt = dt_;
        n  = T_max/dt;
        steps_per_day = 1/dt;
    }

    // calculate price as a function of time, measured in days
    std::vector<double> path(){
        std::vector<double> data;

        s = s_0;
        
        for(int i=0; i<n; i++){

            z = gen_norm(0.0, 1.0);
            s = s*(1.0 + mu*dt+sigma*z*sqrt(dt));

            data.push_back(s);
        }

        return data;
    }
};

// Given the correlations between pairs of variables and the individual standard deviations of each variable, return the correlation matrix. 
// "correlations" should be an upper triangular nxn matrix with zeros on the diagonal where the (i,j)th entry is the correlation coefficient 
// for variables i and j.
// The ith element of "std_deviations" should be the standard deviation for the ith variable.
Eigen::MatrixXd construct_correlation_matrix(std::vector<double> std_deviations, Eigen::MatrixXd correlations){

    int x = correlations.cols();
    int y = correlations.rows();
    Eigen::MatrixXd C(x, y);

    for(int i=0; i<x; i++){
        C(i,i) = std::pow(std_deviations[i], 2);
    }

    for(int i=0; i<y; i++){
        for(int j=i+1; j<x; j++){
            C(i,j) = correlations(i,j) * std_deviations[i] * std_deviations[j];
            C(j,i) = C(i,j);
        }
    }

    return C;
}

// Given a correlation matrix C, generate a vector of correlated random variables. 
// C can be generated using "construct_correlation_matrix" from correlation coefficients and individual standard deviations.
Eigen::VectorXd correlated_samples(Eigen::MatrixXd C){

    // A Hermitian, positive definite matrix M can be decomposed as M = LL^T, where L is lower triangular.
    // This gives the matrix L, which is used in generating correlated samples.
    Eigen::MatrixXd L = C.llt().matrixL();


    // is this broken, or do I just have a bad matrix as input?
    std::cout << L << std::endl << L*L.transpose() << std::endl;

    // Vector of uncorrelated N(0,1) random variables
    Eigen::VectorXd Z(C.cols());

    // Vector of random variables with correlations and standard deviations defined by choice of C. 
    Eigen::VectorXd X;

    for(int i=0; i<C.cols(); i++){
        Z(i) = gen_norm(0,1);
    }

    X = L*Z;

    return X;
}
 
int main(){

    Eigen::MatrixXd mat = Eigen::MatrixXd::Zero(3,3);
    mat(0,1) = 1.0;
    mat(0,2) = 2.0;
    mat(1,2) = 3.0;

    std::vector<double> devs = {4.0, 5.0, 6.0};

    Eigen::MatrixXd C = construct_correlation_matrix(devs, mat);

    Eigen::VectorXd X = correlated_samples(C);

}