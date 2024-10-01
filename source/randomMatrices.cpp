#include <iostream>
#include <fstream>
#include <cmath>
#include <vector>
#include <random>
#include <string>
#include <Eigen/Dense> // library for algebra operations

std::vector<double> readData(const std::string& filePath);
void writeData(const std::string& filePath, Eigen::VectorXd vec, int i);
double semi_circle(double x, int N, double sigma);

int main(){

    std::cout << std::endl;
    // read inputs
    std::vector<double> inputs;
    inputs = readData("./inputs/input.txt");

    int repeat = int(inputs[0]);
    double mean = inputs[1];
    double beta = inputs[2];
    int dim = int(inputs[3]);
    std::cout << "Data load was complete." << std::endl;
    std::cout << std::endl << "#####################################################" << std::endl;

    // Making loop for creating rnd matrices
    std::random_device rd;
    std::mt19937 gen(rd());

    double sigma_diag = sqrt(1/(2*beta));
    double sigma_off = sqrt(1/(4*beta));

    for (int k = 0; k < repeat; k++){
        std::cout << "Start generating of matrix: " << k+1 << std::endl; 
        Eigen::MatrixXd matrix(dim, dim);

        for (int i = 0; i < dim; i++){
            for (int j = i; j < dim; j++){

                if (i == j){
                    std::normal_distribution<double> normalDist(mean, sigma_diag);
                    matrix(i, i) = normalDist(gen);
                }
                else{
                    // Hermitian real matrix: same values of diagonals
                    std::normal_distribution<double> normalDist(mean, sigma_off);
                    double num = normalDist(gen);
                    matrix(i, j) = num;
                    matrix(j, i) = num;
                }
            }
        }

        std::cout << "Start calculation of eigenvalues." << std::endl;

        // Get the eigenvalues as a vector
        Eigen::EigenSolver<Eigen::MatrixXd> solver(matrix);
        Eigen::VectorXd eigenvalues = solver.eigenvalues().real();

        // sort eigenvalues so I can cut the edges
        std::sort(eigenvalues.data(), eigenvalues.data() + eigenvalues.size());
        
        // Now reduce eigenvalues by 50 values on both sides (to stay in semi-circle law region)
        Eigen::VectorXd subvector = eigenvalues.segment(50, eigenvalues.size() - 100);

        std::cout << "Making data files." << std::endl;

        // write those eigenvalues into file
        writeData("./data/eigenvalues.txt", subvector, k);

        // Now semi-circle law
        Eigen::VectorXd deltaE(subvector.size() - 1);
        for(int i = 0; i < subvector.size() - 1; i++){
            deltaE(i) = fabs(subvector(i) - subvector(i + 1));
        }

        Eigen::VectorXd s(deltaE.size());
        for(int j = 0; j < deltaE.size(); j++){
            s(j) = deltaE(j)*semi_circle(subvector(j), dim, sigma_diag);
        }

        writeData("./data/s.txt", s, k);

        std::cout << "#####################################################" << std::endl;

    }


    return 0;
}

std::vector<double> readData(const std::string& filePath) {
    std::vector<double> numbers;

    std::ifstream inputFile(filePath);

    if (!inputFile.is_open()) {
        std::cout << "Error opening file: " << filePath << std::endl;
        return numbers;
    }

    std::string line;
    while (std::getline(inputFile, line)) {
        // Ignore everything after '#' and remove whitespaces
        size_t comment = line.find('#');
        if (comment != std::string::npos) {
            line.erase(comment);
        }

        std::stringstream ss(line);

        double number;
        while (ss >> number) {
            numbers.push_back(number);
        }
    }

    inputFile.close();
    return numbers;
}

void writeData(const std::string& filePath, Eigen::VectorXd vec, int i){

    std::ofstream outputFile;

    if (i == 0){
        outputFile.open(filePath); // when is called first time it overwrites old data
    }
    else{
        outputFile.open(filePath, std::ios_base::app); // append data from all repeting 
    }

    if (!outputFile.is_open()) {
        std::cout << "Error opening file for writing" << std::endl;
    }

    // Write the eigenvalues to the file
    for (int i = 0; i < vec.size(); ++i) {
        outputFile << vec(i) << std::endl;
    }

    // Close the file
    outputFile.close();
}

double semi_circle(double x, int N, double sigma){
    double bl = (1/(2*M_PI*sigma*sigma))*sqrt(4*N - x*x);
    return bl;
}





