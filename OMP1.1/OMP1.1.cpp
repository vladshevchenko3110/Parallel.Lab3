#include <iostream>
#include <omp.h>
#include <time.h>
#include <chrono>
#include <fstream>
#include <sstream>
#include <vector>

using namespace std;

class Timer {
private: 
    using clock_t = chrono::high_resolution_clock;
    using second_t = chrono::duration<double, ratio<1,1000> >;
    chrono::time_point<clock_t> m_beg;
public:
    Timer() : m_beg(clock_t::now()) {}
    void reset() { m_beg = clock_t::now(); }
    double elapsed() const {
        return chrono::duration_cast<second_t>(clock_t::now() - m_beg).count();
    }
};

void mxcComparison(double coefficient, int loadmatrixsize) {

    cout << "======== Comparison with a matrix of " << loadmatrixsize * loadmatrixsize << " elements ========" << endl;
    vector<vector<double>> loadedMatrix(loadmatrixsize, vector<double>(loadmatrixsize));
    ifstream inFile("matrix.bin", ios::binary);
    for (int i = 0; i < loadmatrixsize; ++i) {
        inFile.read(reinterpret_cast<char*>(loadedMatrix[i].data()), loadmatrixsize * sizeof(double));
    }
    inFile.close();

    double sum = 0;
    Timer t1;
    for (int i = 0; i < loadmatrixsize; ++i) {
        for (int j = 0; j < loadmatrixsize; ++j) {
            loadedMatrix[i][j] *= coefficient;
            sum += loadedMatrix[i][j];
        }
    }
    double time1 = t1.elapsed();
    cout << "Sequential processing time : " << time1 << " ms" << '\n';
    cout << "Total sum: " << sum << endl;
    
    ofstream newOutFile("serial_new_matrix_" + to_string(loadmatrixsize) + ".bin", ios::binary);
    ofstream newOutFileResult("serial_new_result_" + to_string(loadmatrixsize) + ".bin", ios::binary);
    for (int i = 0; i < loadmatrixsize; ++i) {
        newOutFile.write(reinterpret_cast<char*>(loadedMatrix[i].data()), loadmatrixsize * sizeof(double));
    }    
    newOutFile.close();

    newOutFileResult.write(reinterpret_cast<char*>(&sum), sizeof(double));
    newOutFileResult.close();


    vector<vector<double>> omploadedMatrix(loadmatrixsize, vector<double>(loadmatrixsize));
    ifstream inFileOMP("matrix.bin", ios::binary);
    for (int i = 0; i < loadmatrixsize; ++i) {
        inFileOMP.read(reinterpret_cast<char*>(omploadedMatrix[i].data()), loadmatrixsize * sizeof(double));
    }
    inFileOMP.close();

    double sum1 = 0;
    int chunk = 100;
    Timer t2;
#pragma omp parallel for schedule(dynamic,chunk) num_threads(10) reduction(+:sum1) 
    for (int i = 0; i < loadmatrixsize; ++i) {
        for (int j = 0; j < loadmatrixsize; ++j) {
            omploadedMatrix[i][j] *= coefficient;
            sum1 += omploadedMatrix[i][j];
        }
    }
    double time2 = t2.elapsed();
    cout << "Parallel processing time : " << time2 << " ms" << '\n';
    cout << "Total sum: " << sum1 << endl;
    cout << "Acceleration of parallel block processing: " << time1 / time2 << endl;

    ofstream ompnewOutFile("parallel_new_matrix_" + to_string(loadmatrixsize) + ".bin", ios::binary);
    ofstream ompnewOutFileResult("parallel_new_result_" + to_string(loadmatrixsize) + ".bin", ios::binary);
    for (int i = 0; i < loadmatrixsize; ++i) {
        ompnewOutFile.write(reinterpret_cast<char*>(omploadedMatrix[i].data()), loadmatrixsize * sizeof(double));
    }
    ompnewOutFile.close();

    ompnewOutFileResult.write(reinterpret_cast<char*>(&sum1), sizeof(double));
    ompnewOutFileResult.close();
    cout << endl;
}


int main(int argc, char* argv[])
{
    srand((unsigned int)time(0));
    double coefficient = stod(argv[1]);

    /*const int matrixSize = 10000;*/

    
    //vector<vector<double>> matrix(matrixSize, std::vector<double>(matrixSize));
    //for (int i = 0; i < matrixSize; ++i) {
    //    for (int j = 0; j < matrixSize; ++j) {
    //        matrix[i][j] = rand() % 100; 
    //    }
    //}
    
    //ofstream outFile("matrix.bin", ios::binary);
    //for (int i = 0; i < matrixSize; ++i) {
    //    outFile.write(reinterpret_cast<char*>(matrix[i].data()), matrixSize * sizeof(double));
    //}
    //outFile.close();
    
    mxcComparison(coefficient, stoi(argv[2]));
    mxcComparison(coefficient, stoi(argv[3]));
    mxcComparison(coefficient, stoi(argv[4]));
    mxcComparison(coefficient, stoi(argv[5]));
    mxcComparison(coefficient, stoi(argv[6])); 
    mxcComparison(coefficient, stoi(argv[7]));
    mxcComparison(coefficient, stoi(argv[8]));
    mxcComparison(coefficient, stoi(argv[9]));
    mxcComparison(coefficient, stoi(argv[10]));
    mxcComparison(coefficient, stoi(argv[11]));

}

