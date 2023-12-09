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
    using second_t = chrono::duration<double, ratio<1, 1000> >;
    chrono::time_point<clock_t> m_beg;
public:
    Timer() : m_beg(clock_t::now()) {}
    void reset() { m_beg = clock_t::now(); }
    double elapsed() const {
        return chrono::duration_cast<second_t>(clock_t::now() - m_beg).count();
    }
};

void mxvComparison(int loadmatrixsize) {

    vector<double> result(loadmatrixsize);
    cout << "======== Comparison with a matrix of " << loadmatrixsize * loadmatrixsize << " elements ========" << endl;
    vector<vector<double>> loadedMatrix(loadmatrixsize, vector<double>(loadmatrixsize));
    ifstream MinFile("matrix.bin", ios::binary);
    for (int i = 0; i < loadmatrixsize; ++i) {
        MinFile.read(reinterpret_cast<char*>(loadedMatrix[i].data()), loadmatrixsize * sizeof(double));
    }
    MinFile.close();

    vector<double> loadedVector(loadmatrixsize);
    ifstream VinFile("vector.bin", ios::binary);    
    VinFile.read(reinterpret_cast<char*>(loadedVector.data()), loadmatrixsize * sizeof(double));    
    VinFile.close();
    
    
    Timer t1;
    for (int i = 0; i < loadmatrixsize; ++i) {
        for (int j = 0; j < loadmatrixsize; ++j) {
            result[i] += loadedMatrix[i][j] * loadedVector[j];
        }
    }
    double time1 = t1.elapsed();
    cout << "Sequential processing time : " << time1 << " ms" << '\n';
    
    
    
    ofstream newOutFileResult("serial_new_result_" + to_string(loadmatrixsize) + ".bin", ios::binary);    

    newOutFileResult.write(reinterpret_cast<char*>(result.data()), loadmatrixsize * sizeof(double));
    newOutFileResult.close();

    vector<double> ompresult(loadmatrixsize);
    vector<vector<double>> omploadedMatrix(loadmatrixsize, vector<double>(loadmatrixsize));
    ifstream inFileOMP("matrix.bin", ios::binary);
    for (int i = 0; i < loadmatrixsize; ++i) {
        inFileOMP.read(reinterpret_cast<char*>(omploadedMatrix[i].data()), loadmatrixsize * sizeof(double));
    }
    inFileOMP.close();

    vector<double> omploadedVector(loadmatrixsize);
    ifstream ompVinFile("vector.bin", ios::binary);
    ompVinFile.read(reinterpret_cast<char*>(omploadedVector.data()), loadmatrixsize * sizeof(double));
    ompVinFile.close();

    int chunk = 100;
    Timer t2;    
#pragma omp parallel for schedule(dynamic,chunk) num_threads(10)           
    for (int i = 0; i < loadmatrixsize; ++i) {
        for (int j = 0; j < loadmatrixsize; ++j) {
            ompresult[i] += loadedMatrix[i][j] * loadedVector[j];
        }
    }   
    
    double time2 = t2.elapsed();
    cout << "Parallel processing time : " << time2 << " ms" << '\n';    
    cout << "Acceleration of parallel block processing: " << time1 / time2 << endl;

    
    ofstream ompnewOutFileResult("parallel_new_result_" + to_string(loadmatrixsize) + ".bin", ios::binary);  

    ompnewOutFileResult.write(reinterpret_cast<char*>(ompresult.data()), loadmatrixsize * sizeof(double));
    ompnewOutFileResult.close();
    cout << endl;
}


int main(int argc, char* argv[])
{
    /*srand((unsigned int)time(0));


    const int matrixSize = 10000;


    vector<vector<double>> matrix(matrixSize, vector<double>(matrixSize));
    for (int i = 0; i < matrixSize; ++i) {
        for (int j = 0; j < matrixSize; ++j) {
            matrix[i][j] = rand() % 10000;
        }
    }

    ofstream MoutFile("matrix.bin", ios::binary);
    for (int i = 0; i < matrixSize; ++i) {
        MoutFile.write(reinterpret_cast<char*>(matrix[i].data()), matrixSize * sizeof(double));
    }
    MoutFile.close();

    vector<double> vector(matrixSize);
    for (int i = 0; i < matrixSize; ++i) {       
        vector[i] = rand() % 10000;        
    }

    ofstream VoutFile("vector.bin", ios::binary);
    VoutFile.write(reinterpret_cast<char*>(vector.data()), matrixSize * sizeof(double));    
    VoutFile.close();*/
   
    mxvComparison(stoi(argv[1]));
    mxvComparison(stoi(argv[2]));
    mxvComparison(stoi(argv[3]));
    mxvComparison(stoi(argv[4]));
    mxvComparison(stoi(argv[5]));
    mxvComparison(stoi(argv[6]));
    mxvComparison(stoi(argv[7]));
    mxvComparison(stoi(argv[8]));
    mxvComparison(stoi(argv[9]));
    mxvComparison(stoi(argv[10]));
    

}