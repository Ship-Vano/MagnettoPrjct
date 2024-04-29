//
// Функции для записи результатов вычислений - векторов и матриц - в файл
// @author: SHIP-VAN0
// @author:
//

#ifndef MAGNETTOPRJCT_FILEIO_H
#define MAGNETTOPRJCT_FILEIO_H

#include<fstream>
#include<istream>
#include<vector>
#include<string>
#include<sstream>
#include<iostream>
#include <iomanip>
using namespace std;

//TODO: вынести шаблоны функций в .cpp файл
//запись вектора в файл
template <typename DT>
void writeVectorToFile(ofstream& file, vector<DT> v) {
    for (int i = 0; i < v.size(); ++i)
        file << setprecision(16) << v[i] << " ";
    file << " " << endl;
}

template <typename DT>
void writeVectorToFile(fstream& file, vector<DT> v) {
    for (int i = 0; i < v.size(); ++i)
        file << setprecision(16) << v[i] << " ";
    file << " " << endl;
}
template <typename DT>
void writeVectorToFile(ofstream& file, DT v_0, vector<DT> v)
{
    file << v_0 << " ";
    for (int i = 0; i < v.size(); ++i)
        file << setprecision(16) << v[i] << " ";
    file << " " << endl;
}
template <typename DT>
void writeVectorToFile(fstream& file, DT v_0, vector<DT> v)
{
    file << v_0 << " ";
    for (int i = 0; i < v.size(); ++i)
        file << setprecision(16) <<  v[i] << " ";
    file << " " << endl;
}

template<typename DT>
void write_data_to_file(string filepath, vector<vector<DT>> data)
{
    ofstream output_data;
    output_data.open(filepath);
    int n = data.size();
    int m = data[0].size();
    cout << m << endl;
    for (int i = 0; i < n; ++i)
    {
        for (int j = 0; j < m; ++j)
        {
            output_data << data[i][j] << "  ";
        }
        output_data << endl;
    }
    output_data.close();
}


#endif //MAGNETTOPRJCT_FILEIO_H