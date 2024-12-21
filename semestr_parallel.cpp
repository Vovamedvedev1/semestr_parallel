#include <iostream>
#include <iomanip>
#include <windows.h>
#include <sstream>
#include <fstream>
#include <string>
#include <time.h>
#include <mpi.h>
#include <vector>
#include <map>
using namespace std;

void permutation(vector<int>& vc)
{
    int i, j;
    for (i = vc.size() - 2; i >= 1; i--)
    {
        if (vc[i] < vc[i + 1])
            break;
    }
    for (j = vc.size() - 2; j > i; j--)
    {
        if (vc[j] > vc[i])
            break;
    }
    if (i == -1 || j == -1)
        return;
    int temp;
    temp = vc[i];
    vc[i] = vc[j];
    vc[j] = temp;
    int k = vc.size() - 2;
    for (int m = i + 1; m <= (k + i + 1) / 2; m++)
    {
        temp = vc[m];
        vc[m] = vc[k - (m - i - 1)];
        vc[k - (m - i - 1)] = temp;
    }
}

int notUsed(vector<int>& used, long long blockNum) {

    int j;
    long long pos = 0;
    for (j = 1; j < used.size(); j++) {
        if (!used[j]) pos++;
        if (blockNum == pos)
            break;
    }
    return j;
}
void permutation_by_num(vector<int>& re, int n, long long num, long long countMarshruts)
{
    vector<int>used(n + 1, 0);
    int temp = n;
    countMarshruts = countMarshruts / temp;
    for (int i = 0; i < n; i++)
    {
        long long blockNum = (num - 1) / countMarshruts + 1;
        int j = notUsed(used, blockNum);
        re[i + 1] = j;
        used[j] = 1;
        num = (num - 1) % countMarshruts + 1;
        temp--;
        if (temp != 0)
            countMarshruts = countMarshruts / temp;
    }
}

void solver(int& res, vector<vector<int>>& re_res, int**& a, int n, long long startPermutation, long long endPermutation, int rank, long long countMarshruts)
{
    res = 0;
    vector<int>re(n + 1, 0);
    re[0] = 0;
    re[n] = 0;
    permutation_by_num(re, n - 1, startPermutation + 1, countMarshruts);
    for (int i = 1; i < re.size(); i++)
        res += a[re[i - 1]][re[i]];
    re_res.push_back(vector<int>(re));
    int res_temp = 0;
    while (startPermutation != endPermutation)
    {
        permutation(re);
        res_temp = 0;
        for (int i = 1; i < re.size(); i++)
            res_temp += a[re[i - 1]][re[i]];
        if (res_temp < res)
        {
            res = res_temp;
            re_res.clear();
            re_res.push_back(re);
        }
        else if (res_temp == res)
        {
            re_res.push_back(re);
        }
        startPermutation++;
    }
}

int main(int argc, char** argv)
{
    setlocale(LC_ALL, "Russian");
    int rank, size;
    double start = 0, stop = 0;
    MPI_Status stat;
    MPI_Request request;
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    int n;
    string file = "ex1.txt";
    if (rank == 0)
    {
        start = MPI_Wtime();
        ifstream ofs(file);
        if (ofs.is_open())
            ofs >> n;
        else
        {
            cout << "----Ошибка открытия файла----" << endl;
            ofs.close();
            MPI_Abort(MPI_COMM_WORLD, 1);
            MPI_Finalize();
            return 0;
        }
        ofs.close();
    }
    MPI_Bcast(&n, 1, MPI_INT, 0, MPI_COMM_WORLD);
    int** a = new int* [n];
    a[0] = new int[n * n];
    for (int i = 1; i < n; i++)
        a[i] = a[i - 1] + n;
    if (rank == 0)
    {
        ifstream ofs(file);
        if (ofs.is_open())
        {
            for (int i = 0; i < n; i++)
                for (int j = 0; j < n; j++)
                    ofs >> a[i][j];
        }
        else
        {
            cout << "----Ошибка открытия файла----" << endl;
            ofs.close();
            MPI_Abort(MPI_COMM_WORLD, 1);
            MPI_Finalize();
            return 0;
        }
        ofs.close();
    }
    MPI_Bcast(&a[0][0], n * n, MPI_INT, 0, MPI_COMM_WORLD);
    vector<vector<int>>re_res(0);
    long long countMarshruts = 1;
    for (int i = 1; i <= n - 1; i++)
        countMarshruts *= i;
    long long startPermutation = rank * (countMarshruts / size + countMarshruts % size);
    long long endPermutation = startPermutation + countMarshruts / size + 1 * int(rank < countMarshruts% size);
    int res = 0;
    int itog = 0;
    solver(res, re_res, a, n, startPermutation, endPermutation, rank, countMarshruts);
    MPI_Allreduce(&res, &itog, 1, MPI_INT, MPI_MIN, MPI_COMM_WORLD);
    stringstream ss;
    ss << "";
    if (res == itog)
    {
        for (int i = 0; i < re_res.size(); i++)
        {
            for (int j = 0; j < re_res[i].size(); j++)
                ss << re_res[i][j] << " ";
            ss << endl;
        }
    }
    string textData = ss.str();
    // 3. Параллельная запись в файл
    MPI_File fl;
    string aa = "OUTPUT.txt";
    MPI_File_open(MPI_COMM_WORLD, (char*)aa.data(), MPI_MODE_CREATE | MPI_MODE_WRONLY, MPI_INFO_NULL, &fl);
    // Вычисляем размер данных (в байтах) на каждом процессе
    int dataSize = textData.size();
    // Вычисляем смещение каждого процесса
    MPI_Offset offset = 0;
    MPI_Exscan(&dataSize, &offset, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
    // Записываем данные в файл
    MPI_File_write_at_all(fl, offset, (char*)textData.c_str(), dataSize, MPI_CHAR, MPI_STATUS_IGNORE);
    MPI_File_close(&fl);
    if (rank == 0)
    {
        stop = MPI_Wtime();
        cout << "SIZE:  " << size << endl;
        cout << setprecision(10) << fixed << "TIME:  " << stop - start << endl;
        cout << "Длина оптимального пути:  " << itog << "   " << endl;
        int t;
        ifstream st("OUTPUT.txt");
        int ty = 0;
        if (st.is_open())
        {
            while (!st.eof())
            {
                if (ty % (n + 1) == 0 && ty != 0)
                    cout << endl;
                st >> t;
                cout << t << " ";
                ty++;
            }
        }
        else
        {
            cout << "----Ошибка открытия файла----" << endl;
            st.close();
            MPI_Abort(MPI_COMM_WORLD, 1);
            MPI_Finalize();
            return 0;
        }
        st.close();
    }
    delete[]a[0];
    delete[]a;
    MPI_Finalize();
    return 0;
}