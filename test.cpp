#include <fstream>
#include <iostream>
#include <cstdlib>
#include <ctime>
using namespace std;

int matrix[128][128];

void gen_matrix(int n, int mod) {
    for (int i = 1; i <= n; ++i) {
        for (int j = i + 1; j <= n; ++j) {
            matrix[i][j] = rand() % mod;
            matrix[j][i] = (mod - matrix[i][j]) % mod;
        }
    }
}


int main() {
    srand(time(NULL));
    for (int test = 1; test <= 1000; ++test) {
        ofstream fo("data.txt");
        int mod = 10;
        int n = 10;

        gen_matrix(n, mod);

        fo << mod << " " << n << endl;
        for (int i = 1; i <= n; ++i)
            for (int j = 1; j <= n; ++j)
                fo << matrix[i][j] << " \n"[j == n];
        fo.close();

        if (system("./main < data.txt > Command > /dev/null")) {
            cout << "Failed test #" << test << endl;
            break;
        }

        if (test % 10 == 0)
            cout << "Passed " << test << " tests" << endl;

    }

    return 0;
}