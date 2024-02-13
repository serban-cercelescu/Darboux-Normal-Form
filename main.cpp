#include <algorithm>
#include <cassert>
#include <iostream>
#include <numeric>
#include <utility>
#include <vector> 
using namespace std;

int MOD, n;

struct Num {
    int num;

    Num() : num(0) {}

    Num(int n) {
        num = n % MOD;
        if (num < 0)
            num += MOD;
    }
};

Num operator+(Num a, Num b) {
    int res;
    res = a.num + b.num;
    if (res > MOD)
        res -= MOD;
    return Num { res };
}

Num operator-(Num a, Num b) {
    int res;
    res = a.num - b.num;
    if (res < 0)
        res += MOD;
    return Num { res };
}

Num operator*(Num a, Num b) {
    return Num { 1LL * a.num * b.num % MOD };
}

bool operator==(Num a, Num b) {
    return a.num == b.num;
}

bool operator!=(Num a, Num b) {
    return a.num != b.num;
}

struct Matrix {
    vector<vector<Num>> mx;
    
    Matrix() {}

    Matrix(int n, int m) {
        mx.resize(n);
        for (int i = 0; i < n; ++i) {
            mx[i].resize(m);
        }
    }

    vector<Num> &operator[](int i) {
        return mx[i];
    }

    void rowop(int i, int j, Num k) {
        for (int l = 0; l < mx[i].size(); ++l) {
            mx[j][l] = mx[j][l] + mx[i][l] * k;
        }
    }

    void colop(int i, int j, Num k) {
        for (int l = 0; l < mx.size(); ++l) {
            mx[l][j] = mx[l][j] + mx[l][i] * k;
        }
    }

    void swap_rows(int i, int j) {
        swap(mx[i], mx[j]);
    }

    void swap_cols(int i, int j) {
        for (int k = 0; k < mx.size(); ++k) {
            swap(mx[k][i], mx[k][j]);
        }
    }
};

Matrix operator * (Matrix a, Matrix b) {
    Matrix res(a.mx.size(), b.mx[0].size());
    for (int i = 0; i < a.mx.size(); ++i) {
        for (int j = 0; j < b.mx[0].size(); ++j) {
            for (int k = 0; k < a.mx[0].size(); ++k) {
                res[i][j] = res[i][j] + a[i][k] * b[k][j];
            }
        }
    }
    return res;
}

Matrix transpose(Matrix a) {
    Matrix res(a.mx[0].size(), a.mx.size());
    for (int i = 0; i < a.mx.size(); ++i) {
        for (int j = 0; j < a.mx[0].size(); ++j) {
            res[j][i] = a[i][j];
        }
    }
    return res;
}

bool operator==(Matrix a, Matrix b) {
    if (a.mx.size() != b.mx.size() || a.mx[0].size() != b.mx[0].size())
        return false;
    for (int i = 0; i < a.mx.size(); ++i) {
        for (int j = 0; j < a.mx[0].size(); ++j) {
            if (a[i][j] != b[i][j])
                return false;
        }
    }
    return true;
}

bool is_skew_symmetric(Matrix a) {
    if (a.mx.size() != a.mx[0].size())
        return false;
    for (int i = 0; i < a.mx.size(); ++i) {
        for (int j = 0; j < i; ++j) {
            if (a[i][j] != 0-a[j][i])
                return false;
        }
    }
    return true;
}

void print(Matrix a) {
    for (int i = 0; i < a.mx.size(); ++i) {
        for (int j = 0; j < a.mx[0].size(); ++j) {
            cout << a[i][j].num << ' ';
        }
        cout << endl;
    }
}


// Utility functions end here
// ------------------------------------------------------------------------------------------------------------------------------------------------
// Algorithm starts here


void basis_swap(Matrix &form, Matrix &u, int i, int j) {
    form.swap_cols(i, j);
    form.swap_rows(i, j);
    u.swap_cols(i, j);
}

/*
    Given basis elements e_i, e_j, e_k, performs elementary row operations on e_i and e_j such that the resulting basis elements e_i', e_j' saitsfy
    mu(e_i', e_k) = gcd(mu(e_i, e_k), mu(e_j, e_k))
    mu(e_j', e_k) = 0
*/
void euclidean(Matrix &form, Matrix &u, int i, int j, int k) {
    swap(i, j);
    int a = form[i][k].num;
    int b = form[j][k].num;
    
    while (b != 0) {
        int r = a / b;

        form.rowop(j, i, Num(-r));
        form.colop(j, i, Num(-r));
        u.colop(j, i, Num(-r));
    
        basis_swap(form, u, i, j);

        a%= b;
        swap(a, b);
    }

    assert(is_skew_symmetric(form));
}

void make_chains(Matrix &form, Matrix &u, vector<int> indices) {
    for (int i = 0; i < indices.size(); ++i) {
        for (int j = int(indices.size()) - 2; j > i; --j) {
            euclidean(form, u, indices[j + 1], indices[j], indices[i]);
        }
    }
}


void assert_chain_form(Matrix &form) {
    assert(is_skew_symmetric(form));
    for (int i = 0; i < n; ++i)
        for (int j = i + 2; j < n; ++j)
            assert(form[i][j] == 0);
}

vector<int> gather_chain(Matrix form) {
    assert_chain_form(form);
    vector<int> res;

    res.push_back(0);
    for (int i = 1; i < n; ++i) {
        if (form[i][i - 1] != 0) {
            res.push_back(i);
        }
        else {
            if (res.size() >= 3)
                return res;
            res.clear();
            res.push_back(i);
        }

    }
    if (res.size() >= 3)
        return res;
    throw -1;
}

bool is_divisible(int a, int b) {
    for (int i = 0; i < MOD; ++i)
        if (1LL * i * a % MOD == b)
            return true;
    return false;
}

int divide(int a, int b) {
    for (int i = 0; i < MOD; ++i)
        if (1LL * i * a % MOD == b)
            return i;
    assert(false);
}

void break_chain(Matrix &form, Matrix &u, vector<int> indices) {
    assert(indices.size() >= 3);

    if (is_divisible(form[indices[0]][indices[1]].num, form[indices[2]][indices[1]].num)) {
        int d = divide(form[indices[0]][indices[1]].num, form[indices[2]][indices[1]].num);
        form.rowop(indices[0], indices[2], Num(-d));
        form.colop(indices[0], indices[2], Num(-d));
        u.colop(indices[0], indices[2], Num(-d));
        return;
    }


    euclidean(form, u, indices[2], indices[0], indices[1]);
    make_chains(form, u, indices);
    break_chain(form, u, indices);
}

Matrix original_form, form, u;

int main() {
    cin >> MOD >> n;

    form = Matrix(n, n);
    u = Matrix(n, n);
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            int t;
            cin >> t;
            form[i][j] = Num(t);
        }
    }

    assert(is_skew_symmetric(form));

    original_form = form;
    for (int i = 0; i < n; ++i)
        u[i][i] = Num { 1 };

    vector<int> indices(n);
    iota(begin(indices), end(indices), 0);
    make_chains(form, u, indices);


    while (true) {
        try {
            vector<int> chain = gather_chain(form);
            break_chain(form, u, chain);
        } catch (int e) {
            break;
        }
    }

    assert(is_skew_symmetric(form));
    assert_chain_form(form);
    for (int i = 1; i < n - 1; ++i)
        assert(!(form[i][i - 1] != 0 && form[i][i + 1] != 0));

    assert(transpose(u) * original_form * u == form);

    cout << "Darboux form: \n";
    print(form);
    cout << "Change of basis matrix: \n";
    print(u);

    return 0;
}
