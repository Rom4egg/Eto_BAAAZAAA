#include <iostream>
#include <vector>
#include <cmath>

#include <functional>
#include <fstream>

using namespace std;

/***
 * j бля всегда столбцы  в рк !!!
 * @param A
 * @param b
 * @return
 */
double D = 120;
double d = 40;
double L = 300;
double V0 = 1e-3;
double rho = 1800;
double nu = 0.4;

vector<double> operator*(vector<vector<double>> A, vector<double> b) {
    vector<double> c(A[0].size(), 0);
    for (int i = 0; i < A[0].size(); i++) {
        for (int j = 0; j < A.size(); j++) {
            c[i] += A[j][i] * b[j];
        }
    }
    return c;
}

vector<double> operator*(double c, vector<double> b) {
    for (int i = 0; i < b.size(); i++) {
        b[i] = c * b[i];
    }
    return b;
}

vector<double> operator+(vector<double> c, vector<double> b) {
    for (int i = 0; i < b.size(); i++) {
        b[i] = c[i] + b[i];
    }
    return b;
}

double max_difference(vector<vector<double>> A, vector<vector<double>> B) {
    int max = 0;
    for (int j = 0; j < A.size(); j++) {
        for (int i = 0; i < A[j].size(); i++) {
            if (max < abs(A[j][i] - B[j][i])) {
                max = abs(A[j][i] - B[j][i]);
            }
        }
    }
    return max;
}

vector<vector<double>> K(double time, double h, vector<double> c, vector<double> y, vector<vector<double>> A,
                         function<vector<double>(double, vector<double>)> f) {
    vector<vector<double>> K0(c.size());
    for (int i = 0; i < K0.size(); i++) {
        K0[i].resize(y.size());
    }
    vector<vector<double>> K(c.size());
    for (int i = 0; i < K.size(); i++) {
        K[i].resize(y.size());
    }
    double error = 1e10;
    while (error > 1e-7) {
        for (int j = 0; j < K.size(); j++) {
            K[j] = f(time + h * c[j], y + h * (K * A[j]));
        }
        error = max_difference(K, K0);
        K0 = K;

    }
    return K;
}

vector<double> step(double time, double h, vector<double> c, vector<double> y, vector<vector<double>> A,
                    function<vector<double>(double, vector<double>)> f, vector<double> b) {
    return y + h * (K(time, h, c, y, A, f) * b);
}

vector<pair<double, vector<double>>>
solve(double start, double end, double h, vector<double> y0, vector<double> c, vector<vector<double>> A,
      vector<double> b,
      function<vector<double>(double, vector<double>)> f) {
    double time = start;
    vector<pair<double, vector<double>>> result;
    vector<double> y = y0;
    result.emplace_back(time, y);
    while (time < end) {
        y = step(time, h, c, y, A, f, b);
        time += h;
        result.emplace_back(time, y);
    }
    return result;
}

double S(vector<double> y, double D, double d, double L) {
    return (2 * M_PI * (D * D - (d + 2 * y[1]) * (d + 2 * y[1])) + 2 * M_PI * (d + 2 * y[1]) * (L - 2 * y[1])) * 1e-6;
}

double V(vector<double> y, double D, double d, double L, double V0) {
    return V0 +
            (M_PI_4 *
           L * (D * D - d * d) - D * D * (L - 2 * y[1]) + (d + 2 * y[1]) * (d + 2 * y[1]) * (L - 2 * y[1])) * 1e-9;
}

double P(vector<double> y, double D, double d, double L, double V0) {
    return 8.31 * 3500 * y[0] / (V(y, D, d, L, V0) * 0.028);
}


vector<double> f(double time, vector<double> y) {
    return {10 / (pow(1e7, nu)) * pow(P(y, D, d, L, V0), nu) * rho * S(y, D, d, L) -
            P(y, D, d, L, V0) * (M_PI * 1e-6 * 22 * 22 / 4) / 265.3190878,
            10 / (pow(1e6, nu)) * pow(P(y, D, d, L, V0), nu)};
}


int main() {
    vector<vector<double>> A = {{0,   0,   0, 0},
                                {0.5, 0,   0, 0},
                                {0,   0.5, 0, 0},
                                {0,   0,   1, 0}};
    vector<double> b = {1. / 6, 1. / 3, 1. / 3, 1. / 6};
    std::vector<double> c = {0, 0.5, 0.5, 1};
    auto result = solve(0, 1, 0.001, {9.626e-5, 0}, c, A, b, f);
//    auto result = solve(0, M_PI, 0.01, {0, 1}, c, A, b, f);
    fstream fout("test.csv", ios::out);
    fout << "time, mass"<<endl;
    for (int i=0; i<result.size();i++){
        fout<<result[i].first<<", "<<result[i].second[1]<<endl;
    }
    cout << result.back().first << ' ' << result.back().second[1] << endl;
    fout.close();
}