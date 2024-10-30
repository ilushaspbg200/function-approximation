
#include <iostream>
#include <iomanip>
#include<vector>
#include<string>
#include <sstream>
#include<algorithm>
#include<cmath>
#include <stdexcept>
#include <numeric>
#include "matplotlibcpp.h"
#include "MatrixClass.h"
#include <limits>  // Для std::numeric_limits
using namespace std;
namespace plt = matplotlibcpp;
const double a = -1.0;
const double b = 1.0;
//изначальная функция
double f(double x) {
    return x * x * sin(x);
}

double* gauss(double** A, double* b, int n)
{
    for (int i = 0; i < n; ++i) {
        // Поиск максимального элемента
        int max_row = i;
        for (int k = i + 1; k < n; ++k) {
            if (std::abs(A[k][i]) > std::abs(A[max_row][i])) {
                max_row = k;
            }
        }
        std::swap(A[i], A[max_row]);
        std::swap(b[i], b[max_row]);
        for (int k = i + 1; k < n; ++k) {
            double factor = A[k][i] / A[i][i];
            for (int j = i; j < n; ++j) {
                A[k][j] -= factor * A[i][j];
            }
            b[k] -= factor * b[i];
        }
    }
    double* x = new double[n];
    for (int i = n - 1; i >= 0; --i) {
        x[i] = b[i];
        for (int j = i + 1; j < n; ++j) {
            x[i] -= A[i][j] * x[j];
        }
        x[i] /= A[i][i];
    }
    return x;
}
double random_double(double min, double max) {
    return min + static_cast<double>(rand()) / (RAND_MAX / (max - min));
}
vector<double> points(int m) {
    vector<double> temp;
    m = m / 3;
    double step = (b - a) / (m - 1);
    for (double i = a; i <= b; i += step) {
        temp.push_back(i);
    }
    temp.push_back(b);
    return temp;
}
Matrix Vandermond(vector<double> x, int n) {
    int m = x.size();
    Matrix E(m, n + 1);
    for (int j = 0; j < m; j++) {
        for (int k = 0; k < n + 1; k++) {
            E.array[j][k] = pow(x[j], k);
        }
    }
    return E;

}
vector<double> normal_equation_method(Matrix A, Matrix B) {
    double* b_temp = new double[B.line];
    vector<double> result;

    for (int i = 0; i < B.line; i++) {
        b_temp[i] = B.array[i][0];
    }
    double* res = gauss(A.array, b_temp, A.line);
    for (int i = 0; i < A.line; i++) {
        result.push_back(res[i]);
    }
    return result;
}
double p_norm(vector<double> coeff, double x, int n) {
    double res = 0;
    for (int i = 0; i < n; i++) {
        res = res + (coeff[i] * pow(x, i));
    }
    return res;
}

double orthogonal_polynomial_method(const std::vector<double>& x, const std::vector<double>& y, int n, double val) {
    size_t m = x.size();
    if (m == 0) {
        throw std::invalid_argument("Input vector 'x' must not be empty");
    }
    if (x.size() != y.size()) {
        throw std::invalid_argument("Input vectors 'x' and 'y' must have the same size");
    }

    std::vector<std::vector<double>> q(n + 1, std::vector<double>(m, 0.0));
    std::vector<double> q_val(n + 1, 0.0);

    // Initialize the first terms
    std::fill(q[0].begin(), q[0].end(), 1.0);  // Fill q[0] with 1s
    q_val[0] = 1.0;

    // Calculate the mean of x
    double x_sum = std::accumulate(x.begin(), x.end(), 0.0);
    double x_mean = x_sum / m;

    // Set the second terms
    for (size_t i = 0; i < m; ++i) {
        q[1][i] = x[i] - x_mean;
    }
    q_val[1] = val - x_mean;

    // Calculate the rest of the orthogonal polynomials
    for (int j = 1; j < n; ++j) {
        double alpha = 0.0;
        double beta = 0.0;

        double qj_sq_sum = 0.0;
        for (auto& elem : q[j]) {
            qj_sq_sum += elem * elem;
        }
        double qj_1sq_sum = 0.0;
        for (auto& elem : q[j - 1]) {
            qj_1sq_sum += elem * elem;
        }

        for (size_t i = 0; i < m; ++i) {
            alpha += x[i] * q[j][i] * q[j][i];
            beta += x[i] * q[j][i] * q[j - 1][i];
        }

        alpha /= qj_sq_sum;
        beta /= qj_1sq_sum;

        for (size_t i = 0; i < m; ++i) {
            q[j + 1][i] = x[i] * q[j][i] - alpha * q[j][i] - beta * q[j - 1][i];
        }

        q_val[j + 1] = val * q_val[j] - alpha * q_val[j] - beta * q_val[j - 1];
    }

    // Calculate coefficients 'a'
    std::vector<double> a(n + 1);
    for (int k = 0; k <= n; ++k) {
        double qk_y_sum = 0.0;
        double qk_sq_sum = 0.0;
        for (auto& elem : q[k]) {
            qk_sq_sum += elem * elem;
        }
        for (size_t i = 0; i < m; ++i) {
            qk_y_sum += q[k][i] * y[i];
        }

        a[k] = qk_y_sum / qk_sq_sum;
    }

    double s = 0.0;
    for (int i = 0; i <= n; ++i) {
        s += a[i] * q_val[i];
    }

    return s;
}



int main() {
    setlocale(LC_ALL, "rus");
    int m = 50; // количество точек
    vector<double> x_val = points(m);
    //double eps = 0.000001;
    double eps = 0.1;
    vector<double> x_vec(m);
    vector<double> y_vec(m);
    int cnt = 0;
    int it = 0;
    while (cnt != m) {
        if (cnt + 3 > m) {
            double x = x_val[it - 1];
            x_vec[cnt] = x;
            y_vec[cnt] = f(x) + random_double(-eps, eps);
            cnt += 1;
        }
        else {
            double x = x_val[it];
            for (int i = 0; i < 3; i++) {
                x_vec[cnt] = x;
                y_vec[cnt] = f(x) + random_double(-eps, eps);
                cnt += 1;
            }
            it += 1;
        }
    }
    Matrix Y_VEC(y_vec.size(), 1);
    for (int j = 0; j < y_vec.size(); j++) {
        Y_VEC.array[j][0] = y_vec[j];
    }

    cout << "n\tR_norm\tR_ort" << endl;
    for (int n = 1; n < 6; n++) {
        Matrix E = Vandermond(x_vec, n);
        Matrix E_T = E.Transp();
        Matrix E_TE = E_T * E;
        Matrix B = E_T * Y_VEC;
        vector<double> check = x_vec;
        double sum_norm = 0;
        double sum_ort = 0;
        vector<double> coeffs = normal_equation_method(E_TE, B);
        for (int j = 0; j < check.size(); j++) {
            sum_norm += (y_vec[j] - p_norm(coeffs, check[j], coeffs.size())) * (y_vec[j] - p_norm(coeffs, check[j], coeffs.size()));
            sum_ort += (y_vec[j] - orthogonal_polynomial_method(x_vec, y_vec, n, check[j])) * (y_vec[j] - orthogonal_polynomial_method(x_vec, y_vec, n, check[j]));
        }
        cout << setprecision(15) << n << "\t" << sum_norm << "\t" << sum_ort;
        cout << endl;
    }
    cout << endl;
    int k = 1000;
    vector<double> check_points;
    double step = (b - a) / (k - 1);
    for (double i = a; i <= b; i += step) {
        check_points.push_back(i);
    }

    for (int n = 1; n < 6; n++) {
        std::ostringstream ostr;
        ostr << n; 
        std::string theNumberString = ostr.str();
        string st = "n = " + theNumberString;
        vector<double> approx;
        for (int i = 0; i < k; i++) {
            approx.push_back(orthogonal_polynomial_method(x_vec, y_vec, n, check_points[i]));
        }
        plt::plot(check_points, approx);
        plt::scatter(x_vec, y_vec);
        plt::title(st);
        plt::show();
    }




}