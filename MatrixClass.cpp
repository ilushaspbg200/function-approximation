#include "MatrixClass.h"
void Matrix::SetMatrix()
{
	array = new double* [line];
	cout << "Введите матрицу: ";
	for (int i = 0; i < line; i++) {
		array[i] = new double[column];
		for (int j = 0; j < column; j++)
			cin >> array[i][j];
	}
}

void Matrix::ShowMyMatrix()
{
	for (int i = 0; i < line; i++) {
		for (int j = 0; j < column; j++)
			cout << array[i][j] << " ";
		cout << setprecision(10) << endl;
	}
	cout << setprecision(10) << endl;


}

void Matrix::operator=(const Matrix& other) // присваивание
{
	for (int i = 0; i < line; i++) {
		delete[] array[i];
	}
	delete[] array;
	line = other.line;
	column = other.column;
	array = new double* [line];
	for (int i = 0; i < line; i++) {
		array[i] = new double[column];
		for (int j = 0; j < column; j++) {
			array[i][j] = other.array[i][j];
		}
	}
}


Matrix Matrix::operator+(const Matrix& other)
{
	if ((this->line == other.line) && (this->column == other.column)) {
		Matrix temp;
		temp.line = other.line;
		temp.column = other.column;
		temp.array = new double* [other.line];
		for (int i = 0; i < other.line; i++) {
			temp.array[i] = new double[temp.column];
			for (int j = 0; j < other.column; j++) {
				temp.array[i][j] = this->array[i][j] + other.array[i][j];
			}
		}
		return temp;
	}

}
Matrix Matrix::operator*(const Matrix& other) {
	if (this->column == other.line) {
		Matrix temp;
		temp.line = this->line;
		temp.column = other.column;
		temp.array = new double* [temp.line];
		for (int i = 0; i < temp.line; i++) {
			temp.array[i] = new double[temp.column];
			for (int j = 0; j < temp.column; j++) {
				double p = 0;
				for (int z = 0; z < this->column; z++) {
					p = p + (this->array[i][z] * other.array[z][j]);
				}
				temp.array[i][j] = p;
			}
		}
		return temp;

	}
}
Matrix Matrix::operator-(const Matrix& other) {
	if ((this->line == other.line) && (this->column == other.column)) {
		Matrix temp;
		temp.line = other.line;
		temp.column = other.column;
		temp.array = new double* [temp.line];
		for (int i = 0; i < temp.line; i++) {
			temp.array[i] = new double[temp.column];
			for (int j = 0; j < temp.column; j++) {
				temp.array[i][j] = this->array[i][j] - other.array[i][j];
			}
		}
		return temp;
	}
}
bool Matrix::operator==(const Matrix& other) {
	if ((this->column == other.column) && (this->line == other.line)) {
		int k = 0;
		for (int i = 0; i < line; i++) {
			for (int j = 0; j < column; j++) {
				if (this->array[i][j] != other.array[i][j])
					k = 1;
			}
		}
		if (k == 0)
			return true;
		else
			return false;
	}
	else
		return false;
}
bool Matrix::operator!=(const Matrix& other) {
	if ((this->column == other.column) && (this->line == other.line)) {
		int k = 0;
		for (int i = 0; i < line; i++) {
			for (int j = 0; j < column; j++) {
				if (this->array[i][j] != other.array[i][j])
					k = 1;
			}
		}
		if (k == 0)
			return false;
		else
			return true;
	}
	else
		return true;
}
Matrix Matrix::operator*(const double num)
{
	Matrix temp(line, column);
	for (int i = 0; i < line; i++) {
		for (int j = 0; j < column; j++) {
			temp.array[i][j] = this->array[i][j] * num;
		}
	}
	return temp;
}
//Matrix Matrix::inverseMatrix() {
//	int n = this->line;
//	double** res;
//	double* y = new double[n];
//	double* gau;
//	int i, j, k;
//	res = new double* [n];
//	for (i = 0; i < n; i++) {
//		res[i] = new double[n];
//	}
//	for (i = 0; i < n; i++) {
//
//		for (j = 0; j < n; j++) {
//			if (i == j) y[j] = 1;
//			else y[j] = 0;
//		}
//		gau = gauss(this->array, y, n);
//		for (k = 0; k < n; k++) {
//			res[k][i] = gau[k];
//		}
//	}
//	Matrix RES(n, n, res);
//	return RES;
//}

Matrix Matrix::inverseMatrix()
{
	size_t n = this->line;
	Matrix result(n, n);
	double** temp = new double* [n];
	for (int i = 0; i < n; i++) {
		temp[i] = new double[n];
		for (int j = 0; j < n; j++) {
			temp[i][j] = this->array[i][j];
		}
	}
	for (int i = 0; i < n; ++i) {
		result.array[i][i] = 1.0;
	}

	// Приводим матрицу к верхнетреугольному виду
	for (int i = 0; i < n; ++i) {
		if (temp[i][i] == 0) {
			throw std::invalid_argument("Матрица вырожденная, обратной матрицы не существует");
		}
		for (int j = i + 1; j < n; ++j) {
			double factor = temp[j][i] / temp[i][i];
			for (int k = 0; k < n; ++k) {
				temp[j][k] -= factor * temp[i][k];
				result.array[j][k] -= factor * result.array[i][k];
			}
		}
	}

	// Приводим матрицу к диагональному виду
	for (int i = 0; i < n; ++i) {
		double factor = temp[i][i];
		for (int j = 0; j < n; ++j) {
			temp[i][j] /= factor;
			result.array[i][j] /= factor;
		}
	}

	// Обратный проход метода Гаусса-Жордана
	for (int i = n - 1; i > 0; --i) {
		for (int j = i - 1; j >= 0; --j) {
			double factor = temp[j][i];
			for (int k = 0; k < n; ++k) {
				temp[j][k] -= factor * temp[i][k];
				result.array[j][k] -= factor * result.array[i][k];
			}
		}
	}

	return result;
}

Matrix Matrix::Transp()
{
	Matrix temp(this->column, this->line);
	for (int i = 0; i < this->column; i++) {
		for (int j = 0; j < this->line; j++) {
			temp.array[i][j] = this->array[j][i];
		}
	}
	return temp;
}


