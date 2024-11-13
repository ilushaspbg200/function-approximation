#pragma once
#include <iostream>
#include <string>
#include <iomanip>
using namespace std;
class Matrix {
public:
	int line; 
	int column;
	double** array;
	Matrix() {
		line = 0;
		column = 0;
		array = nullptr;
	}
	Matrix(int line, int column) {
		this->line = line;
		this->column = column;
		array = new double* [line];
		for (int i = 0; i < line; i++) {
			array[i] = new double[column];
			for (int j = 0; j < column; j++) {
				array[i][j] = 0;
			}
		}
	}
	Matrix(int line, int column, double** temp) {
		this->line = line;
		this->column = column;
		array = new double* [line];
		for (int i = 0; i < line; i++) {
			array[i] = new double[column];
			for (int j = 0; j < column; j++) {
				array[i][j] = temp[i][j];
			}
		}
	}
	Matrix(int line, int column, double* temp) {
		this->line = line;
		this->column = column;
		array = new double* [line];
		for (int i = 0; i < line; i++) {
			array[i] = new double[1];
			array[i][0] = temp[i];
		}
	}
	Matrix(const Matrix& other) { 
		this->line = other.line;
		this->column = other.column;
		this->array = new double* [line];
		for (int i = 0; i < line; i++) {
			array[i] = new double[column];
			for (int j = 0; j < column; j++)
				array[i][j] = other.array[i][j];
		}


	}
	void SetMatrix();
	void ShowMyMatrix();
	void operator = (const Matrix& other);
	~Matrix() {
		for (int i = 0; i < line; i++) {
			delete[] array[i];
		}
		delete[] array;
	}
	Matrix operator +(const Matrix& other);
	Matrix operator *(const Matrix& other);
	Matrix operator -(const Matrix& other);
	bool operator ==(const Matrix& other);
	bool operator !=(const Matrix& other);
	Matrix operator *(const double num);
	Matrix inverseMatrix();
	Matrix Transp();

};


