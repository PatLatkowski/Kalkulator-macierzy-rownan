#include <iostream>
#include <cmath>
#include <algorithm>
#include <iterator>
#include <Windows.h>

#define ACCURACY 1e-9

using namespace std;


// Obs³uga macierzy
double ** createMatrix(int n) {
	double ** matrix = new double *[n];
	for (int i = 0; i < n; i++) {
		matrix[i] = new double[n];
	}
	return matrix;
}

void removeMatrix(double ** matrix, int n) {
	for (int i = 0; i < n; i++) {
		free(matrix[i]);
	}
	free(matrix);
}

void fillMatrix(double ** tab, double a1, int n) {
	double a2 = -1, a3 = -1;
	for (int i = 0; i < n; i++) {
		for (int j = 0; j < n; j++) {
			if (j == i) tab[i][j] = a1;
			else if (j == i - 1 && j >= 0) tab[i][j] = a2;
			else if (j == i - 2 && j >= 0) tab[i][j] = a3;
			else if (j == i + 1 && j < n) tab[i][j] = a2;
			else if (j == i + 2 && j < n) tab[i][j] = a3;
			else tab[i][j] = 0;
		}
	}
}

void copyMatrix(double ** m1, double ** m2, int size) {
	for (int i = 0; i < size; i++) {
		for (int j = 0; j < size; j++) {
			m2[i][j] = m1[i][j];
		}
	}
}

void printMatrix(double ** tab, int n) {
	for (int i = 0; i < n; i++) {
		for (int j = 0; j < n; j++) {
			cout << tab[i][j] << " ";
		}
		cout << endl;
	}
	cout << endl;
}


//Obs³uga wektorów
void fillVector_b(double tab[], int n, int f) {
	for (int i = 0; i < n; i++) {
		tab[i] = sin(i * (f + 1));
	}
}

void fillVector_x(double tab[], int n) {
	for (int i = 0; i < n; i++) {
		tab[i] = 1;
	}
}

void copyVector(double * src, double * dest, int size) {
	for (int i = 0; i < size; i++) {
		dest[i] = src[i];
	}
}

void printVector(double tab[], int n) {
	for (int i = 0; i < n; i++) {
		cout << tab[i] << " ";
	}
	cout << endl;
}


// Metody obliczeñ
double normResiduum(double ** matrix, double * vb, double * vx, int size) {
	double * res = new double[size];
	for (int i = 0; i < size; i++) {
		res[i] = 0;
		for (int j = 0; j < size; j++) {
			res[i] += matrix[i][j] * vx[j];
		}
		res[i] -= vb[i];
	}
	double norm = 0;
	for (int i = 0; i < size; i++) {
		norm += pow(res[i], 2);
	}
	norm = sqrt(norm);
	free(res);
	return norm;
}

int Jacobi(double** matrix, double * vector_b, double * vector_x, double size) {
	double * vector_xCopy = new double[size];
	int counter = 0;
	do {
		copyVector(vector_x, vector_xCopy, size);
		for (int i = 0; i < size; i++) {
			double value = vector_b[i];
			for (int j = 0; j < size; j++) {
				if(i != j) value -= matrix[i][j] * vector_xCopy[j];
			}
			value /= matrix[i][i];
			vector_x[i] = value;
		}
		counter++;
	} while (normResiduum(matrix, vector_b, vector_x, size) > ACCURACY && counter < 1000);
	return counter;
}

int GaussSeidel(double** matrix, double * vector_b, double * vector_x, double size) {
	int counter = 0;
	do {
		for (int i = 0; i < size; i++) {
			double value = vector_b[i];
			for (int j = 0; j < size; j++) {
				if (i != j) value -= matrix[i][j] * vector_x[j];
			}
			value /= matrix[i][i];
			vector_x[i] = value;
		}
		counter++;
	} while (normResiduum(matrix, vector_b, vector_x, size) > ACCURACY && counter < 1200);
	return counter;
}

double factorizationLU(double ** matrix, double * vb, double * vx, int size) {
	// Tworzenie macierzy L i U
	double ** matrixL = createMatrix(size);
	double ** matrixU = createMatrix(size);
	copyMatrix(matrix, matrixU, size);
	for (int i = 0; i < size; i++) {
		for (int j = 0; j < size; j++) {
			if(i == j) matrixL[i][j] = 1;
			else matrixL[i][j] = 0;
		}
	}
	for (int k = 0; k < size - 1; k++) {
		for (int j = k + 1; j < size; j++) {
			matrixL[j][k] = matrixU[j][k] / matrixU[k][k];
			for (int i = k; i < size; i++) {
				matrixU[j][i] -= matrixL[j][k] * matrixU[k][i];
			}
		}
	}
	// Tworzenie wektora pomocniczego y
	double * vy = new double[size];
	for (int i = 0; i < size; i++) {
		vy[i] = 0;
		for (int j = 0; j < size; j++) {
			vy[i] += matrixU[i][j] * vx[j];
		}
	}
	// Rozwiazywanie ukladu rownan Ly = b wprzod
	for (int i = 0; i < size; i++) {
		double value = vb[i];
		for (int j = 0; j < i; j++) {
			value -= matrixL[i][j] * vy[j];
		}
		vy[i] = value / matrixL[i][i];
	}
	// Rozwiazywanie ukladu rownan Ux = y wstecz
	for (int i = size - 1; i >= 0; i--) {
		double value = vy[i];
		for (int j = i + 1; j < size; j++) {
			value -= matrixU[i][j] * vx[j];
		}
		vx[i] = value / matrixU[i][i];
	}
	removeMatrix(matrixL, size);
	removeMatrix(matrixU, size);
	return normResiduum(matrix, vb, vx, size);
}


// MAIN
int main() {
	int index;
	cout << "Prosze podac indeks: ";
	cin >> index;
	int matrixSize = index % 100;
	int e = ((index % 1000) - matrixSize) / 100;
	int f = ((index % 10000) - matrixSize - e) / 1000;
	matrixSize += 900;
	cout << matrixSize << " " << e << " " << f << endl;

	double ** matrix = createMatrix(matrixSize);
	double * vector_b = new double[matrixSize];
	double * vector_x = new double[matrixSize];
	double residuum, time;
	int counter; 
	DWORD start;
	// JACOBI
	fillMatrix(matrix, 5 + e, matrixSize);
	fillVector_b(vector_b, matrixSize, f);
	fillVector_x(vector_x, matrixSize);
	start = GetTickCount();
	counter = Jacobi(matrix, vector_b, vector_x, matrixSize);
	time = (double)(GetTickCount() - start) / 1000.0f;
	/*cout << "Metoda Jacobiego" << endl << "Macierz: " << endl;
	printMatrix(matrix, MATRIX_SIZE);
	cout << "Wektor b:" << endl;
	printVector(vector_b, MATRIX_SIZE);*/
	cout << "Wektor x:" << endl;
	printVector(vector_x, matrixSize);
	cout << "Metoda Jacobiego" << endl;
	cout << "Liczba iteracji: " << counter << endl;
	cout << "Czas wykonania: " << time << endl << endl;

	// GAUSS - SEIDEL
	fillMatrix(matrix, 5 + e, matrixSize);
	fillVector_b(vector_b, matrixSize, f);
	fillVector_x(vector_x, matrixSize);
	start = GetTickCount();
	counter = GaussSeidel(matrix, vector_b, vector_x, matrixSize);
	time = (double)(GetTickCount() - start) / 1000.0f;
	/*cout << "Metoda Jacobiego" << endl << "Macierz: " << endl;
	printMatrix(matrix, MATRIX_SIZE);
	cout << "Wektor b:" << endl;
	printVector(vector_b, MATRIX_SIZE);
	cout << "Wektor x:" << endl;
	printVector(vector_x, matrixSize);*/
	cout << "Metoda Gaussa - Seidel" << endl;
	cout << "Liczba iteracji: " << counter << endl;
	cout << "Czas wykonania: " << time << endl << endl;

	// JACOBI dla a1 = 3
	fillMatrix(matrix, 3, matrixSize);
	fillVector_b(vector_b, matrixSize, f);
	fillVector_x(vector_x, matrixSize);
	counter = Jacobi(matrix, vector_b, vector_x, matrixSize);
	cout << "Jacobi dla a1 = 3 " << endl;
	cout << "Liczba iteracji: " << counter << endl;

	// GAUSS SEIDEL dla a1 = 3
	fillMatrix(matrix, 3, matrixSize);
	fillVector_b(vector_b, matrixSize, f);
	fillVector_x(vector_x, matrixSize);
	counter = GaussSeidel(matrix, vector_b, vector_x, matrixSize);
	cout << "Gauss - Seidel dla a1 = 3 " << endl;
	cout << "Liczba iteracji: " << counter << endl;

	// Faktoryzacja LU
	fillMatrix(matrix, 3, matrixSize);
	fillVector_b(vector_b, matrixSize, f);
	fillVector_x(vector_x, matrixSize);
	start = GetTickCount();
	residuum = factorizationLU(matrix, vector_b, vector_x, matrixSize);
	time = (double)(GetTickCount() - start) / 1000.0f;
	/*cout << "Faktoryzacja LU" << endl << "Macierz: " << endl;
	printMatrix(matrix, MATRIX_SIZE);
	cout << "Wektor b: " << endl;
	printVector(vector_b, MATRIX_SIZE);
	cout << "Wektor x: " << endl;
	printVector(vector_x, MATRIX_SIZE);*/
	cout << "Faktoryzacja LU" << endl;
	cout << "Residuum: " << residuum << endl;
	cout << "Czas wykonania: " << time << endl << endl;

	int sizeTab[] = { 100, 500, 1000, 2000, 3000, 5000 };
	for (int i = 0; i < 6; i++) {
		cout << "ROZMIAR: " << sizeTab[i] << endl;
		matrix = createMatrix(sizeTab[i]);
		vector_b = new double[sizeTab[i]];
		vector_x = new double[sizeTab[i]];
		// JACOBI
		fillMatrix(matrix, 5 + e, sizeTab[i]);
		fillVector_b(vector_b, sizeTab[i], f);
		fillVector_x(vector_x, sizeTab[i]);
		start = GetTickCount();
		counter = Jacobi(matrix, vector_b, vector_x, sizeTab[i]);
		time = (double)(GetTickCount() - start) / 1000.0f;
		cout << "Metoda Jacobiego dla rozmiaru: " << sizeTab[i] << endl;
		cout << "Liczba iteracji: " << counter << endl;
		cout << "Czas wykonania: " << time << endl << endl;
		// GAUSS - SEIDEL
		fillMatrix(matrix, 5 + e, sizeTab[i]);
		fillVector_b(vector_b, sizeTab[i], f);
		fillVector_x(vector_x, sizeTab[i]);
		start = GetTickCount();
		counter = GaussSeidel(matrix, vector_b, vector_x, sizeTab[i]);
		time = (double)(GetTickCount() - start) / 1000.0f;
		cout << "Metoda Gaussa - Seidel dla rozmiaru: " << sizeTab[i] << endl;
		cout << "Liczba iteracji: " << counter << endl;
		cout << "Czas wykonania: " << time << endl << endl;
		// FAKTORYZACJA LU
		fillMatrix(matrix, 5 + e, sizeTab[i]);
		fillVector_b(vector_b, sizeTab[i], f);
		fillVector_x(vector_x, sizeTab[i]);
		start = GetTickCount();
		residuum = factorizationLU(matrix, vector_b, vector_x, sizeTab[i]);
		time = (double)(GetTickCount() - start) / 1000.0f;
		cout << "Faktoryzacja LU dla rozmiaru: " << sizeTab[i] << endl;
		cout << "Residuum: " << residuum << endl;
		cout << "Czas wykonania: " << time << endl << endl;

		free(vector_b);
		free(vector_x);
		removeMatrix(matrix, sizeTab[i]);
	}
	return 0;
}
