#include <iostream>
#include <iomanip>

using namespace std;

class MatrixDimensionMismatch : public exception {
private:
    char *message;

public:
    MatrixDimensionMismatch(char *msg) : message(msg) {}

    char *what() {
        return message;
    }
};


class Matrix {
protected:
    unsigned int n;
    unsigned int m;
    double **data;

public:
    Matrix(int cols, int rows) {
        n = cols;
        m = rows;

        data = new double *[cols];
        for (int i = 0; i < cols; i++) {
            data[i] = new double [rows];
            for (int j = 0; j < rows; j++) {
                data[i][j] = 0;
            }
        }
    }


    int getN() {
        return n;
    }

    int getM() {
        return m;
    }


    void setValue(int i, int j, double value) {
        data[i][j] = value;
    }

    double getValue(int i, int j) {
        return data[i][j];
    }

//    double *getRow(int i) {
//        double ans[n];
//        for (int k = 0; k < n; k++) {
//            ans[k] = data[i][k];
//        }
//
//        return ans;
//    }

    Matrix getCopy() {
        Matrix ans = Matrix(getN(), getM());
        for (int i = 0; i < getN(); i++) {
            for (int j = 0; j < getM(); j++) {
                ans.setValue(i, j, getValue(i, j));
            }
        }
        return ans;
    }

    Matrix operator+(Matrix const secondMatrix) const {
        if (n != secondMatrix.n || m != secondMatrix.m) {
            throw MatrixDimensionMismatch((char *) "Error: the dimensional problem occurred");
        }
        Matrix answer = Matrix(n, m);
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < m; j++) {
                answer.data[i][j] = data[i][j] + secondMatrix.data[i][j];
            }
        }
        return answer;
    }

    Matrix operator-(Matrix const secondMatrix) const {
        if (n != secondMatrix.n || m != secondMatrix.m) {
            throw MatrixDimensionMismatch((char *) "Error: the dimensional problem occurred");
        }
        Matrix answer = Matrix(n, m);
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < m; j++) {
                answer.data[i][j] = data[i][j] - secondMatrix.data[i][j];
            }
        }
        return answer;
    }

    Matrix operator*(Matrix const secondMatrix) const {
        if (m != secondMatrix.n) {
            throw MatrixDimensionMismatch((char *) "Error: the dimensional problem occurred");
        }
        Matrix answer = Matrix(n, secondMatrix.m);

        for (int i = 0; i < n; i++) {
            for (int j = 0; j < secondMatrix.m; j++) {

                for (int k = 0; k < secondMatrix.n; k++) {
                    answer.setValue(i, j, answer.data[i][j] + data[i][k] * secondMatrix.data[k][j]);
                }
            }

        }

        return answer;
    }

    Matrix transpose() {

        Matrix answer = Matrix(m, n);

        for (int i = 0; i < n; i++) {
            for (int j = 0; j < m; j++) {
                answer.setValue(j, i, data[i][j]);
            }
        }

        return answer;

    }

    Matrix &operator=(Matrix const &otherMatrix) {
        data = otherMatrix.data;
        n = otherMatrix.n;
        m = otherMatrix.m;
        return *this;
    }

    friend istream &operator>>(istream &is, Matrix &matrix) {
        for (int i = 0; i < matrix.n; i++) {
            for (int j = 0; j < matrix.m; j++) {
                double t;
                is >> t;
                matrix.setValue(i, j, t);
            }
        }
        return is;
    }

    friend ostream &operator<<(ostream &os, Matrix &matrix) {
        for (int i = 0; i < matrix.n; i++) {
            for (int j = 0; j < matrix.m; j++) {
                os << setprecision (2) << fixed << matrix.getValue(i, j) << " ";
            }
            os << endl;
        }
        return os;
    }

    double getDeterminant();

};

class SquareMatrix : public Matrix {
public:
    SquareMatrix(int n) : Matrix(n, n) {

    }
};

class IdentityMatrix : public SquareMatrix {
public:
    IdentityMatrix(int n) : SquareMatrix(n) {
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < n; j++) {
                if (i == j) {
                    setValue(i, j, 1);
                }
            }
        }
    }

    IdentityMatrix(Matrix matrix) : SquareMatrix(max(matrix.getN(), matrix.getM())) {
        int maxSize = max(matrix.getN(), matrix.getM());
        int minSize = min(matrix.getN(), matrix.getM());
        for (int i = 0; i < maxSize; i++) {
            for (int j = 0; j < maxSize; j++) {
                if (i == j && i < minSize) {
                    setValue(i, j, 1);
                }
            }
        }
    }
};

class EliminationMatrix : public IdentityMatrix {
public:
    EliminationMatrix(Matrix matrix, int i, int j) : IdentityMatrix(matrix.getN()) {
        i--;
        j--;
        double e = -matrix.getValue(i, j) / matrix.getValue(j, j);
        setValue(i, j, e);
        // TODO could be some problems with double type
    }
};

class PermutationMatrix : public IdentityMatrix {
public:
    PermutationMatrix(Matrix matrix, int i, int j) : IdentityMatrix(matrix.getN()) {
        i--;
        j--;
        setValue(i, i, 0);
        setValue(j, j, 0);
        setValue(i, j, 1);
        setValue(j, i, 1);
    }
};

double Matrix::getDeterminant() {
// TODO correct indexes
    Matrix result = getCopy();
    int steps = 1;

    for (int j = 0; j < getM(); j++) {
        if (j != result.getM() - 1){
            double pivot = result.getValue(j, j);
            int pivot_i = j;
            for (int i = 0; i < getN(); i++){
                if (abs(result.getValue(i, j)) > abs(pivot)){
                    pivot = result.getValue(i, j);
                    pivot_i = i;
                }
            }
            if (pivot_i != j){
                Matrix p = PermutationMatrix(result, pivot_i + 1, j + 1);
                result = p * result;
                cout << "step #" << steps << ": permutation" << endl;
                cout << result;
                steps++;
            }
        }
        for (int i = j + 1; i < getN(); i++) {
            if (result.getValue(i, j) != 0){
                Matrix e = EliminationMatrix(result, i + 1, j + 1);
                result = e * result;
                cout << "step #" << steps << ": elimination" << endl;
                cout << result;
                steps++;
            }
        }
    }

    double det = 1;
    for (int i = 0; i < n; i++) {
        det *= result.getValue(i, i);
    }
    cout << "result:" << endl << setprecision (2) << fixed << det << endl;
    return det;


}

int main() {
    int n;

    cin >> n;
    Matrix a = SquareMatrix(n);
    cin >> a;

    a.getDeterminant();

//    cout << "i :" <<  endl << i << endl << "-------" << endl;
//    cout << "e :" <<  endl << e << endl << "-------" << endl;
//    cout << "e * a :" <<  endl << b << endl << "-------" << endl;
//    cout << "p :" <<  endl << p << endl << "-------" << endl;
//    cout << "p * a :" <<  endl << c << endl << "-------" << endl;



}

