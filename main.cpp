#include <iostream>

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
private:
    unsigned int n;
    unsigned int m;
    int **data;

public:
    Matrix(int cols, int rows) {
        n = cols;
        m = rows;

        data = new int *[cols];
        for (int i = 0; i < cols; i++) {
            data[i] = new int[rows];
            for (int j = 0; j < rows; j++) {
                data[i][j] = 0;
            }
        }
    }

    void print() {
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < m; j++) {
                cout << data[i][j] << " ";
            }
            cout << endl;
        }
    }

    Matrix getIdentity(){
        Matrix identity = Matrix(this->n, this->m);
        int size = min(identity.n, identity.m);
        for (int i = 0; i < size; i++){
            for (int j = 0; j < size; j++){
                if (i == j){
                    identity.setValue(i, j, 1);
                }
                else{
                    identity.setValue(i, j, 0);
                }
            }
        }
    }

    Matrix getElimination(int r, int c){
        int pivot = getValue(r, r);
    }

    void setValue(int i, int j, int value) {
        data[i][j] = value;
    }

    int getValue(int i, int j){
        return data[i][j];
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
    }

    friend istream &operator>>(istream &is, Matrix &matrix) {
        for (int i = 0; i < matrix.n; i++) {
            for (int j = 0; j < matrix.m; j++) {
                int t;
                is >> t;
                matrix.setValue(i, j, t);
            }
        }
        return is;
    }

    friend ostream &operator<<(ostream &os, Matrix &matrix) {
        for (int i = 0; i < matrix.n; i++){
            for (int j = 0; j < matrix.m; j++){
                os << matrix.getValue(i, j) << " ";
            }
            os << endl;
        }
        return os;
    }


};

int main() {

    int n1;

    cin >> n1;
    Matrix a = Matrix(n1, n1);
    cin >> a;


}

