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
protected:
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

    int getN() {
        return n;
    }

    int getM() {
        return m;
    }


    void setValue(int i, int j, int value) {
        data[i][j] = value;
    }

    int getValue(int i, int j) {
        return data[i][j];
    }

    int *getRow(int i) {
        int ans[n];
        for (int k = 0; k < n; k++) {
            ans[k] = data[i][k];
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
        for (int i = 0; i < matrix.n; i++) {
            for (int j = 0; j < matrix.m; j++) {
                os << matrix.getValue(i, j) << " ";
            }
            os << endl;
        }
        return os;
    }


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
        double e = -data[i][j] / data[j][j];
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

int main() {

    int n;

    cin >> n;
    SquareMatrix a = SquareMatrix(n);
    cin >> a;

    IdentityMatrix i = IdentityMatrix(a);

    cout << i << endl;


}

