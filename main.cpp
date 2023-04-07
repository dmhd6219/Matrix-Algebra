#include <iostream>
#include <iomanip>
#include <cmath>

using namespace std;

class MatrixException : public exception {
private:
    char *message;

public:
    MatrixException(char *msg) : message(msg) {}

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
            data[i] = new double[rows];
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
            throw MatrixException((char *) "Error: the dimensional problem occurred");
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
            throw MatrixException((char *) "Error: the dimensional problem occurred");
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
            throw MatrixException((char *) "Error: the dimensional problem occurred");
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

    Matrix getTranspose() {

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
                double num = matrix.getValue(i, j);
                if (abs(num) < 0.005) {
                    os << setprecision(4) << fixed << abs(matrix.getValue(i, j));
                } else {
                    os << setprecision(4) << fixed << matrix.getValue(i, j);
                }
                if (j != matrix.m - 1) {
                    cout << " ";
                }


            }
            os << endl;
        }
        return os;
    }

    double getDeterminant();

    Matrix getInverse();

    Matrix solveGaussianEquation(Matrix b);

    Matrix solveJacobianEquation(Matrix b, long double epsilon);

    Matrix solveSeidelEquation(Matrix b, long double epsilon);


    static void printAugmented(Matrix left, Matrix right) {

        for (int i = 0; i < left.getN(); i++) {
            for (int j = 0; j < left.getM(); j++) {
                cout << setprecision(2) << fixed << left.getValue(i, j) << " ";
            }
            for (int j = 0; j < right.getM(); j++) {
                cout << setprecision(2) << fixed << right.getValue(i, j) << " ";
            }
            cout << endl;
        }
    }

};

class ColumnVector : public Matrix {
public:
    ColumnVector(int n) : Matrix(n, 1) {

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
        double e = -matrix.getValue(i, j) / matrix.getValue(j, j);
        setValue(i, j, e);
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
    Matrix result = getCopy();
    int steps = 1;

    for (int j = 0; j < getM(); j++) {
        if (j != result.getM() - 1) {
            double pivot = result.getValue(j, j);
            int pivot_i = j;
            for (int i = j; i < getN(); i++) {
                if (abs(result.getValue(i, j)) > abs(pivot)) {
                    pivot = result.getValue(i, j);
                    pivot_i = i;
                }
            }
            if (pivot_i != j) {
                Matrix p = PermutationMatrix(result, pivot_i + 1, j + 1);
                result = p * result;
                cout << "step #" << steps << ": permutation" << endl;
                cout << result;
                steps++;
            }
        }
        for (int i = j + 1; i < getN(); i++) {
            if (result.getValue(i, j) != 0) {
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
    cout << "result:" << endl << setprecision(2) << fixed << det << endl;
    return det;
}

Matrix Matrix::getInverse() {
    if (getN() != getM()) {
        throw MatrixException((char *) "Error: the dimensional problem occurred");
    }
    Matrix left = getCopy();
    Matrix right = IdentityMatrix(getN());

    for (int j = 0; j < getM(); j++) {
        if (j != left.getM() - 1) {
            double pivot = left.getValue(j, j);
            int pivot_i = j;
            for (int i = j; i < getN(); i++) {
                if (abs(left.getValue(i, j)) > abs(pivot)) {
                    pivot = left.getValue(i, j);
                    pivot_i = i;
                }
            }
            if (pivot_i != j) {
                Matrix p = PermutationMatrix(left, pivot_i + 1, j + 1);
                left = p * left;
                right = p * right;
            }
        }
        for (int i = j + 1; i < getN(); i++) {
            if (left.getValue(i, j) != 0) {
                Matrix e = EliminationMatrix(left, i + 1, j + 1);
                left = e * left;
                right = e * right;
            }
        }
    }

    for (int j = getM() - 1; j >= 0; j--) {

        for (int i = j - 1; i >= 0; i--) {
            if (left.getValue(i, j) != 0) {
                Matrix e = EliminationMatrix(left, i + 1, j + 1);
                left = e * left;
                right = e * right;
            }
        }
    }

    for (int i = 0; i < getN(); i++) {
        double value = 1 / left.getValue(i, i);
        left.setValue(i, i, left.getValue(i, i) * value);

        for (int j = 0; j < getN(); j++) {
            right.setValue(i, j, right.getValue(i, j) * value);
        }
    }


    return right;
}

Matrix Matrix::solveGaussianEquation(Matrix b) {
    if (b.getM() != 1) {
        if (getN() != getM()) {
            throw MatrixException((char *) "Error: the dimensional problem occurred");
        }
    }
    Matrix left = getCopy();
    Matrix right = b.getCopy();


    for (int j = 0; j < getM(); j++) {
        if (j != left.getM() - 1) {
            double pivot = left.getValue(j, j);
            int pivot_i = j;
            for (int i = j; i < getN(); i++) {
                if (abs(left.getValue(i, j)) > abs(pivot)) {
                    pivot = left.getValue(i, j);
                    pivot_i = i;
                }
            }
            if (pivot_i != j) {
                Matrix p = PermutationMatrix(left, pivot_i + 1, j + 1);
                left = p * left;
                right = p * right;
            }
        }
        for (int i = j + 1; i < getN(); i++) {
            if (left.getValue(i, j) != 0) {
                Matrix e = EliminationMatrix(left, i + 1, j + 1);
                left = e * left;
                right = e * right;
            }
        }
    }
    for (int j = getM() - 1; j >= 0; j--) {
        for (int i = j - 1; i >= 0; i--) {
            if (left.getValue(i, j) != 0) {
                Matrix e = EliminationMatrix(left, i + 1, j + 1);
                left = e * left;
                right = e * right;
            }
        }
    }

    for (int i = 0; i < getN(); i++) {
        double value = 1 / left.getValue(i, i);

        left.setValue(i, i, left.getValue(i, i) * value);
        right.setValue(i, 0, right.getValue(i, 0) * value);
    }

    return right;
}

Matrix Matrix::solveJacobianEquation(Matrix b, long double epsilon) {
    if (getM() != b.getN()) {
        throw MatrixException("The method is not applicable!");
    }
    for (int i = 0; i < getN(); i++) {
        if (getValue(i, i) == 0) {
            throw MatrixException("The method is not applicable!");
        }
        long double sum = 0;
        for (int j = 0; j < getM(); j++) {
            if (i != j) {
                sum += getValue(i, j);
            }
        }
        if (getValue(i, i) < sum) {
            throw MatrixException("The method is not applicable!");
        }
    }
    Matrix aD = SquareMatrix(n);
    for (int i = 0; i < getN(); i++) {
        aD.setValue(i, i, getValue(i, i));
    }
    Matrix aDInverse = aD.getInverse();

    Matrix alpha = aDInverse * (aD - *this);
    cout << "alpha:" << endl << alpha;

    Matrix beta = aDInverse * b;
    cout << "beta:" << endl << beta;

    Matrix x = beta.getCopy();
    Matrix newx = IdentityMatrix(beta);
    cout << "x(0):" << endl << x;
    int steps = 1;

    while (true) {
        newx = alpha * x + beta;
        Matrix delta = newx - x;
        long double error = 0;
        for (int i = 0; i < delta.getN(); i++) {
            error += powl(delta.getValue(i, 0), 2);
        }
        error = ::sqrt(error);
        cout << "e: " << error << endl;


        cout << "x(" << steps++ << "):" << endl << newx;

        if (error <= epsilon) {
            return x;
        }

        x = newx;
    }
}

Matrix Matrix::solveSeidelEquation(Matrix b, long double epsilon) {
    if (getM() != b.getN()) {
        throw MatrixException("The method is not applicable!");
    }
    for (int i = 0; i < getN(); i++) {
        if (getValue(i, i) == 0) {
            throw MatrixException("The method is not applicable!");
        }
        long double sum = 0;
        for (int j = 0; j < getM(); j++) {
            if (i != j) {
                sum += fabs(getValue(i, j));
            }
        }
        if (fabs(getValue(i, i)) < sum) {
            throw MatrixException("The method is not applicable!");
        }
    }

    Matrix aD = SquareMatrix(n);
    for (int i = 0; i < getN(); i++) {
        aD.setValue(i, i, getValue(i, i));
    }

    Matrix aDInverse = aD.getInverse();

    Matrix alpha = aDInverse * (aD - *this);
    Matrix beta = aDInverse * b;

    Matrix aL = SquareMatrix(n);
    for (int j = 0; j < getM(); j++){
        for (int i = j + 1; i < getN(); i++){
            aL.setValue(i, j, alpha.getValue(i, j));
        }
    }

    Matrix aU = SquareMatrix(n);
    for (int j = getM() - 1; j >= 0; j--) {
        for (int i = j - 1; i >= 0; i--){
            aU.setValue(i, j, alpha.getValue(i, j));
        }
    }

    Matrix IminusB = IdentityMatrix(aL) - aL;
    Matrix IminusBInverse = IminusB.getInverse();

    cout << "beta:" << endl << beta;
    cout << "alpha:" << endl << alpha;
    cout << "B:" << endl << aL;
    cout << "C:" << endl << aU;
    cout << "I-B:" << endl << IminusB;
    cout << "(I-B)_-1:" << endl << IminusBInverse;
    cout << "x(0):" << endl << beta;

    Matrix x = beta.getCopy();
    int steps = 1;
    while (true){
        Matrix newx = IminusBInverse * (aU * x + beta);

        Matrix delta = newx - x;
        long double error = 0;
        for (int i = 0; i < delta.getN(); i++) {
            error += powl(delta.getValue(i, 0), 2);
        }
        error = ::sqrt(error);
        cout << "e: " << error << endl;

        cout << "x(" << steps++ << "):" << endl << newx;

        if (error <= epsilon) {
            return x;
        }

        x = newx;
    }
}

int main() {
    int n;
    cin >> n;

    Matrix a = SquareMatrix(n);
    cin >> a;

    cin >> n;
    Matrix b = ColumnVector(n);
    cin >> b;

    double epsilon;
    cin >> epsilon;

    try {
        a.solveSeidelEquation(b, epsilon);
    }
    catch (MatrixException e) {
        cout << e.what() << endl;
    }
    return 0;

}



