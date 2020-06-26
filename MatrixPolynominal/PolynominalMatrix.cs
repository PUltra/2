using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Numerics;

namespace UltraMatrixPolynominal
{
    /// <summary>
    /// Operations on matrices
    /// </summary>
    public class Matrix
    {

        private Double[,] coefs;

        /// <summary>
        /// Rows number
        /// </summary>
        public int Rows;

        /// <summary>
        /// Columns number
        /// </summary>
        public int Columns;

        /// <summary>
        /// Creates new matrix from Double array
        /// </summary>
        /// <param name="InputCoefs">Input array</param>
        public Matrix(Double[,] InputCoefs)
        {
            Rows = InputCoefs.GetLength(0);
            Columns = InputCoefs.GetLength(1);

            coefs = new Double[Rows, Columns];

            for (int n = 0; n < Rows; n++)
            {
                for (int m = 0; m < Columns; m++)
                {
                    coefs[n, m] = InputCoefs[n, m];
                }
            }

        }

        /// <summary>
        /// Creates new matrix of given value in given dimensions
        /// </summary>
        /// <param name="rows">Rows number</param>
        /// <param name="columns">Columns number</param>
        /// <param name="value">Each element value</param>
        public Matrix(int rows, int columns, Double value)
        {
            Rows = rows;
            Columns = columns;

            coefs = new Double[Rows, Columns];

            for (int n = 0; n < Rows; n++)
            {
                for (int m = 0; m < Columns; m++)
                {
                    coefs[n, m] = value;
                }
            }

        }

        /// <summary>
        /// Creates new square matrix with given value on diagonal
        /// </summary>
        /// <param name="dim">Rows / Columns number</param>
        /// <param name="value">Diagonal elements value</param>
        public Matrix(int dim, Double value)
        {
            Rows = dim;
            Columns = dim;

            coefs = new Double[Rows, Columns];

            for (int n = 0; n < Rows; n++)
            {
                for (int m = 0; m < Columns; m++)
                {
                    if (n == m) { coefs[n, m] = value; }
                }
            }
        }


        private int MaxCharsValue()
        {
            int max = 0;

            for (int n = 0; n < Rows; n++)
            {
                for (int m = 0; m < Columns; m++)
                {
                    if (coefs[n,m].ToString().Length > max)
                    {
                        max = coefs[n,m].ToString().Length;
                    }
                }
            }

            return max;
        }

        private string valuePadLeft(Double value, int MaxChars)
        {
            string valueStr = value.ToString();
           
            return valueStr.PadLeft(MaxChars, ' ');
        }

        /// <summary>
        /// Prints matrix on console
        /// </summary>
        public void Print()
        {
            int max = MaxCharsValue();

            for (int n = 0; n < Rows; n++)
            {
                for (int m = 0; m < Columns; m++)
                {
                    string show = valuePadLeft(coefs[n, m], max);

                    if (m == Columns-1)
                    {
                        Console.WriteLine(show);
                    }
                    else
                    {
                        Console.Write(show + " ");
                    }
                }
            }
        }

        /// <summary>
        /// Get n, m element value
        /// </summary>
        /// <param name="n">Row</param>
        /// <param name="m">Column</param>
        /// <returns>Value</returns>
        public Double Element(int n, int m)
        {
            if (n > Rows || m > Columns)
            {
                throw new Exception("Element coordinates out of range");
            }

            return coefs[n, m];

        }

        /// <summary>
        /// Set n, m element new value
        /// </summary>
        /// <param name="n">Row</param>
        /// <param name="m">Column</param>
        /// <param name="value">New value</param>
        public void SetElement(int n, int m, Double value)
        {
            if (n > Rows || m > Columns)
            {
                throw new Exception("Element coordinates out of range");
            }

            coefs[n, m] = value;
        }

        /// <summary>
        /// Transpose matrix
        /// </summary>
        /// <param name="Matrix">Input matrix</param>
        /// <returns>Transposed matrix</returns>
        public Matrix Transpose()
        {
            Double[,] _MatrixT = new Double[Columns, Rows];

            for (int m = 0; m < Rows; m++)
            {
                for (int n = 0; n < Columns; n++)
                {
                    _MatrixT[n, m] = coefs[m, n];
                }
            }

            return new Matrix(_MatrixT);
        }

        /// <summary>
        /// Matrix trace
        /// </summary>
        /// <returns>Matrix trace</returns>
        public Double Trace()
        {
            if (Rows == Columns)
            {
                Double Trace = 0;

                for (int m = 0; m < Rows; m++)
                {
                    Trace += coefs[m, m];
                }

                return Trace;
            }
            else
            {
                throw new Exception("Matrix must be square.");
            }
        }

        /// <summary>
        /// Adds new matrix to existing one
        /// </summary>
        /// <param name="MatrixB">Real numbers matrix</param>
        /// <returns>Sum</returns>
        public Matrix Add(Matrix MatrixB)
        {

            if (Rows == MatrixB.Rows && Columns == MatrixB.Columns)
            {
                Matrix result = new Matrix(Rows, Columns, 0);

                for (int n = 0; n < Rows; n++)
                {
                    for (int m = 0; m < Columns; m++)
                    {
                        result.SetElement(n, m, (this.Element(n, m) + MatrixB.Element(n, m)));
                    }
                }

                return result;

            }
            else
            {
                throw new Exception("Matrices dimensions must be equal.");
            }
        }

        /// <summary>
        /// Subtracts new matrix from existing one
        /// </summary>
        /// <param name="MatrixB">Real numbers matrix</param>
        /// <returns>Result of subtraction</returns>
        public Matrix Sub(Matrix MatrixB)
        {

            if (Rows == MatrixB.Rows && Columns == MatrixB.Columns)
            {
                Matrix result = new Matrix(Rows, Columns, 0);

                for (int n = 0; n < Rows; n++)
                {
                    for (int m = 0; m < Columns; m++)
                    {
                        result.SetElement(n, m, (this.Element(n, m) - MatrixB.Element(n, m)));
                    }
                }

                return result;

            }
            else
            {
                throw new Exception("Matrices dimensions must be equal.");
            }
        }

        /// <summary>
        /// Calculates matrix determinant using Laplace expansion and algebraic complement
        /// </summary>
        /// <returns>Matrix determinant</returns>
        public Double Determinant()
        {

            if (Rows == Columns)
            {
                Double Det = 0;

                if (Rows == 1)
                {
                    Det = this.Element(0, 0);
                }
                else if (Rows == 2)
                {
                    Det = this.Element(0, 0) * this.Element(1, 1) - this.Element(0, 1) * this.Element(1, 0);
                }
                else
                {
                    Matrix AlgebraicComplement = new Matrix(Rows - 1, Rows - 1, 0);

                    for (int m = 0; m < Rows; m++)
                    {
                        int a = 0;

                        for (int k = 1; k < Rows; k++)
                        {
                            int b = 0;

                            for (int l = 0; l < Rows; l++)
                            {
                                if (l != m)
                                {
                                    AlgebraicComplement.SetElement(a, b, this.Element(k, l));
                                    b++;
                                }
                            }
                            a++;
                        }
                        Det += Math.Pow(-1, m) * this.Element(0, m) * AlgebraicComplement.Determinant();
                    }

                }

                return Det;
            }
            else
            {
                throw new Exception("Matrix must be square.");
            }
        }

        /// <summary>
        /// Multiplies (right-side) existing matrix (multiplicand) by given one (multiplier)
        /// </summary>
        /// <param name="MatrixB">Multiplier matrix</param>
        /// <returns>Result of multiplification</returns>
        public Matrix Multiply(Matrix MatrixB)
        {
            if (Columns != MatrixB.Rows)
            {
                throw new Exception("Number of columns in A matrix must be equal to number of rows in B matrix.");
            }
            else
            {
                Double[,] _MatrixAB = new Double[Rows, MatrixB.Columns];

                for (int m_a = 0; m_a < Rows; m_a++)
                {
                    for (int n_b = 0; n_b < MatrixB.Columns; n_b++)
                    {
                        _MatrixAB[m_a, n_b] = 0;

                        for (int n_a = 0; n_a < MatrixB.Rows; n_a++)
                        {
                            _MatrixAB[m_a, n_b] = _MatrixAB[m_a, n_b] + this.Element(m_a, n_a) * MatrixB.Element(n_a, n_b);
                        }
                    }
                }


                return new Matrix(_MatrixAB);

            }
        }

        /// <summary>
        /// Multiply each matrix element by scalar
        /// </summary>
        /// <param name="Coefficient">Multiplier scalar</param>
        /// <returns>Scalerd matrix</returns>
        public Matrix Scale(Double Coefficient)
        {
            Matrix result = new Matrix(Rows, Columns, 0);

            for (int m = 0; m < Rows; m++)
            {
                for (int n = 0; n < Columns; n++)
                {
                    result.SetElement(m, n, this.Element(m, n) * Coefficient);

                }
            }

            return result;
        }

        /// <summary>
        /// Calculates the inverse matrix using Faddeev-Leverrier algorithm
        /// </summary>
        /// <returns>Inverse matrix</returns>
        public Matrix Inverse()
        {
            if (Rows == Columns)
            {
                if (this.Determinant() != 0)
                {
                    Matrix LastB;
                    Matrix NextB;
                    Matrix Inv = new Matrix(Rows, Columns, 0);
                    Double LastCoeff;
                    Double NextCoeff;

                    LastB = this;
                    LastCoeff = LastB.Trace();

                    for (int m = 2; m < Rows; m++)
                    {
                        LastB = this.Multiply(LastB.Sub(new Matrix(Rows, 1).Scale(LastCoeff)));
                        LastCoeff = LastB.Trace() / m;
                    }

                    NextB = this.Multiply(LastB.Sub(new Matrix(Rows, 1).Scale(LastCoeff)));
                    NextCoeff = NextB.Trace() / Rows;
                    
                    Inv = LastB.Sub(new Matrix(Rows, 1).Scale(LastCoeff)).Scale(1 / NextCoeff);

                    return Inv;
                }
                else
                {
                    throw new Exception("The matrix is not invertible over commutative ring.");
                }
            }
            else
            {
                throw new Exception("Input data matrix must be square.");
            }
        }

        /// <summary>
        /// Calculates matrix's characteristic polynominal using Faddeev-Leverrier algorithm
        /// </summary>
        /// <returns>Characteristic polynominal</returns>
        public Polynominal CharacteristicPolynominal()
        {
            if (Rows == Columns)
            {
                Double[] Coeffs = new Double[Rows + 1];
                Double[] CoeffsSorted = new Double[Rows + 1];
                Matrix B;
                Double LastCoeff;

                Coeffs[0] = 1;

                B = this;
                LastCoeff = B.Trace();
                Coeffs[1] = -LastCoeff;

                for (int m = 2; m < Rows + 1; m++)
                {
                    B = this.Multiply(B.Sub(new Matrix(B.Rows, 1).Scale(LastCoeff)));

                    LastCoeff = B.Trace() / m;
                    Coeffs[m] = -LastCoeff;
                }

                for (int m = 0; m < Rows + 1; m++)
                {
                    CoeffsSorted[m] = Coeffs[Rows - m];
                }

                return new Polynominal(CoeffsSorted);
            }
            else
            {
                throw new Exception("Input data matrix must be square.");
            }
        }

        /// <summary>
        /// Calculating matrix rank
        /// </summary>
        /// <returns>Rank</returns>
        public int Rank()
        {
            int rank = 0;
            Double[,] mat = coefs;

            for (int k = 0; k < Rows; k++)
            {
                for (int i = k + 1; i < Rows; i++)
                {
                    if (mat[k,k] == 0)
                    {
                        mat = RowPivot(mat, k);

                        Console.WriteLine("Pivot!");

                        Matrix y = new Matrix(mat);
                        y.Print();

                    }

                    Double c = mat[i, k] / mat[k, k];

                    for (int k1 = 0; k1 < Columns; k1++)
                    {
                        mat[i, k1] = - mat[k, k1] * c + mat[i, k1];
                    }
                }

                Matrix x = new Matrix(mat);
                x.Print();

                Double sum = 0;

                for (int i = 0; i < Columns; i++)
                {
                    sum += mat[k, i];
                }

                if (sum != 0) { rank++; }
            }


            return rank;
        }

        private Double[,] RowPivot(Double[,] matrix, int k)
        {
            //int row = 0;

            for (int i = k + 1; i < matrix.GetLength(0); i++ )
            {
                if (matrix[i, i] != 0)
                {
                    Double[] x = new Double[matrix.GetLength(1)];

                    for (int j = 0; j < matrix.GetLength(1); j++)
                    {
                        x[j] = matrix[k, j];
                        matrix[k, j] = matrix[i, j];
                        matrix[i, j] = x[j];
                    }
                    break;
                }
            }

            return matrix;
        }
    }

    /// <summary>
    /// Operations on polynominals
    /// </summary>
    public class Polynominal
    {

        private Complex[] coefs;

        /// <summary>
        /// Creates new polynominal of given coeficients
        /// </summary>
        /// <param name="InputCoefs">Array of coeficients starting with free term</param>
        public Polynominal(Double[] InputCoefs)
        {
            coefs = new Complex[InputCoefs.Length];

            for (int i = 0; i < InputCoefs.Length; i++)
            {
                coefs[i] = new Complex(InputCoefs[i], 0);
            }

        }

        /// <summary>
        /// Creates new polynominal of given coeficients
        /// </summary>
        /// <param name="InputCoefs">Array of coeficients starting with free term</param>
        public Polynominal(Complex[] InputCoefs)
        {
            coefs = new Complex[InputCoefs.Length];

            for (int i = 0; i < InputCoefs.Length; i++)
            {
                coefs[i] = InputCoefs[i];
            }
        }

        /// <summary>
        /// Creates new polynominal (zeros) with given degree
        /// </summary>
        /// <param name="Degree">Degree</param>
        public Polynominal(int Degree)
        {
            coefs = new Complex[Degree+1];
        }

        /// <summary>
        /// Gives polynominal degree
        /// </summary>
        /// <returns>Degree</returns>
        public int Degree()
        {
            return coefs.Length-1;
        }

        /// <summary>
        /// Sets element on given position
        /// </summary>
        /// <param name="Position">Position</param>
        /// <param name="Value">New value</param>
        public void SetElement(int Position, Double Value)
        {

            if (Position > coefs.Length-1)
            {
                throw new Exception("Position out of range of polynominal");
            }

            coefs[Position] = new Complex(Value, 0);
        }

        /// <summary>
        /// Sets element on given position
        /// </summary>
        /// <param name="Position">Position</param>
        /// <param name="Value">New value</param>
        public void SetElement(int Position, Complex Value)
        {

            if (Position > coefs.Length - 1)
            {
                throw new Exception("Position out of range of polynominal");
            }

            coefs[Position] = Value;
        }

        /// <summary>
        /// Returns polynominal coeficient on given position
        /// </summary>
        /// <param name="Position">Position</param>
        /// <returns>Element value</returns>
        public dynamic Element(int Position)
        {
            if (Position > coefs.Length - 1)
            {
                throw new Exception("Position out of range of polynominal");
            }

            if (coefs[Position].Imaginary == 0)
            {
                return coefs[Position].Real;
            }
            else
            {
                return coefs[Position];
            }
        }

        /// <summary>
        /// Prints polynominal in symbolic representation
        /// </summary>
        public void Print()
        {
            string print = "";

            for (int i = coefs.Length - 1; i > -1; i--)
            {
                if (coefs[i].Imaginary == 0)
                {
                    if (coefs[i].Real == 0)
                    {
                    }
                    else if (i == coefs.Length - 1)
                    {
                        print += coefs[i].Real + " x^" + i;
                    }
                    else if (i == 0)
                    {
                        if (coefs[i].Real < 0)
                        {
                            print += " " + coefs[i].Real;
                        }
                        else
                        {
                            print += " + " + coefs[i].Real;
                        }

                    }
                    else
                    {
                        print += " + " + coefs[i].Real + " x^" + i;
                    }
                    
                }
                else
                {

                    if (i == coefs.Length - 1)
                    {
                        print += " (" + coefs[i].Real + " " + coefs[i].Imaginary + "i) x^" + i;
                    }
                    else if (i == 0)
                    {
                        print += " (" + coefs[i].Real + " " + coefs[i].Imaginary + "i)";
                    }
                    else
                    {
                        print += " + (" + coefs[i].Real + " " + coefs[i].Imaginary + "i) x^" + i;
                    }

                    
                }
            }

            Console.WriteLine(print);
        }

        /// <summary>
        /// The method calculates polynominal evaluation with given value
        /// </summary>
        /// <param name="Value">Real value</param>
        /// <returns>Value</returns>
        public dynamic Evaluate(Complex Value)
        {
            Complex sum = 0;

            for (int i = 0; i < coefs.Length; i++)
            {
                sum += coefs[i] * Complex.Pow(Value, i);
            }

            if (sum.Imaginary == 0)
            {
                return sum.Real;
            }
            else
            {
                return sum;
            }
        }

        /// <summary>
        /// The method calculates polynominal evaluation with given value
        /// </summary>
        /// <param name="Value">Complex value</param>
        /// <returns>Value</returns>
        public dynamic Evaluate(Double Value)
        {
            Complex sum = 0;
            Complex ValueC = new Complex(Value, 0);

            for (int i = 0; i < coefs.Length; i++)
            {
                sum += coefs[i] * Complex.Pow(ValueC, i);
            }

            if (sum.Imaginary == 0)
            {
                return sum.Real;
            }
            else
            {
                return sum;
            }
        }

        /// <summary>
        /// The method returns product of polynominal addition
        /// </summary>
        /// <param name="poly">Polynominal coefficients (vector starting with free term)</param>
        /// <returns>Polynominal</returns>
        public Polynominal Add(Polynominal poly)
        {
            int M = this.Degree();
            int N = poly.Degree();
            int K = Math.Max(M, N);
            Polynominal Poly1Ext = new Polynominal(K);
            Polynominal Poly2Ext = new Polynominal(K);
            Polynominal Add = new Polynominal(K);

            for (int m = 0; m < M + 1; m++)
            {
                Poly1Ext.SetElement(m, this.Element(m));
            }

            for (int n = 0; n < N + 1; n++)
            {
                Poly2Ext.SetElement(n, poly.Element(n));
            }

            for (int k = 0; k < K + 1; k++)
            {
                Add.SetElement(k, Poly1Ext.Element(k) + Poly2Ext.Element(k));
            }

            return Add;
        }

        /// <summary>
        /// The method returns product of polynominal multiplication by scalar
        /// </summary>
        /// <param name="Value">Real value</param>
        /// <returns>Polynominal</returns>
        public Polynominal Scale(Double Value)
        {
            Polynominal scale = new Polynominal(this.Degree());

            for (int i = 0; i < this.Degree() + 1; i++)
            {
                scale.SetElement(i, this.Element(i) * Value);
            }

            return scale;
        }

        /// <summary>
        /// The method returns product of polynominal multiplication by scalar
        /// </summary>
        /// <param name="Value">Complex value</param>
        /// <returns>Polynominal</returns>
        public Polynominal Scale(Complex Value)
        {
            Polynominal scale = new Polynominal(this.Degree());

            for (int i = 0; i < this.Degree() + 1; i++)
            {
                scale.SetElement(i, this.Element(i) * Value);
            }

            return scale;
        }

        /// <summary>
        /// Subtracts two polynominals
        /// </summary>
        /// <param name="poly">Subtrahend polynominal</param>
        /// <returns>Difference</returns>
        public Polynominal Sub(Polynominal poly)
        {
            return this.Add(poly.Scale(-1));
        }

        /// <summary>
        /// The method calculates polynominal derivative
        /// </summary>
        /// <returns>Polynominal</returns>
        public Polynominal Derivative()
        {
            int Degree = this.Degree();
            int NumCoefs = Degree + 1;
            Polynominal Derivative = new Polynominal(this.Degree()-1);

            if (Degree > 0)
            {
                for (int i = 1; i < NumCoefs; i++)
                {
                    Derivative.SetElement(i - 1, i * this.Element(i));
                }
            }
            else
            {
                return new Polynominal(new Double[1] { 0 });
            }

            return Derivative;
        }

        /// <summary>
        /// The method returns product of polynominal by binominal division
        /// </summary>
        /// <param name="BinominalRoot">Binominal complex root</param>
        /// <returns>Polynominal</returns>
        public Polynominal ByBinominalDivision(Complex BinominalRoot)
        {
            int M = this.Degree() + 1;

            if (M > 1)
            {
                int N = M - 1;
                Complex[] Quotient = new Complex[N];
                Complex[] QuotientSorted = new Complex[N];
                Complex[] CoeffsSorted = new Complex[M];
                Complex Last;

                for (int m = 0; m < M; m++)
                {
                    CoeffsSorted[m] = this.Element(M - m - 1);
                }

                Last = CoeffsSorted[0];
                Quotient[0] = Last;

                for (int m = 1; m < N; m++)
                {
                    Last = CoeffsSorted[m] + BinominalRoot * Last;
                    Quotient[m] = Last;
                }

                Complex Remainder = CoeffsSorted[M - 1] + BinominalRoot * Last;

                for (int n = 0; n < N; n++)
                {
                    QuotientSorted[n] = Quotient[N - n - 1];
                }

                return new Polynominal(QuotientSorted);
            }
            else
            {
                throw new Exception("Given coefficients number is not corresponding to correct polynominal structure.");
            }
        }

        /// <summary>
        /// The method calculates polynominal roots
        /// </summary>
        /// <param name="Seed">Seed root (if unknown any complex number can be given)</param>
        /// <param name="Accuracy">Result accuracy</param>
        /// <returns>Polynominal roots</returns>
        public Complex[] Roots(Complex Seed, Double Accuracy)
        {
            int N = this.Degree();
            Complex[] Roots = new Complex[N];
            int degree = N;

            Polynominal tmpoly = this;

            for (int n = 0; n < N; n++)
            {
                while (true)
                {
                    Complex _tmp0 = tmpoly.Derivative().Evaluate(Seed);
                    Complex _tmp1 = (degree - 1) * Complex.Pow(tmpoly.Derivative().Evaluate(Seed), 2);
                    Complex _tmp2 = degree * tmpoly.Evaluate(Seed) * tmpoly.Derivative().Derivative().Evaluate(Seed);
                    Complex _tmp3 = Complex.Pow((degree - 1) * (_tmp1 - _tmp2), 0.5);
                    Complex _tmp4 = MaximizeResultAbsolute(_tmp0, _tmp3);

                    Seed = Seed - (degree * tmpoly.Evaluate(Seed)) / _tmp4;
                    Complex result = tmpoly.Evaluate(Seed);

                    if (Complex.Abs(result) < Math.Abs(Accuracy))
                    {
                        Roots[n] = Seed;
                        tmpoly = tmpoly.ByBinominalDivision(Seed);
                        degree--;
                        Random _Random = new Random();
                        Seed = -10 + (_Random.NextDouble() * 15);
                        break;
                    }
                }
            }

            return Roots;
        }


        private Complex MaximizeResultAbsolute(Complex Value1, Complex Value2)
        {
            if (Complex.Abs(Value1 + Value2) > Complex.Abs(Value1 - Value2))
            {
                return Value1 + Value2;
            }
            else
            {
                return Value1 - Value2;
            }
        }
    }

}
