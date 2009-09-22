using System;
using System.Collections;
using System.Collections.Generic;
using ComplexMath;

namespace SbB.Diploma
{
    public class Matrix: List<Vector>,ICloneable
    {
        #region Constructors
        public Matrix() {}
        public Matrix(int m, int n) : this(m, n,new Complex(0.0,0.0)) {}
        public Matrix(int m, int n, Complex element)
            : base(new Vector[m])
        {
            for (int i = 0; i < m; i++)
                this[i] = new Vector(n, element);
        }
        public Matrix(PairInt dim) : this(dim,new Complex(0.0,0.0)) { }
        public Matrix(PairInt dim, Complex element) : this(dim.m, dim.n, element) { }
        public Matrix(Vector[] array): base(array) {}
        public Matrix(Vector vector) : this(new Vector[] {vector}) {}
        public Matrix(Complex[][] elements)
            : base(new Vector[elements.Length]) 
        {
            int m = elements.Length;
            if (m == 0) return;
            int n = elements[0].Length;
            for (int i = 0; i < m; i++)
                if (elements[i].Length != n)
                    throw new Exception("Can't create this matrix");
                else this[i] = new Vector(elements[i]);
        }
        #endregion

        #region Properties
        public string toText
        {
            get
            {
                return this.ToString();
            }
        }
        public PairInt Size
        {
            get 
            {
                if (this.Count == 0)
                    return new PairInt(0, 0);
                return new PairInt(this.Count, (this[0]).Length);
            }
            set 
            {
                if (value.m < 0 || value.n < 0)
                    throw new ArgumentOutOfRangeException();
                for (int i = 0; i < this.Count; i++)
                    (this[i]).Length = value.n;
                for (int i = this.Count; i <= value.m; i++)
                    this.Add(new Vector(value.n));
                for (int i = this.Count; i > value.m; i--)
                    this.RemoveAt(i - 1);
            }
        }
        #endregion

        #region Operators
        public static Matrix operator +(Matrix mLeft, Matrix mRight)
        {
            if ((mLeft.Size.m != mRight.Size.m) || (mLeft.Size.n != mRight.Size.n))
                throw new ArgumentOutOfRangeException();
            Matrix m = new Matrix(mLeft.Size);
            for (int i = 0; i < mLeft.Count; i++)
                m[i] = mLeft[i] + mRight[i];
            return m;
        }
        public static Matrix operator -(Matrix mLeft, Matrix mRight)
        {
            if ((mLeft.Size.m != mRight.Size.m) || (mLeft.Size.n != mRight.Size.n))
                throw new ArgumentOutOfRangeException();
            Matrix m = new Matrix(mLeft.Size);
            for (int i = 0; i < mLeft.Count; i++)
                m[i] = mLeft[i] - mRight[i];
            return m;
        }
        public static Matrix operator *(Matrix mLeft, Matrix mRight)
        {
            if (mLeft.Size.n != mRight.Size.m)
                throw new ArgumentOutOfRangeException();
            Matrix m = new Matrix(new PairInt(mLeft.Size.m, mRight.Size.n));
            for (int i = 0; i < m.Size.m; i++)
                for (int j = 0; j < m.Size.n; j++)
                    m[i][j] = mLeft[i] * mRight.Column(j);
            return m;
        }
        public static Vector operator *(Vector vLeft, Matrix mRight)
        {
            if (vLeft.Length != mRight.Size.m)
                throw new ArgumentOutOfRangeException();
            Vector v = new Vector(mRight.Size.n);
            for (int i = 0; i < v.Length; i++)
                v[i] = vLeft * mRight.Column(i);
            return v;
        }
        public static Vector operator *(Matrix mRight, Vector vLeft)
        {
            if (vLeft.Length != mRight.Size.n)
                throw new ArgumentOutOfRangeException();
            Vector v = new Vector(mRight.Size.m);
            for (int i = 0; i < v.Length; i++)
                v[i] = vLeft * mRight[i];
            return v;
        }
        public static Matrix operator *(Complex k, Matrix mRight)
        {
            Vector[] array = new Vector[mRight.Count];
            for (int i = 0; i < mRight.Count; i++)
                array[i] = k*(mRight[i]);
            return new Matrix(array);
        }
        public static Matrix operator *(Matrix mRight, Complex k)
        {
            Vector[] array = new Vector[mRight.Count];
            for (int i = 0; i < mRight.Count; i++)
                array[i] = k*(mRight[i]);
            return new Matrix(array);
        }
        public static Matrix operator /(Matrix mRight, Complex k)
        {
            Vector[] array = new Vector[mRight.Count];
            for (int i = 0; i < mRight.Count; i++)
                array[i] = (mRight[i])/k;
            return new Matrix(array);
        }
        #endregion

        #region Methods
        public Matrix Transposed()
        {
            Matrix m = new Matrix(this.Size.n, this.Size.m);
            for (int i = 0; i < this.Size.n; i++)
                m[i] = this.Column(i);
            return m;
        }
        public Vector Column(int index)
        {
            Vector v = new Vector(Count);
            for (int i = 0; i < v.Length; i++)
                v[i] = this[i][index];
            return v;
        }
        public void ZeroIn()
        {
            for (int i = 0; i < this.Count; i++)
                this[i].ZeroIn();
        }
        public void RemoveColumn(int index)
        {
            if (index < 0 || index > Size.n) throw new ArgumentOutOfRangeException();
            for(int i=0; i<Size.m; i++)
                (this[i]).RemoveAt(index);
        }
        public void RemoveRow(int index)
        {
            if (index < 0 || index > Size.m) throw new ArgumentOutOfRangeException();
            this.RemoveAt(index);
        }

        public override string ToString()
        {
            string str = "";
            for (int i = 0; i < this.Count; i++)
                str += (this[i]).ToString() + "\r\n";
            return str;
        }
        #endregion

        public object Clone()
        {
            Matrix a= new Matrix(this.Size);
            for (int i = 0; i < this.Size.n; i++)
            {
                for (int j = 0; j < this.Size.m; j++)
                {
                    a[i][j] = this[i][j];
                }
            }
            return a;
        }
    }
}
