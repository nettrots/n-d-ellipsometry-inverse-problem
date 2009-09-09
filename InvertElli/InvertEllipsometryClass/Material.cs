using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using ComplexMath;
using SbB.Diploma;

namespace InvertEllipsometryClass
{
    public class Material
    {
        private int number;
        private double d;
        private Complex angle;
        private Complex n;
        private Complex beta;
        private Matrix l;
        private Complex incAngle;
        private Complex incN;
        private double lambda;

        public Material(Complex n)
        {
            this.n = n;
        }


        public Complex Beta
        {
            get
            {
                beta = getBeta();
                return beta;
            }
        }

        public Complex N
        {
            get { return n; }
            set { n = value; }
        }

        public double D
        {
            get { return d; }
            set { d = value; }
        }

        public Complex Angle
        {
            get { return angle; }
            set { angle = value; }
        }


        public Matrix calcL(Complex incAngle, Complex incN, double d, double lambda)
        {
            this.d = d;
            Func<Complex, Complex> arcsin = (x) =>
            {
                Complex ci = new Complex(0, 1);
                Complex c1 = new Complex(1, 0);//(new Complex(Math.PI, 0)) +
                return - ci * Complex.Log(ci * x + Complex.Sqrt(c1 - x * x));
            };
            //    -i ln(iz +sqrt 1-z^2)}
            this.angle = arcsin(incN / N * Complex.Sin(incAngle));
            Complex c = arcsin(Complex.Sin(incAngle));
            this.incAngle = incAngle;
            this.incN = incN;
            this.lambda = lambda;

            Complex i = new Complex(0, 1);
            return l = new Matrix(new Complex[2][]
                                      {
                                          new[]
                                              {
                                                  Complex.Exp(i*Beta), new Complex(0.0, 0.0)
                                              }
                                          , new[]
                                               {
                                                new Complex(0.0, 0.0), Complex.Exp(-i*Beta)
                                               }
                                      });
        }



        public Complex getBeta()
        {
            return (2 * Math.PI * (d*N / lambda) * Complex.Cos(angle));
        }

    }
    public class Ambient : Material
    {
        public Ambient(Complex n)
            : base(n)
        {
        }
    }
    public class Substract : Material
    {
        public Substract(Complex n)
            : base(n)
        {
        }
    }
    public class Film : Material
    {
        private bool recalc = false;

        public Film(Complex n)
            : base(n)
        {
        }

        public bool Recalc
        {
            get { return recalc; }
            set { recalc = value; }
        }
    }
   
}
