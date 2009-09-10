using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using ComplexMath;
using InvertEllipsometryClass;

namespace InvertEllipsometryClass
{
    public class Functional
    {
        private Experiment expData;
        private FilmsModel calsData;

        public Functional(double psi, double delta, double incidentAngle, Complex[] N, double[] d, double lambda)
        {
            expData = new Experiment(psi, delta, incidentAngle);
            calsData = new FilmsModel(N, d, incidentAngle, lambda);
        }
        public double functional(double n, double d, ref double psi, ref double delta)
        {
            //if (n <= 0 || d <= 0){ Random r=new Random(1); return r.Next()* 10e30;}
            calsData.changeParams(new Complex(n, 0), d);
            Complex x1 = expData.Rho;
            Complex x2 = calsData.PhoExp();
            double t = Math.Atan(x1.Modulus);
            double t1 = Math.Atan(x2.Modulus);
            double c = x1.Argument;
            double c1 = x2.Argument;
            psi = t1 *180/Math.PI;
            delta = c1 * 180 / Math.PI;
            double k =  (t - t1) * (t - t1), k1 = (c - c1) * (c - c1);
            return
                k + k1;
        }
        public double functional(double n, double d)
        {
            //if (n <= 0 || d <= 0){ Random r=new Random(1); return r.Next()* 10e30;}
            calsData.changeParams(new Complex(n, 0), d);
            Complex x1 = expData.Rho;
            Complex x2 = calsData.PhoExp();
            double t = Math.Atan(x1.Modulus);
            double t1 = Math.Atan(x2.Modulus);
            double c = (Complex.Log(x1 / x1.Modulus)).Argument;
            double c1 = (Complex.Log(x2 / x2.Modulus)).Argument;
           
            double k = (t - t1) * (t - t1), k1 = (c - c1) * (c - c1);
            return
                k + k1;
        }
        public double functional(double n, double d, int psi, int delta)
        {
            //if (n <= 0 || d <= 0){ Random r=new Random(1); return r.Next()* 10e30;}
            calsData.changeParams(new Complex(n, 0), d);
            Complex x1 = expData.Rho;
            Complex x2 = calsData.PhoExp();
            double t = Math.Atan(x1.Modulus);
            double t1 = Math.Atan(x2.Modulus);
            double c = (Complex.Log(x1 / x1.Modulus)).imag;
            double c1 = (Complex.Log(x2 / x2.Modulus)).imag;
            return psi * 100 * (t - t1) * (t - t1) + delta * (c - c1) * (c - c1);
        }

    }
}
