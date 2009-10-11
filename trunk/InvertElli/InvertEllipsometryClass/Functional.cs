using System;
using ComplexMath;

namespace InvertEllipsometryClass
{
    public class Functional
    {
        public readonly static double isParametr ;
        static Functional()
        {
            isParametr = 0;
        }
        
        private Experiment expData;
        private FilmsModel calsData;
        private double psi;
        private double delta;
        private double incidentAngle;
        private double lambda;
        private Complex[] n;
        private double[] d;

        public Functional(double psi, double delta, double incidentAngle, Complex[] N, double[] d, double lambda)
        {
            this.n = N;
            this.d = d;
            this.psi = psi * Math.PI / 180;
            this.delta = delta * Math.PI / 180;
            this.incidentAngle = incidentAngle*Math.PI/180;
            this.lambda = lambda;
            initData();
        }
        private void initData()
        {
            expData = new Experiment(psi, delta, incidentAngle);
            calsData = new FilmsModel(n, d, incidentAngle , lambda);
        }

        public double Psi
        {
            get { return psi; }
            set { psi = value * Math.PI / 180; initData(); }
        }

        public double Delta
        {
            get { return delta; }
            set { delta = value * Math.PI / 180; initData(); }
        }

        public double IncidentAngle
        {
            get { return incidentAngle; }
            set { incidentAngle = value * Math.PI / 180; initData(); }
        }

        public double functional(double n, double d, ref double psi, ref double delta)
        {
            //if (n <= 0 || d <= 0){ Random r=new Random(1); return r.Next()* 10e30;}
            calsData.changeParams(new Complex(n, 0), d);
            Complex x1 = expData.Rho;
            Complex x2 = calsData.PhoExp();
            psi = Math.Atan(x2.Modulus) * 180 / Math.PI;
            delta = x2.Argument* 180 / Math.PI;
            double k = (Math.Atan(x1.Modulus) - Math.Atan(x2.Modulus)) * (Math.Atan(x1.Modulus) - Math.Atan(x2.Modulus)), 
                k1 = (x1.Argument - x2.Argument) * (x1.Argument - x2.Argument);
            return
              Math.Log(k + k1);
        }
        public double functional(double n, double d)
        {
            return functional( n,  d, ref psi, ref delta);
        }
      

    }
}
