using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using ComplexMath;

namespace InvertEllipsometryClass
{
    public class Experiment
    {
        private double psi;
        private double delta;
        private double incidentAngle;
        private Complex rho;
        public Experiment(double psi, double delta, double incidentAngle)
        {
            this.psi = psi;
            this.delta = delta;
            this.incidentAngle = incidentAngle;
            calcPho();
        }
        public double Psi
        {
            get { return psi; }
            set { psi = value; }
        }

        public double Delta
        {
            get { return delta; }
            set { delta = value; }
        }

        public double IncidentAngle
        {
            get { return incidentAngle; }
            set { incidentAngle = value; }
        }

        private void calcPho()
        {
            rho = Math.Tan(psi) * Complex.Exp(new Complex(0, delta));
        }
        public Complex Rho
        {
            get
            {
                if (rho == (null)) calcPho();
                return rho;
            }

        }
    }
}
