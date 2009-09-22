using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using alglib;

namespace InvertEllipsometryClass.Optimisation_Algorithms.PrincipalAxis
{
    class PrinccipalAxisFMW:OptimisationAlgorythm_FMW
    {
        double f(double[] x)
        {
            return func.functional(x[1], x[2]);
        }
        private void principleAxis()
        {
            func = func;
            double[] aprx = new double[] { 0, 1, 10 };
            principalaxis.f = f;
            principalaxis.principalaxisminimize(2, ref aprx, 0.00000000001, 0.01);
            //textBox1.Text += "n = " + aprx[1].ToString() + "\td =  " + aprx[2].ToString() + "\tf  = " + func.functional(aprx[1], aprx[2]) + "\r\n";

        }
        public override OptimizeResult Optimize()
        {
            throw new NotImplementedException();
        }
    }
}
