using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using alglib;
using InvertEllipsometryClass.Optimisation_Algorithms;

namespace InvertElli.LBFGS
{
    class LBFGS_FMW:OptimisationAlgorythm_FMW

    {
        void gradfunc(ref double[] x, ref double f, ref double[] g)
        {
            f = func.functional(x[1], x[2]);
            g = new double[]
                  {0,
              (func.functional(x[1]+0.001, x[2])-func.functional(x[1]-0.001, x[2]))/0.002,       
              (func.functional(x[1], x[2]+0.001)-func.functional(x[1], x[2]-0.001))/0.002,
        };
        }
        private void LBFGSalg()
        {
            func = func;
            double[] aprx = new double[] { 0, 1.5, 100 };
            lbfgs.lbfgsstate state = new lbfgs.lbfgsstate();
            lbfgsb.funcgrad = gradfunc;
            int[] nb = new int[] { 0, 2, 2 };
            double[] l = new double[] { 0, 1.3, 0 };
            double[] u = new double[] { 0, 1.7, 100 };
            int info = 0;
            lbfgsb.lbfgsbminimize(2, 2, ref aprx, 0.00000000001, 0.00000000001, 0.00000000001, 10000, ref nb, ref l, ref u, ref info);
           // textBox1.Text += "n = " + aprx[1].ToString() + "\td =  " + aprx[2].ToString() + "\tf = " + func.functional(aprx[1], aprx[2]) + "\r\n";
            //lbfgs.minlbfgs(2,2,aprx,0.0000001,0.0000001,0.0000000001,100,0,state);
        }
        public override OptimizeResult Optimize()
        {
            throw new NotImplementedException();
        }
    }
}
