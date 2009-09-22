using alglib;
using InvertEllipsometryClass.Optimisation_Algorithms;

namespace InvertElli
{
    public class LevenbergMarquardtAlg:OptimisationAlgorythm_FMW
    {

        double[] aprx = new double[] { 1.4, 90 };

        public override OptimizeResult Optimize()
        {

            minlm.lmstate state = new minlm.lmstate();
            minlm.lmreport rep = new minlm.lmreport();
            //state.
            double hN = 0.1;
            double hD = 0.1;

            minlm.minlmfgh(2, ref aprx, 0.0, 0.0, 100, ref state);

            bool referror;
            while (minlm.minlmiteration(ref state))
            {
                state.f = func.functional(aprx[0], aprx[1]);
                state.g[0] = (func.functional(aprx[0] + hN, aprx[1]) - func.functional(aprx[0] - hN, aprx[1])) /
                                  (2 * hN);
                state.g[1] = (func.functional(aprx[0], aprx[1] + hD) - func.functional(aprx[0], aprx[1] - hD)) /
                                 (2 * hD);
                state.h[0, 0] = (func.functional(aprx[0] + hN, aprx[1]) - func.functional(aprx[0], aprx[1]) + func.functional(aprx[0] - hN, aprx[1])) /
                                  (hN * hN);
                state.h[1, 1] = (func.functional(aprx[0], aprx[1] + hD) - func.functional(aprx[0], aprx[1]) + func.functional(aprx[0], aprx[1] - hD)) /
                                  (hD * hD);
                state.h[1, 0] = state.h[0, 1] = (func.functional(aprx[0] + hN, aprx[1] + hD) - func.functional(aprx[0] - hN, aprx[1] + hN) - func.functional(aprx[0] + hN, aprx[1] - hD) + func.functional(aprx[0] - hN, aprx[1] - hD)) /
                                                  (4 * hN * hD);
                minlm.minlmresults(ref state, ref aprx, ref rep);


            }
            //minlm.minlmresults(ref state, ref aprx, ref rep);
            //textBox1.Text += "n = " + (aprx[0]).ToString() + "\td = " + (aprx[1]).ToString() + "\tf= " + func.functional(aprx[0], aprx[1]) + " \r\n";
            //textBox1.Text += "terminationtype: " + rep.terminationtype.ToString() + "\r\n";
            return new OptimizeResult();
        }

        double f(double[] x)
        {
            return func.functional(x[1], x[2]);
        }
        delegate string SampleDelegate();
        private void LevebAlg()
        {

        }
    }
}