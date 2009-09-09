using System;
using System.Collections;
using System.Collections.Generic;
using System.ComponentModel;
using System.Data;
using System.Drawing;
using System.Linq;
using System.Text;
using System.Threading;
using System.Windows.Forms;
using alglib;
using AP;
using InvertEllipsometryClass;
using Complex = ComplexMath.Complex;
using Math=System.Math;

namespace InvertElli
{
    public partial class Form1 : Form
    {
        public Form1()
        {
            InitializeComponent();
        }


        private string perebor()
        {
         
            double nmin = 1.333, nmax = 1.7,
                   dmin = 80,
                   dmax = 100
                   ,
                   dn = 0.01,
                   dd = 0.5;
            List<double[]> res = new List<double[]>();
            for (double n = nmin; n < nmax; n += dn)
                for (double d = dmin; d < dmax; d += dd)
                    res.Add( new double[] {func.functional(n, d),n, d});
            res.Sort(new ArrComparer());
            string ans = "";
            foreach (var xx in res)
            {
                //ans += "f = " + xx[0] + " n,d = " + xx[1] + "," + xx[2] + "\r\n";
                ans += xx[1] + "\t" + xx[2] + "\t"+xx[0] + "\r\n";

            }
            textBox1.Text += "n\td\tf\t\r\n";
            textBox1.Text += ans;
            return ans;
        }
        public class ArrComparer : IComparer<double[]>
        {
            #region IComparer<double[]> Members

            public int Compare(double[] x, double[] y)
            {
                return x[0].CompareTo(y[0]);
            }

            #endregion
        }
        IComparer<double[]> Compdoublearr()
        {
            return new ArrComparer();
        }
        delegate string SampleDelegate();

        private void LevebAlg()
        {

            double[] aprx = new double[] { 1.4, 90 };
            minlm.lmstate state = new minlm.lmstate();
            minlm.lmreport rep = new minlm.lmreport();
            //state.
            double hN = 0.1;
            double hD = 0.1;

            minlm.minlmfgh(2,  ref aprx, 0.0, 0.0, 100, ref state);

            bool referror;
            while (minlm.minlmiteration(ref state))
            {
                state.f = func.functional(aprx[0], aprx[1]);
                state.g[0] = (func.functional(aprx[0] + hN, aprx[1]) - func.functional(aprx[0] - hN, aprx[1])) /
                                  (2 * hN);
                state.g[1] =  (func.functional(aprx[0], aprx[1] + hD) - func.functional(aprx[0], aprx[1] - hD)) /
                                 (2 * hD);
                state.h[0, 0] = (func.functional(aprx[0] + hN, aprx[1]) - func.functional(aprx[0], aprx[1]) + func.functional(aprx[0] - hN, aprx[1])) /
                                  (hN * hN);
                state.h[1, 1] = (func.functional(aprx[0], aprx[1] + hD) - func.functional(aprx[0], aprx[1]) + func.functional(aprx[0] , aprx[1] - hD)) /
                                  (hD * hD);
                state.h[1, 0] = state.h[0, 1] = (func.functional(aprx[0] + hN, aprx[1] + hD) - func.functional(aprx[0] - hN, aprx[1] + hN) - func.functional(aprx[0] + hN, aprx[1] - hD) + func.functional(aprx[0] - hN, aprx[1] - hD)) /
                                                  (4*hN * hD);
                minlm.minlmresults(ref state, ref aprx, ref rep);

               
            }
            //minlm.minlmresults(ref state, ref aprx, ref rep);
            textBox1.Text += "n = " + (aprx[0]).ToString() + "\td = " + (aprx[1]).ToString() + "\tf= " + func.functional(aprx[0], aprx[1]) + " \r\n";
            //textBox1.Text += "terminationtype: " + rep.terminationtype.ToString() + "\r\n";

        }

        void gradfunc(ref double[] x,ref double f,ref double[] g)
        {
            f=currF.functional(x[1], x[2]);
            g=new double[]
                  {0,
              (currF.functional(x[1]+0.001, x[2])-currF.functional(x[1]-0.001, x[2]))/0.002,       
              (currF.functional(x[1], x[2]+0.001)-currF.functional(x[1], x[2]-0.001))/0.002,
        };
        }

        private Functional currF;
        private void LBFGSalg()
        {
             currF = func;
            double[] aprx = new double[] {0, 1.5, 100 };
            lbfgs.lbfgsstate state=new lbfgs.lbfgsstate();
            lbfgsb.funcgrad = gradfunc;
            int[] nb=new int[]{0,2,2};
            double [] l=new double[]{0,1.3,0};
            double [] u=new double[]{0,1.7,100};
            int info=0;
            lbfgsb.lbfgsbminimize(2, 2, ref aprx, 0.00000000001, 0.00000000001, 0.00000000001, 10000, ref nb, ref l, ref u, ref info);
            textBox1.Text += "n = " + aprx[1].ToString() + "\td =  " + aprx[2].ToString() + "\tf = " + func.functional(aprx[1],aprx[2])+"\r\n";
            //lbfgs.minlbfgs(2,2,aprx,0.0000001,0.0000001,0.0000000001,100,0,state);
        }

        double f(double[] x)
        {
            return currF.functional(x[1],x[2]);
        }
     private void principleAxis()
     {
         currF = func;
         double[] aprx = new double[] {0, 1, 10 };
         principalaxis.f = f;
         principalaxis.principalaxisminimize(2,ref aprx, 0.00000000001, 0.01);
         textBox1.Text += "n = " + aprx[1].ToString() + "\td =  " + aprx[2].ToString() + "\tf  = " + func.functional(aprx[1], aprx[2]) + "\r\n";

     }
     static ComplexMath.Complex[] N = new ComplexMath.Complex[] { new Complex(1, 0), new Complex(1.59, 0), new Complex(1.5, 0), new Complex(1.4598, 0), new Complex(3.8858, -0.02), };
     Functional func = new Functional(23.6688315 * Math.PI / 180, 164.1178038 * Math.PI / 180, Math.PI / 3, N, new double[] { 1.3, 30, 15 }, 6328);
        
        private void button1_Click(object sender, EventArgs e)
        {
       /*     FilmsModel calsData;
            calsData = new FilmsModel(N, new double[] { 95, 30, 15 }, Math.PI / 3, 6328);
            

            Experiment expData;
            //164.1225614, 23.6909338   ||| 23.643, 164.159


            expData = new Experiment(23.6909338 * Math.PI / 180, 164.1225614 * Math.PI / 180, Math.PI / 3);
           
            Complex x1 = expData.Rho;
            Complex x2 = calsData.PhoExp();
            double t = Math.Atan(x1.Modulus);
            double t1 = Math.Atan(x2.Modulus);
            double c = (Complex.Log(x1 / x1.Modulus)).Argument;
            double c1 = (Complex.Log(x2 / x2.Modulus)).Argument;*/
            textBox1.Text += "principleAxis\r\n";
            principleAxis();
            
            textBox1.Text += "LevebAlg\r\n";
            LevebAlg();
            textBox1.Text += "LBFGSalg\r\n";
            LBFGSalg();
            textBox1.Text += "Simple Search";
            perebor();
          //  testlm.testlm_test();
            //SampleDelegate sd=new SampleDelegate(make);
            //textBox1.Text = sd.Invoke();

        }

        private void textBox1_TextChanged(object sender, EventArgs e)
        {

        }
    }
}

