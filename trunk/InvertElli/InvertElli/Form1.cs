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

        #region Methods
        private void perebor()
        {
            vf1 vfcl = new vf1(wcl);
            vf1 vf = new vf1(wr);
            vf1 vf23 = new vf1(sign);
            double nmin = Convert.ToDouble(nmin_b.Text),
                   nmax = Convert.ToDouble(nmax_b.Text),
                   dmin = Convert.ToDouble(dmin_b.Text),
                   dmax = Convert.ToDouble(dmax_b.Text)
                   ,
                   dn = Convert.ToDouble(nstep_b.Text),
                   dd = Convert.ToDouble(dstep_b.Text);
            double psi=0, delta=0;
            List<double[]> res = new List<double[]>();
            double num = (nmax - nmin) / (dn) * (dmax - dmin) / (dd);
            int i=0;
            double per=0;
            for (double n = nmin; n <= nmax; n += dn)
                for (double d = dmin; d <= dmax; d += dd)
                {
                    res.Add(new double[] { func.functional(n, d, ref psi, ref delta), n, d, delta, psi });
                    textBox1.BeginInvoke(vfcl, new object[] { "" });
                    per=i++/num*100;
                    textBox1.BeginInvoke(vf, new object[] { (per).ToString()+"%"});
                }

            res.Sort(new ArrComparer());
            string ans = "";
            
            textBox1.BeginInvoke(vfcl, new object[] { ""});
            
            
            foreach (var xx in res)
            {
                //ans += "f = " + xx[0] + " n,d = " + xx[1] + "," + xx[2] + "\r\n";
                //ans += xx[1] + "\t" + xx[2] + "\t"+xx[0] + "\r\n";
                textBox1.BeginInvoke(vf, new object[] { xx[1] + "\t" + xx[2] + "\t" + xx[0] + "\t" + xx[3] + "\t" + xx[4] + "\r\n" });
                
                //ans += xx[1] + "\t" + xx[2] + "\t" + xx[0] + "\t" + xx[3] + "\t" + xx[4]+ "\r\n";

            }
            textBox1.BeginInvoke(vf23, new object[] { "" });
            //textBox1.Text += "n\td\tf\t\r\n";
            //textBox1.Text += ans;
            //Clipboard.SetText(textBox1.Text);
        }
        private void sign(string text)
        {
            done_l.Text = "DONE";
            return;
        }
        private void wr(string text) {
            textBox1.Text += text;
            return;
        }
        private void wcl(string text)
        {
            textBox1.Text = "";
            return;
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
        #endregion

     private Functional currF;
     Functional func;
     delegate void vf();
     delegate void vf1(string s);
     //23.6688315 164.1178038
        private void button1_Click(object sender, EventArgs e)
        {
            ComplexMath.Complex[] N = new ComplexMath.Complex[] { 
                new Complex(Convert.ToDouble(amb_n.Text), Convert.ToDouble(amb_k.Text)), 
                new Complex(1, 0), 
                new Complex(Convert.ToDouble(l2_n.Text), Convert.ToDouble(l2_k.Text)), 
                new Complex(Convert.ToDouble(l1_n.Text), Convert.ToDouble(l2_k.Text)), 
                new Complex(Convert.ToDouble(subst_n.Text), Convert.ToDouble(subst_k.Text)), };
   
            func = new Functional(
                Convert.ToDouble(psi_tb.Text) * Math.PI / 180, 
                Convert.ToDouble(delta_tb.Text) * Math.PI / 180, 
                Convert.ToDouble(aoi.Text) * Math.PI / 180, 
                N, 
                new double[] { 
                    0, 
                    Convert.ToDouble(l2_d.Text), 
                    Convert.ToDouble(l1_d.Text) }
                , 6328);
            done_l.Text = "-";
           /*   textBox1.Text += "principleAxis\r\n";
            principleAxis();
            
            textBox1.Text += "LevebAlg\r\n";
            LevebAlg();
            textBox1.Text += "LBFGSalg\r\n";
            LBFGSalg();
            textBox1.Text += "Simple Search";*/
            vf vf = new vf(perebor);
            vf.BeginInvoke(null, null);
            
        }

        private void textBox13_TextChanged(object sender, EventArgs e)
        {

        }

        private void anb_n_TextChanged(object sender, EventArgs e)
        {

        }

    }
}

