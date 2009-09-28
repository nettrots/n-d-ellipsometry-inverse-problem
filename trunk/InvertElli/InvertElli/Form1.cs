using System;
using System.Collections;
using System.Collections.Generic;
using System.ComponentModel;
using System.Data;
using System.Drawing;
using System.Globalization;
using System.Linq;
using System.Text;
using System.Threading;
using System.Windows.Forms;
using alglib;
using AP;
using ChartDirector;
using InvertEllipsometryClass;
using Complex = ComplexMath.Complex;
using Math=System.Math;

namespace InvertElli
{
    public partial class Form1 : Form
    {
        public Form1()
        {
            Thread.CurrentThread.CurrentCulture = new CultureInfo("en-US");
            InitializeComponent();
        }

        List<double[]> res;
        double[] dataX;
        double[] dataY;
        double[] dataZ;
        delegate void chartdel(WinChartViewer viewer, string img);

        private double dn;
        private double dd;
        private double min = 0;
        private double max = 0;
        private void perebor()
        {
            vf1 vfcl = new vf1(wcl);
            vf1 vf = new vf1(wr);
            vf1 vf23 = new vf1(sign);
          

            double psi=0, delta=0;
            res = new List<double[]>();

            int num = (int)((dn) * (dd));
            
            int i=0;
            double per=0;
            for (double n = nmin; n < nmax; n += (nmax-nmin)/(dn))
                for (double d = dmin; d < dmax; d += (dmax-dmin)/(dd))
                {
                    double f = func.functional(n, d, ref psi, ref delta);

                   
                    res.Add(new[] {f, n, d, delta, psi});
                    per = i++*1.0/num*100;
                    double tempd;
                    if ((int)(tempd = Math.Round(per)) % 2 == 0)
                    {
                        textBox1.BeginInvoke(vfcl, new object[] { "" });

                        textBox1.BeginInvoke(vf, new object[] {(tempd).ToString() + "%"});
                    }


                }

            res.Sort(new ArrComparer());
            min = res[0][0];
            max = res[res.Count-1][0];
            string ans = "";
            
            textBox1.BeginInvoke(vfcl, new object[] { ""});

            dataX = new double[i];
            dataY = new double[i];
            dataZ = new double[i];
        
            i = 0;
            foreach (var xx in res)
            {
                dataX[i] = xx[1];
                dataY[i] = xx[2];
                dataZ[i] = xx[0];
                i++;
                //ans += "f = " + xx[0] + " n,d = " + xx[1] + "," + xx[2] + "\r\n";
                //ans += xx[1] + "\t" + xx[2] + "\t"+xx[0] + "\r\n";
            //    textBox1.BeginInvoke(vf, new object[] { xx[1] + "\t" + xx[2] + "\t" + xx[0] + "\t" + xx[3] + "\t" + xx[4] + "\r\n" });
                
                //ans += xx[1] + "\t" + xx[2] + "\t" + xx[0] + "\t" + xx[3] + "\t" + xx[4]+ "\r\n";

            }
            textBox1.BeginInvoke(vf23, new object[] { "" });
            ChartViewer.BeginInvoke(new chartdel(createChart), new object[] { ChartViewer, "" });
            textBox1.BeginInvoke(vf, new object[] { res[0][1] + "\t" + res[0][2] + "\t" + res[0][0] + "\t" + res[0][3] + "\t" + res[0][4] + "\r\n" });

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
  

     private Functional currF;
     Functional func;
        private double nmin;
        private double nmax;
        private double dmin;
        private double dmax;

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
            nmin = Convert.ToDouble(nmin_b.Text);
            nmax = Convert.ToDouble(nmax_b.Text);
            dmin = Convert.ToDouble(dmin_b.Text);
            dmax = Convert.ToDouble(dmax_b.Text);

            dn = Convert.ToDouble(nstep_b.Text);
            dd = Convert.ToDouble(dstep_b.Text);
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

        private XYChart chart;
        //Main code for creating chart.
        //Note: the argument img is unused because this demo only has 1 chart.
        public void createChart(WinChartViewer viewer, string img)
        {
            chart = new XYChart(viewer.Width, viewer.Height, 0xffffff, 0x888888);
            chart.setSize(viewer.Width, viewer.Height);
            chart.setRoundedFrame();


            chart.setPlotArea(30, 30, viewer.Width-60 , viewer.Height -60, 0);

            ContourLayer layer = chart.addContourLayer(dataX, dataY, dataZ);
            chart.setClipping();
            chart.xAxis().setTickDensity(10);
            chart.yAxis().setTickDensity(10);
            // Move the grid lines in front of the contour layer
            chart.getPlotArea().moveGridBefore(layer);

            // Add a color axis (the legend) in which the top center is anchored at
            // (245, 455). Set the length to 330 pixels and the labels on the top
            // side.
            int w = chart.getWidth();
            int h = chart.getHeight();
            ColorAxis cAxis = layer.setColorAxis(0, 0, Chart.TopLeft, h - w - 40,
                Chart.Top);

            // Add a bounding box to the color axis using the default line color as
            // border.
            cAxis.setBoundingBox(Chart.Transparent, Chart.LineColor);


            // Set the color axis range as 0 to 20, with a step every 2 units
            cAxis.setLinearScale(min, max, Math.Abs(max - min) / 100);

            viewer.Image = chart.makeImage();
        }

      

    }
}

