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

       

     static ComplexMath.Complex[] N = new ComplexMath.Complex[] { new Complex(1.333, 0), new Complex(1.59, 0), new Complex(1.5, 0), new Complex(1.4598, 0), new Complex(3.8858, -0.02), };
     Functional func = new Functional(17.472 * Math.PI / 180, 165.747 * Math.PI / 180, Math.PI / 3, N, new double[] { 1.3, 30, 15 }, 6328);
        private void button1_Click(object sender, EventArgs e)
        {
  
        }

    }
}

