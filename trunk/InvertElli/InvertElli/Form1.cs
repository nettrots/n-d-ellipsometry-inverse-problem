using System;
using System.Linq;
using System.Collections.Generic;
using System.Globalization;
using System.Threading;
using System.Windows.Forms;
using ChartDirector;
using FieldsConnector;
using Graphics;
using InvertEllipsometryClass;
using InvertEllipsometryClass.Optimisation_Algorithms;
using InvertEllipsometryClass.Optimisation_Algorithms.SimpleSearch;
using Complex = ComplexMath.Complex;
using Math=System.Math;
using TextBox=System.Windows.Forms.TextBox;

namespace InvertElli
{
    public partial class Form1 : Form
    {

        private FieldsList fieldsList = new FieldsList();
        private Functional func;
        private SimpleSearchFMW algorythm;
        private NetChart netChart;
        public Form1()
        {
            Thread.CurrentThread.CurrentCulture = new CultureInfo("en-US");
            InitializeComponent();
            netChart=new NetChart(ChartViewer);
            firtsInit();
            netChart.Init();
        }
        private void firtsInit()
        {
       

            ComplexMath.Complex[] N = new []
                                          {
                                              new Complex(amb_n.Text, amb_k.Text),
                                              new Complex(Functional.isParametr, Functional.isParametr),
                                              new Complex(l2_n.Text, l2_k.Text),
                                              new Complex(l1_n.Text, l2_k.Text),
                                              new Complex(subst_n.Text, subst_k.Text)
                                              ,
                                          };


            func = new Functional(
               Convert.ToDouble(psi_tb.Text),
               Convert.ToDouble(delta_tb.Text),
               Convert.ToDouble(aoi.Text),
               N,
               new []
                    {
                        Functional.isParametr,
                        Convert.ToDouble(l2_d.Text),
                        Convert.ToDouble(l1_d.Text)
                    }
               , 6328);
            fieldsList.addBond(FieldConnector.Tie(psi_tb, func, "Text", "Psi", "Psi"));
            fieldsList.addBond(FieldConnector.Tie(delta_tb, func, "Text", "Delta", "Delta"));
            fieldsList.addBond(FieldConnector.Tie(aoi, func, "Text", "IncidentAngle", "IncidentAngle"));

            algorythm = new SimpleSearchFMW(func);


            fieldsList.addBond(FieldConnector.Tie(nmax_b, algorythm, "Text", "Nmax", "Nmax"));
            fieldsList.addBond(FieldConnector.Tie(nstep_b, algorythm, "Text", "Dn", "Dn"));
            fieldsList.addBond(FieldConnector.Tie(nmin_b, algorythm, "Text", "Nmin", "Nmin"));
            fieldsList.addBond(FieldConnector.Tie(dmax_b, algorythm, "Text", "Dmax", "Dmax"));
            fieldsList.addBond(FieldConnector.Tie(dstep_b, algorythm, "Text", "Dd", "Dd"));
            fieldsList.addBond(FieldConnector.Tie(dmin_b, algorythm, "Text", "Dmin", "Dmin"));
            fieldsList.SyncFromAll();
            
        }

        private void button1_Click(object sender, EventArgs e)
        {
            runSimpleSearch();
        }

  
        private void runSimpleSearch()
        {
            OptimizeResult rez=algorythm.Optimize();
            netChart.ChangeData((List<double[]>)rez.Pack);
            netChart.ChangeState();
            netChart.Draw();
            
        }

        private void onTextChanged(object sender, EventArgs e)
        {
            if ((sender as TextBox).Text!="")
            fieldsList[fieldsList.indexOfInterface(sender)].syncFrom();
        }

        private void ChartViewer_Click(object sender, EventArgs e)
        {
            
        }


        


    }
}

