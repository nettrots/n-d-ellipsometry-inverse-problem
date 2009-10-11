using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;

namespace InvertEllipsometryClass.Optimisation_Algorithms.SimpleSearch
{
    public class SimpleSearchFMW : OptimisationAlgorythm_FMW
    {
        #region fields
        private double nmin = 1.40;

        private double nmax ;

        private double dmin ;

        private double dmax ;

        private int dn ;

        private int dd ;

        private double percent = 100;

        private bool sortFlag;
        #endregion

        //Constuctor
        public SimpleSearchFMW(Functional f)
        {
            func = f;
        }

        #region Comparer
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
        #endregion

        #region Properties
        public double Nmin
        {
            get { return nmin; }
            set { nmin = value; }
        }

        public double Dmin
        {
            get { return dmin; }
            set { dmin = value; }
        }

        public double Nmax
        {
            get { return nmax; }
            set { nmax = value; }
        }

        public double Dmax
        {
            get { return dmax; }
            set { dmax = value; }
        }

        public int Dn
        {
            get { return dn; }
            set { dn = value; }
        }

        public int Dd
        {
            get { return dd; }
            set { dd = value; }
        }

        public double Percent
        {
            get { return percent; }
            set { percent = value; }
        }

        public bool SortFlag
        {
            get { return sortFlag; }
            set { sortFlag = value; }
        }

        #endregion


        public override OptimizeResult Optimize()
        {
            List<double[]> res = new List<double[]>();
            double deln = (Nmax - Nmin)/dn;
            double deld = (Dmax-Dmin) / dd;
            for (double n = Nmin; n <= Nmax; n += deln)
                for (double d = Dmin; d <= Dmax; d += deld)
                    res.Add(new double[] { n, d, func.functional(n, d, ref psi, ref delta), delta, psi });
            if(sortFlag)
            {
                res.Sort(new ArrComparer());
                if (percent != 100)
                    res.RemoveRange((int)(res.Count * percent / 100), res.Count - 1);
            }
            OptimizeResult pack=new OptimizeResult();
            pack.Pack = res;
            pack.OType = typeof (List<double[]>);
//            string ans = "";
//            foreach (var xx in res)
//            {
//                ans += xx[1] + "\t" + xx[2] + "\t" + xx[0] + "\t" + xx[3] + "\t" + xx[4] + "\r\n";
//
//            }
            //     textBox1.Text += "n\td\tf\t\r\n";
            //    textBox1.Text += ans;
            return pack;
        }
    }
}
