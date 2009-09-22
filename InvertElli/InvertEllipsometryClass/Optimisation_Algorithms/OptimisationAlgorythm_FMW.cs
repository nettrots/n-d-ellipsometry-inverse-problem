using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;

namespace InvertEllipsometryClass.Optimisation_Algorithms  
{
    public abstract class OptimisationAlgorythm_FMW
    {
        protected double psi = 0,
              delta = 0;

        protected Functional func;

        public abstract OptimizeResult Optimize();
    }
    public struct OptimizeResult
    {
        private Type objectType; 
        private object o;

        

        public object Pack
        {
            set { 
                o = value;
                OType = o.GetType();
            }
            get
            {
                return o;
            }
        }

        public Type OType
        {
            get { return objectType; }
            set { objectType = value; }
        }
    }
}
