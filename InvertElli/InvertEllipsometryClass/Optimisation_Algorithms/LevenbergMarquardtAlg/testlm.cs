
using System;

namespace alglib
{
    public class testlm
    {
        public static bool testminlm(bool silent)
        {
            bool result = new bool();
            bool waserrors = new bool();
            bool referror = new bool();
            bool lin1error = new bool();
            bool lin2error = new bool();
            bool eqerror = new bool();
            bool converror = new bool();
            bool scerror = new bool();
            int rkind = 0;
            int nset = 0;
            int n = 0;
            int m = 0;
            double[] x = new double[0];
            double[] xe = new double[0];
            double[] b = new double[0];
            int i = 0;
            int j = 0;
            double v = 0;
            double[,] a = new double[0,0];
            minlm.lmstate state = new minlm.lmstate();
            minlm.lmreport rep = new minlm.lmreport();

            waserrors = false;
            referror = false;
            lin1error = false;
            lin2error = false;
            eqerror = false;
            converror = false;
            scerror = false;
            
            //
            // Reference problem.
            // RKind is a algorithm selector:
            // * 0 = FJ
            // * 1 = FGJ
            // * 2 = FGH
            //
            x = new double[2+1];
            n = 3;
            m = 3;
            for(rkind=0; rkind<=2; rkind++)
            {
                x[0] = 100*AP.Math.RandomReal()-50;
                x[1] = 100*AP.Math.RandomReal()-50;
                x[2] = 100*AP.Math.RandomReal()-50;
                if( rkind==0 )
                {
                    minlm.minlmfj(n, m, ref x, 0.0, 0.0, 0, ref state);
                }
                if( rkind==1 )
                {
                    minlm.minlmfgj(n, m, ref x, 0.0, 0.0, 0, ref state);
                }
                if( rkind==2 )
                {
                    minlm.minlmfgh(n, ref x, 0.0, 0.0, 0, ref state);
                }
                while( minlm.minlmiteration(ref state) )
                {
                    
                    //
                    // (x-2)^2 + y^2 + (z-x)^2
                    //
                    state.f = AP.Math.Sqr(state.x[0]-2)+AP.Math.Sqr(state.x[1])+AP.Math.Sqr(state.x[2]-state.x[0]);
                    if( state.needfg | state.needfgh )
                    {
                        state.g[0] = 2*(state.x[0]-2)+2*(state.x[0]-state.x[2]);
                        state.g[1] = 2*state.x[1];
                        state.g[2] = 2*(state.x[2]-state.x[0]);
                    }
                    if( state.needfij )
                    {
                        state.fi[0] = state.x[0]-2;
                        state.fi[1] = state.x[1];
                        state.fi[2] = state.x[2]-state.x[0];
                        state.j[0,0] = 1;
                        state.j[0,1] = 0;
                        state.j[0,2] = 0;
                        state.j[1,0] = 0;
                        state.j[1,1] = 1;
                        state.j[1,2] = 0;
                        state.j[2,0] = -1;
                        state.j[2,1] = 0;
                        state.j[2,2] = 1;
                    }
                    if( state.needfgh )
                    {
                        state.h[0,0] = 4;
                        state.h[0,1] = 0;
                        state.h[0,2] = -2;
                        state.h[1,0] = 0;
                        state.h[1,1] = 2;
                        state.h[1,2] = 0;
                        state.h[2,0] = -2;
                        state.h[2,1] = 0;
                        state.h[2,2] = 2;
                    }
                    scerror = scerror | !rkindvsstatecheck(rkind, ref state);
                }
                minlm.minlmresults(ref state, ref x, ref rep);
                referror = referror | rep.terminationtype<=0 | Math.Abs(x[0]-2)>0.001 | Math.Abs(x[1])>0.001 | Math.Abs(x[2]-2)>0.001;
            }
            
            //
            // end
            //
            waserrors = referror | lin1error | lin2error | eqerror | converror | scerror;
            if( !silent )
            {
                System.Console.Write("TESTING LEVENBERG-MARQUARDT OPTIMIZATION");
                System.Console.WriteLine();
                System.Console.Write("REFERENCE PROBLEM:                        ");
                if( referror )
                {
                    System.Console.Write("FAILED");
                    System.Console.WriteLine();
                }
                else
                {
                    System.Console.Write("OK");
                    System.Console.WriteLine();
                }
                System.Console.Write("1-D PROBLEM #1:                           ");
                if( lin1error )
                {
                    System.Console.Write("FAILED");
                    System.Console.WriteLine();
                }
                else
                {
                    System.Console.Write("OK");
                    System.Console.WriteLine();
                }
                System.Console.Write("1-D PROBLEM #2:                           ");
                if( lin2error )
                {
                    System.Console.Write("FAILED");
                    System.Console.WriteLine();
                }
                else
                {
                    System.Console.Write("OK");
                    System.Console.WriteLine();
                }
                System.Console.Write("LINEAR EQUATIONS:                         ");
                if( eqerror )
                {
                    System.Console.Write("FAILED");
                    System.Console.WriteLine();
                }
                else
                {
                    System.Console.Write("OK");
                    System.Console.WriteLine();
                }
                System.Console.Write("CONVERGENCE PROPERTIES:                   ");
                if( converror )
                {
                    System.Console.Write("FAILED");
                    System.Console.WriteLine();
                }
                else
                {
                    System.Console.Write("OK");
                    System.Console.WriteLine();
                }
                System.Console.Write("STATE FIELDS CONSISTENCY:                 ");
                if( scerror )
                {
                    System.Console.Write("FAILED");
                    System.Console.WriteLine();
                }
                else
                {
                    System.Console.Write("OK");
                    System.Console.WriteLine();
                }
                if( waserrors )
                {
                    System.Console.Write("TEST FAILED");
                    System.Console.WriteLine();
                }
                else
                {
                    System.Console.Write("TEST PASSED");
                    System.Console.WriteLine();
                }
                System.Console.WriteLine();
                System.Console.WriteLine();
            }
            result = !waserrors;
            return result;
        }


        /*************************************************************************
        Asserts that State fields are consistent with RKind.
        Returns False otherwise.
        *************************************************************************/
        private static bool rkindvsstatecheck(int rkind,
            ref minlm.lmstate state)
        {
            bool result = new bool();
            int nset = 0;

            nset = 0;
            if( state.needf )
            {
                nset = nset+1;
            }
            if( state.needfg )
            {
                nset = nset+1;
            }
            if( state.needfij )
            {
                nset = nset+1;
            }
            if( state.needfgh )
            {
                nset = nset+1;
            }
            if( nset!=1 )
            {
                result = false;
                return result;
            }
            if( rkind==0 & (state.needfg | state.needfgh) )
            {
                result = false;
                return result;
            }
            if( rkind==1 & state.needfgh )
            {
                result = false;
                return result;
            }
            if( rkind==2 & state.needfij )
            {
                result = false;
                return result;
            }
            result = true;
            return result;
        }


        /*************************************************************************
        Silent unit test
        *************************************************************************/
        public static bool testlm_test_silent()
        {
            bool result = new bool();

            result = testminlm(true);
            return result;
        }


        /*************************************************************************
        Unit test
        *************************************************************************/
        public static bool testlm_test()
        {
            bool result = new bool();

            result = testminlm(false);
            return result;
        }
    }
}
