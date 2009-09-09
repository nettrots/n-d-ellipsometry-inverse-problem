using System;
using System.Collections.Generic;
using System.Text;
using ComplexMath;

namespace TestDLL
{
    class Program
    {
        static void Main(string[] args)
        {
            Program.Test();
        }

        public static void Test()
        {
            Console.WriteLine("Normal Beginning\n");

            //static double srnd = new Random();

            Complex z1, z2, z3;

            // generate some random complex numbers
            z1 = Complex.Random();
            z2 = Complex.Random();
            z3 = Complex.Random();
            Console.WriteLine(string.Format("z1 = {0}", z1));
            Console.WriteLine(string.Format("z2 = {0}", z2));
            Console.WriteLine(string.Format("z3 = {0}", z3));
            /* 
            z1 = ( 0.0170055739660773,  0.475431244110424 )
            z2 = ( 0.562085224111604,  0.720434193369203 )
            z3 = ( 0.982858480877643,  0.213304637099292 )
             * 
            z1 = ( 0.38160862465464,  0.282567815986726 )
            z2 = ( 0.493835404745226,  0.842003089302221 )
            z3 = ( 0.592624074589752,  0.0907069067893116 )
             * 
            z1 = ( 0.0618562023443432,  0.0622684792905433 )
            z2 = ( 0.960061205066769,  0.147320054540094 )
            z3 = ( 0.885829960408541,  0.961100250464445 )
             */

            Console.WriteLine("\nNormal Termination\n");

        }

    }
}
