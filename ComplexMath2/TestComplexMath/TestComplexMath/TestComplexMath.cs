//////////////////////////////////////
// Test ComplexMath.cs
//
// Complex.cs - Complex Class
// Copyright © PliaTech Software 2007
// Author: Michael B. Pliam
// Originated October 27, 2007
// Revised October 31, 2007
// Revised November 10, 2007


using System;
using System.Collections.Generic;
using System.Text;
using ComplexMath;

namespace TestComplexMath
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

            Complex z1, z2, z3;

            // instantiation
            z1 = new Complex(1.0, -2.0);
            z2 = new Complex(-3.0, -5.0);
            z3 = new Complex(0.292342, -0.394875);

            Console.WriteLine("instantiated complex numbers:");
            Console.WriteLine(string.Format("z1 = {0}", z1));
            Console.WriteLine(string.Format("z2 = {0}", z2));
            Console.WriteLine(string.Format("z3 = {0}", z3));

            // random complex numbers
            z1 = Complex.Random();
            z2 = Complex.Random();
            z3 = Complex.Random();

            Console.WriteLine("random generated:");
            Console.WriteLine(string.Format("z1 = {0}", z1));
            Console.WriteLine(string.Format("z2 = {0}", z2));
            Console.WriteLine(string.Format("z3 = {0}", z3));

            Console.WriteLine("basic properties:");
            Console.WriteLine(string.Format("z1.real = {0}  z1.imag = {1}", z1.real, z1.imag));
            Console.WriteLine(string.Format("complex conjugate of z1 = {0}.", z1.Conjugate));
            Console.WriteLine(string.Format("complex modulus of z1 = {0}.", z1.Modulus));
            Console.WriteLine(string.Format("complex argument of z1 = {0}.", z1.Argument));

            // operators
            Console.WriteLine("operators:");
            z3 = z1 + z2; Console.WriteLine(string.Format("z1 + z2 = {0}", z3));
            z3 = z1 - z2; Console.WriteLine(string.Format("z1 - z2 = {0}", z3));
            z3 = z1 * z2; Console.WriteLine(string.Format("z1 * z2 = {0}", z3));
            z3 = z1 / z2; Console.WriteLine(string.Format("z1 / z2 = {0}", z3));

            // transcendental functions
            Console.WriteLine("transcendental functions:");
            z3 = Complex.Sin(z1); Console.WriteLine(string.Format("Sin(z1) = {0}", z3));
            z3 = Complex.Cos(z1); Console.WriteLine(string.Format("Cos(z1) = {0}", z3));
            z3 = Complex.Tan(z1); Console.WriteLine(string.Format("Tan(z1) = {0}", z3));
            z3 = Complex.Sinh(z1); Console.WriteLine(string.Format("Sinh(z1) = {0}", z3));
            z3 = Complex.Cosh(z1); Console.WriteLine(string.Format("Cosh(z1) = {0}", z3));
            z3 = Complex.Tanh(z1); Console.WriteLine(string.Format("Tanh(z1) = {0}", z3));
            z3 = Complex.Exp(z1); Console.WriteLine(string.Format("Exp(z1) = {0}", z3));
            z3 = Complex.Log(z1); Console.WriteLine(string.Format("Log(z1) = {0}", z3));
            z3 = Complex.Log10(z1); Console.WriteLine(string.Format("Log10(z1) = {0}", z3));

            // roots
            Console.WriteLine("roots:");
            z3 = Complex.Sqrt(z1); Console.WriteLine(string.Format("Sqrt(z1) = {0}", z3));

            Console.WriteLine("k 4th roots:");
            for (long k = 0; k < 4; k++)
            {
                z3 = Complex.NthRoot(z1, 4, k);
                Console.WriteLine(string.Format("NthRoot(z1, 4, {0}) = {1}",k, z3));
            }

            Console.WriteLine("the five 5th roots of unity:");
            Complex z = new Complex(1.0, 0.0);
            for (long k = 0; k < 5; k++)
            {
                z3 = Complex.NthRoot(z, 5, k);
                Console.WriteLine(string.Format("{0}th 5th root of unity = {1}", k, z3));

            }

            // powers
            Console.WriteLine("powers:");
            int n = 7;
            z3 = Complex.Pow(z1, z2); Console.WriteLine(string.Format("Complex.Pow(z1, z2) = {0}", z3));
            z3 = Complex.Pow(z1, n); Console.WriteLine(string.Format("Complex.Pow(z1, {0}) = {1}", n, z3));

            // display
            Console.WriteLine("display:");
            z1.Show();
            z1.ShowPolar();
            Console.WriteLine(string.Format("polar form of z1 = {0}.", z1.Polar()));

            // debug tracing
            z1.Dump();


            Console.WriteLine("\nNormal Termination\n");

            //Console.ReadLine();

        }

    }
}
