//////////////////////////////////////
// Complex.cs - Complex Class
// Copyright © PliaTech Software 2007
// Author: Michael B. Pliam
// Originated October 27, 2007
// Revised October 31, 2007
// Revised November 09, 2007
// Revised November 14, 2007

using System;
using System.Collections.Generic;
using System.Text;
using System.Diagnostics;


namespace ComplexMath
{
    /// <summary>
    /// <remarks>
    /// Complex struct encapsulates fields double real and double imag.  
    /// Struct objects are designed to behave as complex numbers.
    /// </remarks>
    /// </summary>
    public struct Complex
    {
        #region fields

        /// <summary>
        /// <remarks>
        /// The real double component of the complex number.
        /// </remarks>
        /// </summary>
        private double m_real;

        /// <summary>
        /// <remarks>
        /// The imaginary double component of the complex number.
        /// </remarks>
        /// </summary>
        private double m_imag;

        /// <summary>
        /// <remarks>
        /// Static random number generator a Random class object which uses current time seed.
        /// </remarks>
        /// </summary>
        private static Random rnd = new Random();

        #endregion

        #region initializers

        /// <summary>
        /// <remarks>
        /// Complex struct initializer uses real and imaginary input parameters of 
        /// type double.
        /// </remarks>
        /// </summary>
        /// <param name="real"></param>
        /// <param name="imag"></param>
        public Complex(double real, double imag)
        {
            this.m_real = real;
            this.m_imag = imag;
        }

        /// <summary>
        /// <remarks>
        /// Complex struct creates a new Complex struct from an existing one.
        /// </remarks>
        /// <example>
        /// <code>
        ///     // Create a new Complex struct using an existing Complex struct.
        ///     Complex z1 = new Complex(0.123, -0.234);
        ///     Complex z2 = new Complex(z1);
        ///     Console.WriteLine(string.Format("z1 = {0}", z1));
        ///     Console.WriteLine(string.Format("z2 = {0}", z2));
        /// 
        ///     Output:
        ///     z1 = ( 0.123,  -0.234 )
        ///     z2 = ( 0.123,  -0.234 )
        /// </code>
        /// </example>
        /// </summary>
        /// <param name="z"></param>
        public Complex(Complex z)
        {
            this.m_real = z.real;
            this.m_imag = z.imag;
        }

        /// <summary>
        /// <remarks>
        /// Set an existing Complex object with new components.
        /// </remarks>
        /// <example>
        /// <code>
        ///     // Create a Complex object with no parameters.
        ///     Complex z = new Complex();
        ///     Console.WriteLine(string.Format("z = {0}", z));
        /// 
        ///     // Reset the real and imaginary components using the complex method.
        ///     z.complex(12.0, 24.3);
        ///     Console.WriteLine(string.Format("z = {0}", z));
        /// 
        ///     Output:
        ///     z = (0, 0)
        ///     z = (12.0, 24.3)
        /// </code>
        /// </example>
        /// </summary>
        /// <param name="re"></param>
        /// <param name="im"></param>
        /// <returns></returns>
        public void complex(double re, double im)
        {
            //Complex cx = new Complex(re, im);
            //return cx;
            this.m_real = re;
            this.m_imag = im;
        }

        /// <summary>
        /// <remarks>
        /// Create a random Complex number with real and imaginary components between 0 and 1.
        /// </remarks>
        /// <example>
        /// <code>
        ///     // Create a random complex number
        ///     Complex z = Complex.Random();
        ///     Console.WriteLine(string.Format("z = {0}", z));
        ///     
        ///     Output:  ( 0.400457408465658,  0.729048190977913 )
        /// </code>
        /// </example>
        /// </summary>
        /// <returns>Complex</returns>
        public static Complex Random()
        {
            Complex z = new Complex();

            z.real = rnd.NextDouble();
            z.imag = rnd.NextDouble();

            return z;
        }

        /// <summary>
        /// <remarks>
        /// Create a random Complex number with real and imaginary components within
        /// and including the selected range from lo to hi
        /// </remarks>
        /// <example>
        /// <code>
        ///     // Create a random complex number
        ///     Complex z = Complex.Random(1, 9);
        ///     Console.WriteLine(string.Format("z = {0}", z));
        ///     
        ///     Output:  z = ( 6.82405618849399,  7.08496072426669 )
        /// </code>
        /// </example>
        /// </summary>
        /// <returns>Complex</returns>
        public static Complex Random(double lo, double hi)
        {
            Complex z = new Complex();
            double range = hi - lo;
            z.real = lo + range * rnd.NextDouble();
            z.imag = lo + range * rnd.NextDouble();

            return z;
        }

        #endregion

        #region accessors

        /// <summary>
        /// <remarks>
        /// Access the real component property of the Complex number.
        /// This accessor is both read and write.
        /// </remarks>
        /// <example>
        /// <code>
        ///     // use the accessor to read the real component
        ///     // given a complex number z = (123.234, -345.456)
        ///     double re = z.real;
        ///     Console.WriteLine(string.Format("z.real = {0}", re));
        /// 
        ///     Output: z.real = 123.234
        /// 
        ///     // use the accessor to write the real component
        ///     // given the complex number z = (123.234, -345.456)
        ///     z.real = 678.345;
        ///     Console.WriteLine(string.Format("z.real = {0}", re));
        /// 
        ///     Output: z.real = 678.345
        /// </code>
        /// </example>
        /// </summary>
        public double real
        {
            get { return m_real; }
            set { m_real = value; }
        }

        /// <summary>
        /// <remarks>
        /// Access the imaginary component property of the Complex number.
        /// This accessor is both read and write.
        /// </remarks>
        /// <example>
        /// <code>
        ///     // use the accessor to read the imaginary component
        ///     // given a complex number z = (123.234, -345.456)
        ///     double im = z.imag;
        ///     Console.WriteLine(string.Format("z.imag = {0}", im));
        /// 
        ///     Output: z.imag = -345.456
        /// 
        ///     // use the accessor to write the imaginary component
        ///     // given the complex number z = (123.234, -345.456)
        ///     im = z.imag = 678.345;
        ///     Console.WriteLine(string.Format("z.imag = {0}", im));
        /// 
        ///     Output: z.imag = 678.345
        /// </code>
        /// </example>
        /// </summary>
        public double imag
        {
            get { return m_imag; }
            set { m_imag = value; }
        }

        #endregion

        #region operators

        /// <summary>
        /// <remarks>
        /// Complex addition operator facilitates addition of two complex numbers.
        /// Note that overriding this operator automatically overrides the += operator as well.
        /// </remarks>
        /// <example>
        /// <code>
        ///     // Add the Complex numbers A and B
        ///     Complex A = new Complex(1.0, 2.0);
        ///     Complex B = new Complex(4.0, 5.0);
        ///     Complex C = A + B;
        ///     Console.WriteLine(string.Format("A + B = {0}", C));
        ///     
        ///     Output: A + B = (5.0, 7.0)
        /// </code>
        /// </example>
        /// </summary>
        /// <param name="z1"></param>
        /// <param name="z2"></param>
        /// <returns>Complex</returns>
        public static Complex operator +(Complex z1, Complex z2)
        {
            return new Complex(z1.real + z2.real, z1.imag + z2.imag);
        }

        /// <summary>
        /// <remarks>
        /// Complex subtraction operator facilitates subtraction of two complex numbers.
        /// Note that overriding this operator automatically overrides the -= operator as well.
        /// </remarks>
        /// <example>
        /// <code>
        ///     // Subtract the Complex number B from A
        ///     Complex A = new Complex(1.0, 2.0);
        ///     Complex B = new Complex(4.0, 5.0);
        ///     Complex C = A - B;
        ///     Console.WriteLine(string.Format("A - B = {0}", C));
        ///     
        ///     Output: A + B = (-3.0, -3.0)
        /// </code>
        /// </example>
        /// </summary>
        /// <param name="z1"></param>
        /// <param name="z2"></param>
        /// <returns>Complex</returns>
        public static Complex operator -(Complex z1, Complex z2)
        {
            return new Complex(z1.real - z2.real, z1.imag - z2.imag);
        }

        /// <summary>
        /// <remarks>
        /// The unary negation operator placed before a complex number negates both real
        /// and imaginary components of that number.
        /// </remarks>
        /// <example>
        /// <code>
        ///     // use the unary negation operator to negate complex number
        ///     Complex z = Complex(32.45, -65.43);
        ///     Console.WriteLine(string.Format("z = {0}", z));
        ///     
        ///     z = -z;
        ///     Console.WriteLine(string.Format("z = {0}", z));
        /// 
        ///     Output: 
        ///     z = (32.45, -65.43)
        ///     z = (-32.45, 65.43)
        /// </code>
        /// </example>
        /// </summary>
        /// <param name="z"></param>
        /// <returns>Complex</returns>
        public static Complex operator -(Complex z)
        {
            z.real = -z.real;
            z.imag = -z.imag;
            return z;
        }

        /// <summary>
        /// <remarks>
        /// Complex multiplication operator facilitates multiplication of two complex numbers.
        /// Note that overriding this operator automatically overrides the *= operator as well.
        /// </remarks>
        /// <example>
        /// <code>
        ///     // Multiply the Complex numbers A and B
        ///     Complex A = new Complex(1.0, 2.0);
        ///     Complex B = new Complex(4.0, 5.0);
        ///     Complex C = A * B;
        ///     Console.WriteLine(string.Format("A * B = {0}", C));
        ///     
        ///     Output: A * B = (-6.0, 13.0)
        /// </code>
        /// </example>
        /// </summary>
        /// <param name="z1"></param>
        /// <param name="z2"></param>
        /// <returns>Complex</returns>
        public static Complex operator *(Complex z1, Complex z2)
        {
            return new Complex(
                z1.real * z2.real - z1.imag * z2.imag,
                z1.real * z2.imag + z1.imag * z2.real);
        }

        /// <summary>
        /// <remarks>
        /// This version of the * operator allows you to multiply a complex number
        /// by a double value.  The result is a complex number.
        /// </remarks>
        /// <example>
        /// <code>
        ///     // multiply a complex number by a scalar double
        ///     Complex z = new Complex(1.0, 2.0);
        ///     double x = 3.0;
        ///     Complex z2 = z * x;
        ///     Console.WriteLine(string.Format("{0} * {1} = {2}", z1, x, z2));
        /// 
        ///     Output: ( 3,  6 ) * 3 = ( 4,  5 )
        /// </code>
        /// </example>
        /// </summary>
        /// <param name="z"></param>
        /// <param name="dval"></param>
        /// <returns>Complex</returns>
        static public Complex operator *(Complex z, double dval)
        {
            Complex z2 = new Complex();
            z2.real = (z.real * dval);
            z2.imag = (z.imag * dval);
            return z2;
        }

        /// <summary>
        /// <remarks>
        /// This version of the * operator allows you to multiply a double value
        /// by a complex number.  The result is a complex number.
        /// </remarks>
        /// <example>
        /// <code>
        ///     // multiply a scalar double value by a complex number
        ///     Complex z = new Complex(1.0, 2.0);
        ///     double x = 3.0;
        ///     Complex z2 = x * z;
        ///     Console.WriteLine(string.Format("{0} * {1} = {2}", x, z1, z2));
        /// 
        ///     Output: 3 * ( 3,  6 ) = ( 4,  5 )
        /// </code>
        /// </example>
        /// </summary>
        /// <param name="dval"></param>
        /// <param name="z"></param>
        /// <returns>Complex</returns>
        static public Complex operator *(double dval, Complex z)
        {
            Complex z2 = new Complex();
            z2.real = (z.real * dval);
            z2.imag = (z.imag * dval);
            return z2;
        }

        /// <summary>
        /// <remarks>
        /// Complex division operator facilitates division of one complex number by another.
        /// Note that overriding this operator automatically overrides the /= operator as well.
        /// </remarks>
        /// <example>
        /// <code>
        ///     // Calculate the Complex quotient of A / B where
        ///     Complex A = new Complex(1.0, 2.0);
        ///     Complex B = new Complex(4.0, 5.0);
        ///     Complex C = A / B;
        ///     Console.WriteLine(string.Format("A / B = {0}", C));
        ///     
        ///     Output: A / B = ( 0.341463414634146,  0.0731707317073171 )
        /// </code>
        /// </example>
        /// </summary>
        /// <param name="z1"></param>
        /// <param name="z2"></param>
        /// <returns>Complex</returns>
        public static Complex operator /(Complex z1, Complex z2)
        {
            double value = z2.real * z2.real + z2.imag * z2.imag;

            return new Complex(
                (z1.real * z2.real + z1.imag * z2.imag) / value,
                (z1.imag * z2.real - z1.real * z2.imag) / value);
        }

        /// <summary>
        /// <remarks>
        /// This version of the / operator allows you to divide a complex number
        /// by a real scalar double value.  The result is a complex number.
        /// </remarks>
        /// <example>
        /// <code>
        ///     // Divide a complex number by a scalar double
        ///     Complex z1 = new Complex(1.0, 2.0);
        ///     double x = 3.125;
        ///     Complex z2 = z1 / x;
        ///     Console.WriteLine(string.Format("{0} / {1} = {2}", z1, x, z2));
        /// 
        ///     Output: ( 1,  2 ) / 3.125 = ( 0.32,  0.64 )
        /// </code>
        /// </example>
        /// </summary>
        /// <param name="z"></param>
        /// <param name="dval"></param>
        /// <returns>Complex</returns>
        static public Complex operator /(Complex z, double dval)
        {
            Complex z2 = new Complex();
            z2.real = z.real / dval;
            z2.imag = z.imag / dval;
            return z2;
        }

        /// <summary>
        /// <remarks>
        /// The comparison operator '==' facilitates the comparison of two distinct complex numbers.  The method
        /// returns true if complex numbers are identical, false if they are not.  The comparison is made between 
        /// the two real doubles and the two imaginary doubles respectively.
        /// </remarks>
        /// <example>
        /// <code>
        ///     // Compare two complex numbers, z1 and z2
        ///     Complex z1 = new Complex(0.2134345345, -0.2938749238);
        ///     Complex z2 = new Complex(0.2343454564, -0.2349893454);
        ///     if(z1 == z2) Console.WriteLine("z1 and z2 are identical"); 
        ///         else Console.WriteLine("z1 is not identical to z2");
        /// 
        ///     z1 = z1.complex(0.2134345345, -0.2938749238);
        ///     z2 = z2.complex(0.2134345345, -0.2938749238);
        ///     if(z1 == z2) Console.WriteLine("z1 and z2 are identical"); 
        ///         else Console.WriteLine("z1 is not identical to z2");
        /// 
        ///     Output: 
        ///     z1 is not identical to z2
        ///     z1 and z2 are identical
        /// </code>
        /// </example>
        /// </summary>
        /// <param name="z1"></param>
        /// <param name="z2"></param>
        /// <returns>bool</returns>
        public static bool operator ==(Complex z1, Complex z2)
        {
            return (z1.real == z2.real && z1.imag == z2.imag);
        }

        /// <summary>
        /// <remarks>
        /// The comparison operator '!=' facilitates the comparison of two distinct complex numbers.  
        /// The method returns true if complex numbers are not identical, false if they are identical.  
        /// The comparison is made between the two real doubles and the two imaginary doubles respectively.  
        /// </remarks>
        /// <example>
        /// <code>
        ///     // Compare two complex numbers, z1 and z2
        ///     Complex z1 = new Complex(0.2134345345, -0.2938749238);
        ///     Complex z2 = new Complex(0.2343454564, -0.2349893454);
        ///     if(z1 != z2) Console.WriteLine("z1 is not identical to z2");
        ///         else Console.WriteLine("z1 and z2 are identical"); 
        /// 
        ///     z1 = z1.complex(0.2134345345, -0.2938749238);
        ///     z2 = z2.complex(0.2134345345, -0.2938749238);
        ///     if(z1 != z2) Console.WriteLine("z1 is not identical to z2");
        ///         else Console.WriteLine("z1 and z2 are identical"); 
        /// 
        ///     Output: 
        ///     z1 is not identical to z2
        ///     z1 and z2 are identical
        /// </code>
        /// </example>
        /// </summary>
        /// <param name="z1"></param>
        /// <param name="z2"></param>
        /// <returns>bool</returns>
        public static bool operator !=(Complex z1, Complex z2)
        {
            return (z1.real != z2.real || z1.imag != z2.imag);
        }

        #endregion

        #region algebra

        /// <summary>
        /// <remarks>
        /// The complex conjugate of a complex number is obtained by changing the sign of the imaginary part. 
        /// Thus, the conjugate of the complex number z = (a, b) (where a and b are real numbers, a is the real 
        /// part and b the imaginary part) is (a, -b).  In polar form, the complex conjugate of z = r*e^(i*t) is
        /// r*e^(-i*t).
        /// </remarks>
        /// <example>
        /// <code>
        ///     // obtain the complex conjugate of a complex number
        ///     Complex z1 = new Complex(3, 4);
        ///     Complex z2 = z1.Conjugate;
        ///     Console.WriteLine(string.Format("The complex conjugate of {0} is {1}.", z1, z2));
        /// 
        ///     Output:
        ///     The complex conjugate of (3, 4) is (3, -4).
        /// </code>
        /// </example>
        /// </summary>
        public Complex Conjugate
        {
            get { return new Complex(m_real, -m_imag); }
        }

        /// <summary>
        /// <remarks>
        /// The complex norm of a complex number z = (a, ib), also called the  modulus, is denoted |z| 
        /// and defined by sqrt(real*real + imag*imag). If z is expressed as a complex exponential 
        /// (i.e., a phasor), then the norm can be expressed as |r*e^(i*t)| = |r|.  For some reason,
        /// some authors refer to the complex norm as simply real*real + imag*imag (without taking the
        /// square root of the sum). This makes little sense from a geometric point of view and is
        /// probably an error.
        /// </remarks>
        /// <example>
        /// <code>
        ///     // find the norm of a complex number
        ///     Complex z = new Complex(1.2345, -2.3456);
        ///     double nrm = z.Norm;
        ///     Console.WriteLine(string.Format("The norm of {0} is {1}.", z, nrm));
        ///     double mod = z.Modulus;
        ///     Console.WriteLine(string.Format("The modulus of {0} is {1}", z, mod));
        ///
        ///     Output:
        ///     The norm of ( 1.2345,  -2.3456 ) is 2.65062815385335.
        ///     The modulus of ( 1.2345,  -2.3456 ) is 2.65062815385335.
        /// </code>
        /// </example>
        /// </summary>
        public double Norm
        {
            get { return Math.Sqrt(m_real * m_real + m_imag * m_imag); }
        }

        /// <summary>
        /// <remarks>
        /// The square |z|^2 of |z| is sometimes called the absolute square. The absolute square of 
        /// a complex number z, also known as the squared norm, is defined as |z|^2 = z * z_, where
        /// z_ is the complex conjugate of z and |z| is the complex modulus. If the complex number is 
        /// written z = (a, b), with a and b real, then the absolute square can be written |z| = a^2 + b^2.
        /// Note that even though the absolute square is the product of two complex numbers, the result is
        /// a real number of type double because the imaginary part is always zero.
        /// </remarks>
        /// <example>
        /// <code>
        ///     // find the absolute square (squared norm) of a complex number
        ///     Complex z = new Complex(1.2345, -2.3456);
        ///     double absq = z.AbsoluteSquare;
        ///     Console.WriteLine(string.Format("The absolute square (squared norm) of {0} is {1}.", z, absq));
        /// 
        ///     Output:
        ///     The absolute square (squared norm) of ( 1.2345,  -2.3456 ) is 7.02582961.
        /// </code>
        /// </example>
        /// </summary>
        public double AbsoluteSquare
        {
            get { return (this * this.Conjugate).real; }
            
        }

        /// <summary>
        /// <remarks>
        /// Provides the modulus (norm, absolute value) of a complex number.  Recall that 
        /// any complex number z = |z| * e^(i * t) where t = argument and |z| = modulus.
        /// The modulus is the value of the radial vector in the Argand Plane that extends
        /// from the origin to a point representing the complex number [x(real), y(imag)].
        /// </remarks>
        /// <example>
        /// <code>
        ///     // obtain the modulus of a complex number
        ///     z1 = Complex.Polar(0.2938, -0.5434);
        ///     Console.WriteLine(string.Format("{0} = {1}", z1.Polar(), z1.ToString()));
        ///     Console.WriteLine(string.Format("z1.Modulus = {0}", z1.Modulus));
        ///
        ///     Output:
        ///     0.2938 * e^(i * -0.5434) = ( 0.251479774320806,  -0.151909061966548 )
        ///     z1.Modulus = 0.2938
        /// </code>
        /// </example>
        /// </summary>
        public double Modulus
        {
            get { return System.Math.Sqrt(m_real * m_real + m_imag * m_imag); }
        }

        /// <summary>
        /// <remarks>
        /// Provides the argument or angle t in radians of a complex number.  Recall that 
        /// any complex number z = |z| * e^(i * t) where t = argument and |z| = modulus.
        /// Recall that a complex number has arguments between 0 and 2PI.  This algorithm
        /// automatically compensates for angles greater than 2PI by using Math.Atan2.
        /// </remarks>
        /// <example>
        /// <code>
        ///     // obtain the argument t of a complex number
        ///     z1 = Complex.Polar(0.2938, -0.5434);
        ///     Console.WriteLine(string.Format("{0} = {1}", z1.Polar(), z1.ToString()));
        ///     Console.WriteLine(string.Format("z1.Argument = {0}", z1.Argument));
        ///
        ///     Output:
        ///     0.2938 * e^(i * -0.5434) = ( 0.251479774320806,  -0.151909061966548 )
        ///     z1.Argument = -0.5434
        /// </code>
        /// </example>
        /// </summary>
        public double Argument
        {
            get { return System.Math.Atan2(m_imag, m_real); }
        }

        /// <summary>
        /// <remarks>
        /// Create a new complex number by inputing the modulus and argument.
        /// </remarks>
        /// <example>
        /// <code>
        ///     // create a new complex number using a modulus and argument
        ///     z1 = Complex.Polar(0.2938, -0.5434);
        ///     Console.WriteLine(string.Format("{0} = {1}", z1.Polar(), z1.ToString()));
        /// 
        ///     Output:
        ///     0.2938 * e^(i * -0.5434) = ( 0.251479774320806,  -0.151909061966548 )
        /// </code>
        /// </example>
        /// </summary>
        /// <param name="modulus"></param>
        /// <param name="argument"></param>
        /// <returns>Complex</returns>
        public static Complex Polar(double modulus, double argument)
        {
            return new Complex(
               modulus * System.Math.Cos(argument),
               modulus * System.Math.Sin(argument));
        }

        /// <summary>
        /// <remarks>
        /// Calculates the complex value of a complex number z raised to a complex power n.
        /// </remarks>
        /// <example>
        /// <code>
        ///     // calculate the value of a complex number raised to a complex power
        ///     Complex a, b, c;
        ///     a = new Complex(3.2, -4.1);
        ///     b = new Complex(5.2, 1.8);
        ///     c = Complex.Pow(a, b);
        ///     Console.WriteLine(string.Format("Pow( {0}, {1} ) = {2}", a, b, c));
        ///     
        ///     Output:
        ///     Pow( ( 3.2,  -4.1 ), ( 5.2,  1.8 ) ) = ( -4943.88551718946,  -26678.8012334349 )
        /// </code>
        /// </example>
        /// </summary>
        /// <param name="baseNumber"></param>
        /// <param name="index"></param>
        /// <returns>Complex</returns>
        public static Complex Pow(Complex baseNumber, Complex index)
        {
            return Exp(index * Log(baseNumber));
        }

        /// <summary>
        /// <remarks>
        /// Calculates the complex value of a complex number z raised to an integral power n.
        /// </remarks>
        /// <example>
        /// <code>
        ///     // calculate the value of a complex number raised to an integral power
        ///     long n = 4;
        ///     Complex a, b;
        ///     a = new Complex(3.2, -4.1);
        ///     b = Complex.Pow(b, n);
        ///     Console.WriteLine(string.Format("Pow( {0}, {1} ) = {2}", a, n, b));
        ///     
        ///     Output:
        ///     Pow( ( 3.2,  -4.1 ), 4 ) = ( -5,  -12 )
        /// </code>
        /// </example>
        /// </summary>
        /// <param name="baseNumber"></param>
        /// <param name="index"></param>
        /// <returns>Complex</returns>
        public static Complex Pow(Complex baseNumber, long index)
        {
            Complex cindex = new Complex(index, 0.0);
            return Exp(cindex * Log(baseNumber));
        }

        /// <summary>
        /// <remarks>
        /// Find the complex square root of a complex number.
        /// </remarks>
        /// <example>
        /// <code>
        ///     // find the square root of complex number
        ///     Complex a = new Complex(2.0, -3.0);
        ///     b = Complex.Sqrt(a);
        ///     Console.WriteLine(string.Format("a = {0}", a));
        ///     Console.WriteLine(string.Format("Complex.Sqrt(a) = {0}", b));
        ///
        ///     Output:
        ///     a = ( 2,  -3 )
        ///     Complex.Sqrt(a) = ( 1.67414922803554,  -0.895977476129838 )
        /// </code>
        /// </example>
        /// </summary>
        /// <param name="z"></param>
        /// <returns>Complex</returns>
        public static Complex Sqrt(Complex z)
        {
            double value = System.Math.Sqrt(z.real * z.real + z.imag * z.imag) + z.real;

            return new Complex(
               System.Math.Sqrt(0.5 * value),
               System.Math.Sqrt(0.5 / value) * z.imag);
        }

        /// <summary>
        /// <remarks>
        /// The reciprocal of a complex number, z = (a, b), sometimes referred to as the 
        /// 'multiplicative inverse', can be computed as follows:
        ///     z^(-1) = Complex(1.0, 0.0)/ z.
        /// Using the Recip method is essentially the same as using the division operator
        /// on complex(1.0, 0.0).
        /// </remarks>
        /// <example>
        /// <code>
        ///     // use both Recip and the division operator to get a complex reciprocal.
        ///     a = a.complex(2.0, -3.0);
        ///     b = Complex.Recip(a);
        ///     Complex cone = new Complex(1.0, 0.0);
        ///     c = cone / a;
        ///     Console.WriteLine(string.Format("a = {0}", a));
        ///     Console.WriteLine(string.Format("recip1(a) = {0}", b));
        ///     Console.WriteLine(string.Format("recip2(a) = {0}", c));
        /// 
        ///     Output:
        ///     a = ( 2,  -3 )
        ///     recip1(a) = ( 0.153846153846154,  0.230769230769231 )
        ///     recip2(a) = ( 0.153846153846154,  0.230769230769231 )
        /// </code>
        /// </example>
        /// </summary>
        /// <param name="z"></param>
        /// <returns>Complex</returns>
        public static Complex Recip(Complex z)
        {
            if (z.real == 0.0 && z.imag == 0.0)
            {
                Console.WriteLine("multiplicative inverse (reciprocal) not defined for complex zero");
                return z;
            }

            Complex zr = new Complex(1.0, 0.0);
            zr.real = z.real / (z.real * z.real + z.imag * z.imag);
            zr.imag = -z.imag / (z.real * z.real + z.imag * z.imag);

            return zr;

        }// Recip(ComplexNumbers z)

        /// <summary>
        /// <remarks>
        /// According to Euler's formula, for any real number x,
        ///     e^(ix) = cos(x) + i sin(x)
        /// where
        ///     e is the base of the natural logarithm 
        ///     i is the imaginary unit 
        ///     cos and sin are trigonometric functions
        ///
        /// Complex numbers z = (a, ib) can be expressed in polar form:
        ///     z = r * e^(i*t) or mod * e^(i*arg) or z = |z|* e^(i*t)
        /// Using the Euler relationship, The polar form can also be expressed
        /// in terms of trigonometric functions:
        ///     z = r*e^(i*t) = r * (cos(t) + i sin(t))
        ///       = |z|* e^(i*t) = |z| * (cos(t) + i sin(t))
        /// The SetPolar function resets the current complex number 
        /// with Modulus mod and Argument arg.
        /// </remarks>
        /// <example>
        /// <code>
        ///    // use the SetPolar method to reset an existing complex number
        ///    a = a.complex(1.0, -2.0);
        ///    Console.WriteLine(string.Format("a = {0}", a));
        ///    a.ShowPolar();
        ///    a = a.SetPolar(3.60555127546399, 0.982793723247329);
        ///    Console.WriteLine(string.Format("a = {0}", a));
        ///     
        ///     Output:
        ///     a = ( 1,  -2 )
        ///     ( 2.23606797749979, -1.10714871779409 )
        ///     a = ( 2,  3 )
        /// </code>
        /// </example>
        /// </summary>
        /// <param name="mod"></param>
        /// <param name="arg"></param>
        /// <returns></returns>
        public Complex SetPolar(double mod, double arg)
        {
            Complex z = new Complex();
            double t = arg;
            z.real = mod * Math.Cos(t);
            z.imag = mod * Math.Sin(t);

            return z;
        }

        /// <summary>
        /// <remarks>
        /// De Moivre's formula states that for any complex number x 
        /// and any integer n it holds that
        ///     cos(x)+ i*sin(x))^n = cos(n*x) + i*sin(n*x)
        /// Since z = r*e^(i*t) = r * (cos(t) + i sin(t))
        ///     where
        ///     z = (a, ib)
        ///     r = modulus of z
        ///     t = argument of z
        ///     i = sqrt(-1.0)
        /// thence, one can calculate the nth root of z by the formula:  
        ///    z^(1/n) = r^(1/n) * (cos(x/n) + i sin(x/n))
        /// by using log division which makes this relatively simple.
        /// Note that there exists a multiplicity of complex roots which must be n - 1 in number.
        /// The k nth roots of z = r cis t are: ( where 'cis' stands for  cos(t) + i sin(t) )
        /// z^(1/n) = r^(1/n) cis ( t/n + k * 2 * PI /n) for k = 0, 1, 2, 3, ... n-1.
        /// This method returns the kth Nth root where k = 0, 1, 2, 3,...,n - 1.
        /// Notice that when k > n-1 the roots simply recycle since there are only n-1.
        /// unique roots. Do not confuse the "nth root" with the "nth root of unity".  Recall 
        /// that an "nth root of unity" is just another name for an nth root of one (or minus 
        /// one).  Thus the 5th root of unity is: b = Complex.NthRoots(a, 5, k) where a = (1.0, 0.0).
        /// Similar calculations can be made for a = (0.0, 1.0).
        /// </remarks>
        /// <example>
        /// <code>
        ///     // loop through to find the n-1 roots of a complex number
        ///     Complex a = new Complex(2.3, -4.5); Console.WriteLine(string.Format("a = {0}", a));
        ///     long n = 5; k = 0;
        ///     for (long i = 0; i -lt n; i++)
        ///     {
        ///         b = Complex.NthRoots(a, n, k);
        ///         c = Complex.Pow(b, n);  // check if root
        ///         Console.WriteLine(string.Format("NthRoots(a, {0}, {1} ) = {2} : {3}", n, k, b, c.Round(2)));
        ///         k++;
        ///     }
        ///
        ///     Output:
        ///     a = ( 2.3,  -4.5 )
        ///     NthRoots(a, 5, 0 ) = ( 1.34945770488324,  -0.301283056467905 ) : ( 2.3,  -4.5 )
        ///     NthRoots(a, 5, 1 ) = ( 0.703542578102255,  1.19030895912809 ) : ( 2.3,  -4.5 )
        ///     NthRoots(a, 5, 2 ) = ( -0.914644479083315,  1.03693445032258 ) : ( 2.3,  -4.5 )
        ///     NthRoots(a, 5, 3 ) = ( -1.26882395379819,  -0.549448224723052 ) : ( 2.3,  -4.5 )
        ///     NthRoots(a, 5, 4 ) = ( 0.130468149896011,  -1.37651212825971 ) : ( 2.3,  -4.5 )
        /// </code>
        /// </example>
        /// </summary>
        /// <param name="z"></param>
        /// <param name="n"></param>
        /// <param name="k"></param>
        /// <returns>Complex</returns>
        public static Complex NthRoot(Complex z, long n, long k)
        {
            double mod, arg, t;
            mod = z.Modulus;
            arg = z.Argument;

            double lmod = Math.Log(mod);
            double rlmod = lmod / n;
            double rmod = Math.Exp(rlmod);

            Complex zrn = new Complex();
            t = arg / n + 2 * k * Math.PI / n;
            zrn.real = rmod * Math.Cos(t);
            zrn.imag = rmod * Math.Sin(t);

            return zrn;

        }

        #endregion

        #region transcendentals

        /// <summary>
        /// <remarks>
        /// Calculates the complex cosine of the complex number z through use of the formula 
        /// cos(z) = ( exp(i*z) + exp(-i*z) ) / 2.
        /// </remarks>
        /// <example>
        /// <code>
        ///     // calculate the complex cosine of the complex number z
        ///     Complex z = new Complex(1.2345, -2.3456);
        ///     Complex cosz = Complex.Cos(z);
        ///     Console.WriteLine(string.Format("The complex cosine of {0} is {1}.", z, cosz));
        /// 
        ///     Output:
        ///     The complex cosine of ( 1.2345,  -2.3456 ) is ( 1.73829246326284,  4.88216132677658 ).
        /// </code>
        /// </example>
        /// </summary>
        /// <param name="z"></param>
        /// <returns>Complex</returns>
        public static Complex Cos(Complex z)
        {
            Complex z1 = Exp(new Complex(-z.imag, z.real));
            Complex z2 = Exp(new Complex(z.imag, -z.real));

            return new Complex(0.5 * (z1.real + z2.real), 0.5 * (z1.imag + z2.imag));
        }

        /// <summary>
        /// <remarks>
        /// Calculates the complex hyperbolic cosine of the complex number z through use of the formula 
        /// cosh(z) = ( exp(z) + exp(-z) ) / 2.
        /// </remarks>
        /// <example>
        /// <code>
        ///     // calculate the hyperbolic cosine of a complex number
        ///     Complex z = new Complex(1.2345, -2.3456);
        ///     Complex coshz = Complex.Cosh(z);
        ///     Console.WriteLine(string.Format("The complex hyperbolic cosine of {0} is {1}.", z, coshz));
        /// 
        ///     Output:
        ///     The complex hyperbolic cosine of ( 1.2345,  -2.3456 ) is ( -1.3038833630954,  -1.12388586077941 ).
        /// </code>
        /// </example>
        /// </summary>
        /// <param name="z"></param>
        /// <returns>Complex</returns>
        public static Complex Cosh(Complex z)
        {
            Complex z1 = Exp(z);
            Complex z2 = Exp(new Complex(-z.real, -z.imag));

            return new Complex(0.5 * (z1.real + z2.real), 0.5 * (z1.imag + z2.imag));
        }

        /// <summary>
        /// <remarks>
        /// Calculates the complex sine of the complex number z through use of the formula 
        /// sin(z) = ( exp(i*z) - exp(-i*z) ) / (2*i).
        /// </remarks>
        /// <example>
        /// <code>
        ///     // calculate the sine of a complex number
        ///     Complex z = new Complex(1.2345, -2.3456);
        ///     Complex sinz = Complex.Sin(z);
        ///     Console.WriteLine(string.Format("The complex sine of {0} is {1}.", z, sinz));
        /// 
        ///     Output:
        ///     The complex sine of ( 1.2345,  -2.3456 ) is ( 4.97258521662104,  -1.706682514036 ).
        /// </code>
        /// </example>
        /// </summary>
        /// <param name="z"></param>
        /// <returns>Complex</returns>
        public static Complex Sin(Complex z)
        {
            Complex z1 = Exp(new Complex(-z.imag, z.real));
            Complex z2 = Exp(new Complex(z.imag, -z.real));

            return new Complex(0.5 * (z1.imag - z2.imag), 0.5 * (z2.real - z1.real));
        }

        /// <summary>
        /// <remarks>
        /// Calculates the complex hyperbolic sine of the complex number z through use of the formula 
        /// sinh(z) = ( exp(z) - exp(-z) ) / 2.
        /// </remarks>
        /// <example>
        /// <code>
        ///     // calculate the hyperbolic sine of a complex number
        ///     Complex z = new Complex(1.2345, -2.3456);
        ///     Complex sinhz = Complex.Sinh(z);
        ///     Console.WriteLine(string.Format("The complex hyperbolic sine of {0} is {1}.", z, sinhz));
        /// 
        ///     Output:
        ///     The complex hyperbolic sine of ( 1.2345,  -2.3456 ) is ( -1.10032064508213,  -1.33180821648497 ).
        /// </code>
        /// </example>
        /// </summary>
        /// <param name="z"></param>
        /// <returns>Complex</returns>
        public static Complex Sinh(Complex z)
        {
            Complex z1 = Exp(z);
            Complex z2 = Exp(new Complex(-z.real, -z.imag));

            return new Complex(0.5 * (z1.real - z2.real), 0.5 * (z1.imag - z2.imag));
        }
        
        /// <summary>
        /// <remarks>
        /// Calculates the complex tangent of the complex number z through use of the formula 
        /// tan(z) = sin(z)/cos(z).
        /// </remarks>
        /// <example>
        /// <code>
        ///     // calculate the complex tangent of a complex number
        ///     Complex z = new Complex(1.2345, -2.3456);
        ///     Complex tanz = Complex.Tan(z);
        ///     Console.WriteLine(string.Format("The complex tangent of {0} is {1}.", z, tanz));
        /// 
        ///     Output:
        ///     The complex tangent of ( 1.2345,  -2.3456 ) is ( 0.0115986961735915,  -1.01439156942989 ).
        /// </code>
        /// </example>
        /// </summary>
        /// <param name="z"></param>
        /// <returns></returns>
        public static Complex Tan(Complex z)
        {
            return Sin(z) / Cos(z);
        }

        /// <summary>
        /// <remarks>
        /// Calculates the complex hyperbolic tangent of the complex number z through use of the formula 
        /// tanh(z) = sinh(z)/cosh(z) = (exp(2*z) - 1) / (exp(2*z) + 1) = -i * tan(i*z).
        /// </remarks>
        /// <example>
        /// <code>
        ///     // calculate the complex hyperbolic tangent of a complex number
        ///     Complex z = new Complex(1.2345, -2.3456);
        ///     Complex tanhz = Complex.Tanh(z);
        ///     Console.WriteLine(string.Format("The complex hyperbolic tangent of {0} is {1}.", z, tanhz));
        /// 
        ///     Output:
        ///     The complex hyperbolic tangent of ( 1.2345,  -2.3456 ) is ( 0.989288367008642,0.168696844208687 ).
        /// </code>
        /// </example>
        /// </summary>
        /// <param name="z"></param>
        /// <returns>Complex</returns>
        public static Complex Tanh(Complex z)
        {
            return Sinh(z) / Cosh(z);
        }
        
        /// <summary>
        /// <remarks>
        /// Calculates the complex exponential of the complex number z = (a, b) through use of the formula 
        /// exp(z) = exp(a) * (cos(b) + i*sin(b)).  In power form, it can be expressed e^(z) where e is
        /// Euler's number which is 2.7182818284590452... 
        /// </remarks>
        /// <example>
        /// <code>
        ///     // calculate the complex exponential of a complex number
        ///     Complex z = new Complex(1.2345, -2.3456);
        ///     Complex expz = Complex.Exp(z);
        ///     Console.WriteLine(string.Format("The complex exponential of {0} is {1}.", z, expz));
        /// 
        ///     Output:
        ///     The complex exponential of ( 1.2345,  -2.3456 ) is ( -2.40420400817753,  -2.45569407726438 ).
        /// </code>
        /// </example>
        /// </summary>
        /// <param name="z"></param>
        /// <returns>Complex</returns>
        public static Complex Exp(Complex z)
        {
            double value = System.Math.Exp(z.real);

            return new Complex(
               value * Math.Cos(z.imag),
               value * Math.Sin(z.imag));
        }

        /// <summary>
        /// <remarks>
        /// If the non-zero complex number z is expressed in polar coordinates as z = r * e^(i*t)
        /// with r > 0 and t is between -pi and pi, then log(z) = ln(r) + i*t, where ln(r) is the 
        /// usual natural logarithm of a real number. Recall that r is the complex modulus and
        /// t is the complex argument of z. So defined, log (z) is holomorphic for all complex 
        /// numbers which are not real numbers less than or equal to 0, and it has the property
        /// e^(log(z)) = z.  Recall that some familiar properties of the real-valued natural 
        /// logarithm are no longer valid for this complex extension. For example, log(e^(z)) does 
        /// not always equal z, and log(z*w) does not always equal log(z) + log(w) – in either case, 
        /// the result might have to be adjusted modulo 2*pi*i to stay within the range of this 
        /// principal branch of the complex log function.  A somewhat more natural definition of 
        /// log(z) interprets it as a multi-valued function. Thus for z = r*e^(i*t), it would be 
        /// possible to choose log(z) = ln(r) + i*(t + 2*pi*k) for any integer k.
        /// </remarks>
        /// <example>
        /// <code>
        ///     // calculate the complex exponential of a complex number
        ///     Complex z = new Complex(1.2345, -2.3456);
        ///     Complex logz = Complex.Log(z);
        ///     Console.WriteLine(string.Format("The complex natural logarithm of {0} is {1}.", z, logz));
        /// 
        ///     Output:
        ///     The complex natural logarithm of ( 1.2345,  -2.3456 ) is ( 0.974796651098723,  -1.08632718334389 ).
        /// </code>
        /// </example>
        /// </summary>
        /// <param name="z"></param>
        /// <returns>Complex</returns>
        public static Complex Log(Complex z)
        {
            return new Complex(Math.Log(z.Modulus), z.Argument);
        }

        /// <summary>
        /// <remarks>
        /// Given any real number greater than zero, the logarithm of that number to base (a) 
        /// is related to the logarithm of that number base (b) by the scaling formula
        /// log[b](x) = log[a](x)/log(b).  Consider the logarithm of a complex number z is
        /// log(z) = ln(r) + i*t, where ln(r) is the usual natural logarithm of a real number, 
        /// r is the complex modulus, and t is the complex argument of z.  Applying the above
        /// logarithmic scaling factor to the natural complex logarithm, we can obtain the 
        /// complex logarithm to virtually any desired base.  Using ln(10) = 2.30258... as
        /// the scaling factor for ln(z).
        /// </remarks>
        /// <example>
        /// <code>
        ///     // find log10 of a complex number z
        ///     Complex z = new Complex(0.12345, -0.23456);
        ///     Complex log10z = Complex.Log10(z);
        ///     Console.WriteLine(string.Format("The complex logarithm (base 10) of {0} is {1}.", z, log10z));
        /// 
        ///     Output:
        ///     The complex log10 of ( 1.2345,  -2.3456 ) is ( 0.423348806549945, -0.471785901267753 ).
        /// </code>
        /// </example>
        /// </summary>
        /// <param name="z"></param>
        /// <returns>Complex</returns>
        public static Complex Log10(Complex z)
        {
            const double log10 = 2.3025850929940459;

            Complex value = Log(z);

            value.real /= log10;
            value.imag /= log10;

            return value;
        }

        #endregion

        #region display

        /// <summary>
        /// <remarks>
        /// Provides a string representation of the complex number z = (x, y), where x and y are real
        /// numbers.  The representation of a complex number in this way is cleaner and less confusing
        /// than the more traditional z = (x + iy).  Provided the imaginary part always appears on the
        /// right, there is no need for the additional 'i'.  And since x and y are really components of
        /// a two element vector, the addition sign is superfluous and inaccurate unless the assumption
        /// is made that it signifies vector addition.  For example, with higher dimensional vectors, 
        /// it seems sensible to use V = (a1, a2, a3,...,an).
        /// </remarks>
        /// <example>
        /// <code>
        ///     // get the string representation of a complex number z
        ///     Complex z = new Complex(1.2345, -2.3456);
        ///     Console.WriteLine(z.ToString());
        /// 
        ///     Output:  (1.2345, -2.3456)
        /// </code>
        /// </example>
        /// </summary>
        /// <returns></returns>
        public override string ToString()
        {
            return (String.Format("( {0},  {1} )", real, imag));
        }

        /// <summary>
        /// <remarks>
        /// Allows comparison of two Complex numbers
        /// </remarks>
        /// <example>
        /// <code>
        ///     // compare two complex numbers
        ///     Complex z1 = new Complex(0.123234, -0.345456);
        ///     Complex z2 = new Complex(0.123234, -0.345455);
        ///     if (z1.Equals(z2)) Console.WriteLine("z1 equals z2"); else Console.WriteLine("z1 not equals z2");
        /// 
        ///     Output: z1 not equals z2
        /// </code>
        /// </example>
        /// </summary>
        /// <param name="obj"></param>
        /// <returns></returns>
        public override bool Equals(object obj)
        {
            return base.Equals(obj);
        }

        /// <summary>
        /// <remarks>
        /// Attempts to provide a unique integer code for this complex object. The method employed here
        /// is to use the exclusive or (XOR) of the two hash code integers obtained from the real and
        /// imaginary parts of this number.  It should be obvious that this will not always provide a 
        /// unique hash code and so should be used reservedly.
        /// </remarks>
        /// <example>
        /// <code>
        ///     // compare the hash codes of two different Complex objects
        ///     Complex z1 = new Complex(0.123234, -0.345456);
        ///     Complex z2 = new Complex(0.123234, -0.345455);
        ///     int hc1, hc2;
        ///     hc1 = z1.GetHashCode();
        ///     hc2 = z2.GetHashCode();
        ///     Console.WriteLine(string.Format("hc1 = {0}  hc2 = {1}", hc1, hc2));
        /// 
        ///     Output:  hc1 = -1809452973  hc2 = -1503159876        
        /// </code>
        /// </example>
        /// </summary>
        /// <returns></returns>
        public override int GetHashCode()
        {
            return m_real.GetHashCode() ^ m_imag.GetHashCode();
        }

        /// <summary>
        /// <remarks>
        /// Uses the string form of the rectangular representation of a complex number to 
        /// display the number in Console mode.  z = (a, b) where a = real part, b = imaginary part.
        /// Notice that the imaginary 'i' is to be understood to be a multiplier of b and has
        /// been left out. Also, the '+' sign between a and b which is sometimes used to suggest 
        /// vector addition has been omitted.
        /// </remarks>
        /// <example>
        /// <code>
        ///     // display a complex number to the screen from a Console application.
        ///     a = a.complex(2.34, -5.67);
        ///     a.Show();
        ///     
        ///     Output: ( 2.34,  -5.67 )
        /// </code>
        /// </example>
        /// </summary>
        public void Show()
        {
            Console.WriteLine(this.ToString());
        }

        /// <summary>
        /// <remarks>
        /// Shows the polar representation of the complex number z in the form in Console mode.
        /// z = |z| * e^(i * t)
        ///    where |z| = modulus of z
        ///           t  = argument of z
        ///           i  = complex(0.0, 1.0)
        /// </remarks>
        /// <example>
        /// <code>
        ///     // use ShowPolar to demonstrate the polar form of the complex number.
        ///     a = a.complex(2.34, -5.67);
        ///     Console.WriteLine(string.Format("a = {0}", a));
        ///     Console.WriteLine(a.Polar());
        ///
        ///     Output:
        ///     a = ( 2.34,  -5.67 )
        ///     6.13388131609995 * e^(i * -1.17939119866969)
        /// </code>
        /// </example>
        /// </summary>
        public void ShowPolar()
        {
            Console.WriteLine(string.Format("{0} * e^(i * {1})", this.Modulus, this.Argument));
        }

        /// <summary>
        /// <remarks>
        /// Provides a string representation of the polar form of the complex number.
        /// Compare with ShowPolar method which shows the polar representation only
        /// in Console mode.
        /// </remarks>
        /// <example>
        /// <code>
        ///     // polar verses rectangular representations
        ///     a = a.complex(2.34, -5.67);
        ///     Console.WriteLine(string.Format("a = {0}", a));
        ///     Console.WriteLine(a.Polar());
        ///
        ///     Output:
        ///     a = ( 2.34,  -5.67 )
        ///     6.13388131609995 * e^(i * -1.17939119866969)
        /// </code>
        /// </example>
        /// </summary>
        /// <returns></returns>
        public string Polar()
        {
            Complex z = new Complex(this);
            string s = z.Modulus.ToString() + " * e^(i * " + z.Argument.ToString() + ")";
            return s;
        }

        /// <summary>
        /// <remarks>
        /// Writes the rectangular form of a complex number z to the Output window in Debug mode.
        /// z = (a, b) where a = real part, b = imaginary part.
        /// </remarks>
        /// <example>
        /// <code>
        ///     // dump a complex number to the Output window in Debug mode.
        ///     a = a.complex(2.34, -5.67);
        ///     a.Dump();
        ///     
        ///     Output: a = ( 2.34,  -5.67 ) in Output window
        /// </code>
        /// </example>
        /// </summary>
        public void Dump()
        {
            Debug.WriteLine(this.ToString());
        }

        /// <summary>
        /// <remarks>
        /// Shows the polar representation of the complex number z in the form in the
        /// Output window when in Debug mode.
        /// z = |z| * e^(i * t)
        ///    where |z| = modulus of z
        ///           t  = argument of z
        ///           i  = complex(0.0, 1.0)
        /// </remarks>
        /// <example>
        /// <code>
        ///     // use DumpPolar to debug the polar form of the complex number.
        ///     a = a.complex(2.34, -5.67);
        ///     Debug.WriteLine(string.Format("a = {0}", a));
        ///     Debug.WriteLine(a.Polar());
        ///
        ///     Output:  (output window)
        ///     a = ( 2.34,  -5.67 )
        ///     6.13388131609995 * e^(i * -1.17939119866969)
        /// </code>
        /// </example>
        /// </summary>
        public void DumpPolar()
        {
            Debug.WriteLine(string.Format("{0} * e^(i * {1})", this.Modulus, this.Argument));

        }

        /// <summary>
        /// <remarks>
        /// This algorithm truncates a double type value to
        /// a specified number of decimal digits (those digits
        /// to the right of the decimal point).  It is used by
        /// Complex::Round(int ndigs).  Warning: if more than
        /// 17 significant digits are present, regardless of 
        /// their position with respect to the decimal point,
        /// overflow could occur and/or significant accuracy
        /// be lost.
        /// </remarks>
        /// </summary>
        /// <param name="dval"></param>
        /// <param name="ndigs"></param>
        /// <returns>double</returns>
        private double Truncate(double dval, int ndigs)
        {
            // MBP rounding algorithm
            // ndigs = number of decimal digits
            double rndfac = Math.Pow(10, ndigs);
            dval *= rndfac;
            dval = Math.Round(dval);
            dval /= rndfac;
            return dval;
        }

        /// <summary>
        /// <remarks>
        /// This method used the Complex::(double dval, int ndigs)
        /// method to round to a specified (ndigs) number of
        /// decimal digits (those digits to the right of the decimal 
        /// point). Warning: if more than 17 significant digits are 
        /// present, regardless of their position with respect to the 
        /// decimal point, overflow could occur and/or significant accuracy
        /// could be lost.  If that is the case, don't use this method.
        /// </remarks>
        /// </summary>
        /// <param name="ndigs"></param>
        /// <returns></returns>
        public Complex Round(int ndigs)
        {
            Complex z = new Complex(this);
            z.real = Truncate(z.real, ndigs);
            z.imag = Truncate(z.imag, ndigs);
            return z;
        }

        #endregion

    }// class Complex

}// namespace ComplexMath

