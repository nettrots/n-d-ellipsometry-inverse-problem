A C# Class to Perform Arithmetic on Complex Numbers
By Mike Pliam. 

Introduction
Although the C# language has been used quite successfully in many kinds of projects, it still remains to be proven proficient for scientific computing. Questions remain as to whether or not C# can measure up to the likes of FORTRAN and C++ for scientific and mathematical projects.

Complex Arithmetic
Arithmetic and algebra on complex numbers is important in many fields of science and engineering. Even though C# does not support an intrinsic complex number data type, it doesn't stop you from creating your own. You can do so using a struct, which has value semantics. The purpose of this article is to present my attempt to fashion a complex arithmetic class for C#.

Using the Code
The class Complex is completely contained within the namespace ComplexMath in the file Complex.cs. There are three separate program files included in the download as well as a documented help file for the project.  You can run the program either either by using the ComplexMath.dll reference, or by using the ComplexMath library code directly by simply adding the Complex.cs file to a project. You must, in either case, add the header using ComplexMath; to the main project file. If you decide to go the DLL route, first build the library (DLL) by running the CompexMath.sln file from which you can generate two ComplexMath.dll files, one each in bin/Release and bin/Debug.

Steps to use the library as a DLL:

Start a new C# console application. 
Build both the Debug and Release versions. 
Place the Debug and Release versions of the ComplexMath.dll in their respective bin folders. 
From Project/Add Reference menus, find the ComplexMath.dll for the version (debug or release) and add it to the project. 
Add the header: using ComplexMath; to your program file. 
Use the accompanying TestComplexMath.cs file as a prototype for you own program. 

Your test program should look something like the one below which briefly demonstrates many of the capabilities of the library.

[code]
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


[/code]

The accompanying TestLib program demonstrates many of the capabilities of the Complex class. One point worth making is that complex numbers can be expressed both in a rectangular form that most high-school algebra students are familiar with and in the more sophisticated polar form which requires some more advanced understanding of trigonometry, De Moivre's Theorem and Euler's identity. If you are interested, check them out on the Web. 

Suffice it to say that all of the constructor input assumes that you are using the rectangular form of the complex number, which implies that real is measured along the x-axis of the Argand plane and imaginary is measured along the y-axis of that plane. However, the program has the capacity to output the components of the polar form -- that is, the Modulus |z| and the Argument (theta) -- by using the show_polar() and Polar() methods.  Please see the accompanying ComplexMath.chm help file for the details.

Since complex numbers can be added, subtracted, multiplied and divided, it is convenient to have operators to accomplish those ends. Thus, we have incorporated those operators. 

History
There appears to be no shortage of postings offering C# code for dealing with complex numbers. These include:

Complex math library for C# and VB.NET by Karl Tarbet 
Boxing Wang's Complex struct, a comment to the Tarbet article 
Sharp3D.Math - A 3D math library for .NET by Eran Kampf 
Class to handle complex numbers natively by Dan2178, which alleges to include Over 100 methods, operators and overloads 
There are software companies selling extensive collections of C# scientific and engineering software for around $1000.00 
Although I have been involved in scientific mathematical and statistical programming for over 10 years, this article is the first that I have attempted to post. I quickly learned that there are many smarter and better programmers out there who have done a much better job at using C#. For me, this was an excellent learning experience, both in terms of learning some of the pitfalls of scientific computing with C# and having to explain to others what my code does. 

I plan to update and improve this code for my own use and for any others who wish to use it, and to do a better job in the future. You may use this code in any manner that you wish. It is always polite and professional to give the original author credit. Please report any comments, criticisms or bugs to the forum below.


Michael B. Pliam
mbpliam@pliatech.com