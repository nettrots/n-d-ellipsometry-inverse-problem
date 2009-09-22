
using System;
using System.Collections.Generic;
using ComplexMath;
using SbB.Diploma;


namespace InvertEllipsometryClass
{
    public class Experiment
    {
        private double psi;
        private double delta;
        private double incidentAngle;
        private Complex rho;
        public Experiment(double psi,double delta,double incidentAngle)
        {
            this.psi = psi;
            this.delta = delta;
            this.incidentAngle = incidentAngle;
            calcPho();
        }
        public double Psi
        {
            get { return psi; }
            set { psi = value; }
        }

        public double Delta
        {
            get { return delta; }
            set { delta = value; }
        }

        public double IncidentAngle
        {
            get { return incidentAngle; }
            set { incidentAngle = value; }
        }

        private void calcPho()
        {
            rho = Math.Tan(psi)*Complex.Exp(new Complex(0, delta));
        }
        public Complex Rho
        {
            get
            {
                if (rho == (null))calcPho();
                return rho;
            }
           
        }
    }
    public class Material
    {
        private int number;
        private double d;
        private Complex angle;
        private Complex n;
        private Complex beta;
        private Matrix l;
        private Complex incAngle;
        private Complex incN;
        private double lambda;

        public Material(Complex n)
        {
            this.n = n;
        }


        public Complex Beta
        {
            get
            {
                beta=getBeta();
                return beta;
            }
        }

        public Complex N
        {
            get { return n; }
            set { n = value; }
        }

        public double D
        {
            get { return d; }
            set { d = value; }
        }

        public Complex Angle
        {
            get { return angle; }
            set { angle = value; }
        }


        public Matrix calcL(Complex incAngle, Complex incN, double d, double lambda)
        {
            this.d = d;
            Func<Complex, Complex> arcsin = (x) =>
            {
                Complex ci = new Complex(0, 1);
                Complex c1 = new Complex(1, 0);
                return (new Complex(Math.PI,0))+ci * Complex.Log(ci * x + Complex.Sqrt(c1 - x * x));
            };
            //    -i ln(iz +sqrt 1-z^2)}
            this.angle = arcsin(incN / N * Complex.Sin(incAngle));
            Complex c = arcsin(Complex.Sin(incAngle));
            this.incAngle = incAngle;
            this.incN = incN;
            this.lambda = lambda;

            Complex i = new Complex(0, 1);
            return l = new Matrix(new Complex[2][]
                                      {
                                          new[]
                                              {
                                                  Complex.Exp(i*Beta), new Complex(0.0, 0.0)
                                              }
                                          , new[]
                                               {
                                                new Complex(0.0, 0.0), Complex.Exp(-i*Beta)
                                               }
                                      });
        }

    

        public Complex getBeta()
        {
            return new Complex(2*Math.PI*(d/lambda)*(n*n - incN*incN*Complex.Pow(Complex.Sin(incAngle), 2)));
        }

    }
    public class Ambient:Material
    {
        public Ambient(Complex n) : base(n)
        {
        }
    }
    public class Substract:Material
    {
        public Substract(Complex n) : base(n)
        {
        }
    }
    public class Film:Material
    {
        private bool recalc = false;

        public Film(Complex n) : base(n)
        {
        }

        public bool Recalc
        {
            get { return recalc; }
            set { recalc = value; }
        }
    }
    public class Functional
    {
        private Experiment expData;
        private FilmsModel calsData;

        public Functional(double psi,double delta,double incidentAngle, Complex[] N,double[] d,double lambda)
        {
            expData=new Experiment( psi, delta, incidentAngle);
            calsData = new FilmsModel(N, d, incidentAngle, lambda);
        }
        public double functional(double n, double d)
        {
            //if (n <= 0 || d <= 0){ Random r=new Random(1); return r.Next()* 10e30;}
            calsData.changeParams(new Complex(n,0),d );
            Complex x1 = expData.Rho;
            Complex x2 = calsData.PhoExp();
            double t = Math.Atan(x1.Modulus);
            double t1 = Math.Atan(x2.Modulus);
            double c = (Complex.Log(x1 / x1.Modulus)).imag;
            double c1 = (Complex.Log(x2 / x2.Modulus)).imag;
            double k =(t - t1)*(t - t1), k1 = (c - c1)*(c - c1);
            return 
                k+k1 ;
        }
        public double functional(double n, double d, int psi, int delta)
        {
            //if (n <= 0 || d <= 0){ Random r=new Random(1); return r.Next()* 10e30;}
            calsData.changeParams(new Complex(n, 0), d);
            Complex x1 = expData.Rho;
            Complex x2 = calsData.PhoExp();
            double t = Math.Atan(x1.Modulus);
            double t1 = Math.Atan(x2.Modulus);
            double c = (Complex.Log(x1 / x1.Modulus)).imag;
            double c1 = (Complex.Log(x2 / x2.Modulus)).imag;
            return psi*100 * (t - t1) * (t - t1) + delta*(c - c1) * (c - c1);
        }

    }
    public class FilmsModel
    {
        List<Material> stack=new List<Material>();
        private List<Matrix> interfaceMatrixPs = new List<Matrix>();
        private List<Matrix> interfaceMatrixSs = new List<Matrix>();
        private List<Matrix> lMatricies = new List<Matrix>();
        private Matrix scaterringMatrixP, scaterringMatrixS;
        private Complex incedentAngle;
        private double lambda;
        public FilmsModel(Complex[] N, double[] D, double incedentAngle, double lambda)
        {
            this.lambda = lambda;
           stack=new List<Material>();
           stack.Add(new Ambient(N[0]));
            stack[0].Angle =new Complex(incedentAngle,0);
            int i = 1,j=0;
            foreach (var d in D)
            {
                stack.Add(new Film(N[i++]));
                stack[stack.Count - 1].D = d;
                if(j++==0) (stack[stack.Count - 1] as Film).Recalc = true;

            }
            stack.Add(new Substract(N[i]));
            this.incedentAngle =new Complex(incedentAngle,0);
        }

        public void changeParams(Complex n,double d)
        {
            stack[1].N = n;
            stack[1].D = d;
        }

        public Complex PhoExp()
        {
            calcMatrix();
            return scaterringMatrixP[1][0] / scaterringMatrixP[0][0] * scaterringMatrixS[0][0] / scaterringMatrixS[1][0];

        }
        public Matrix ScaterringMatrixP
        {
            get
            {
                if (scaterringMatrixP == null)
                    calcMatrix();
                return scaterringMatrixP;
            }
        }
        public Matrix ScaterringMatrixS
        {
            get
            {
                if (scaterringMatrixS == null)
                    calcMatrix();
                return scaterringMatrixS;
            }
        }
        private void calcMatrix()
        {
            scaterringMatrixP = new Matrix(2, 2);
            scaterringMatrixP[0][0]=new Complex(1,0);
            scaterringMatrixP[1][1] = new Complex(1, 0);
            scaterringMatrixS = new Matrix(2, 2);
            scaterringMatrixS[0][0]=new Complex(1,0);
            scaterringMatrixS[1][1] = new Complex(1, 0);


           for(int i=0;i<=stack.Count-2;i++)
           {
               
               if ((stack[i] is Film) && ((stack[i] as Film).Recalc || lMatricies.Count-1<i))
               {
                   
                    if ((stack[i] as Film).Recalc && !(lMatricies.Count - 1 < i))
                   {
                       //TODO: some thing with incedentAngle...

                       lMatricies[i] = (stack[i].calcL(stack[i-1].Angle, stack[i - 1].N, stack[i].D, lambda));
                       interfaceMatrixPs[i]=(interfaceMatrixS(stack[i], stack[i + 1]));
                       interfaceMatrixSs[i]=(interfaceMatrixS(stack[i], stack[i + 1]));
                   }
                   else
                   {
                       lMatricies.Add(stack[i].calcL(stack[i - 1].Angle, stack[i].N, stack[i].D, lambda));
                       interfaceMatrixPs.Add(interfaceMatrixS(stack[i], stack[i + 1]));
                       interfaceMatrixSs.Add(interfaceMatrixS(stack[i], stack[i + 1]));
                       
                   }
               }
            
               if ((stack[i] is Ambient) )
                   if (interfaceMatrixPs.Count - 1 < i)
                   {
                       interfaceMatrixPs.Add(interfaceMatrixP(stack[i], stack[i + 1]));
                       Matrix te = new Matrix(2, 2);
                       te[0][0] = new Complex(1, 0);
                       te[1][1] = new Complex(1, 0);
                       lMatricies.Add(te);
                   }
                   else
                       interfaceMatrixPs[i] = (interfaceMatrixP(stack[i], stack[i + 1]));

               if ((stack[i] is Ambient))
                   if (interfaceMatrixSs.Count - 1 < i)
                       interfaceMatrixSs.Add(interfaceMatrixS(stack[i], stack[i + 1]));
                   else
                       interfaceMatrixSs[i] = (interfaceMatrixS(stack[i], stack[i + 1]));
               
               scaterringMatrixP *= lMatricies[i] * interfaceMatrixPs[i];
               scaterringMatrixS *= lMatricies[i] * interfaceMatrixSs[i];

           }
            
        }

        private Matrix interfaceMatrixP(Material material, Material material_2)
        {
            Matrix temp = new Matrix(2, 2);
            temp[0][0] = new Complex(1, 0);
            temp[1][1] = new Complex(1, 0);
            Complex t1 = material.N * Complex.Cos(material.Angle)
                    ,
                    t2 = material_2.N*
                         Complex.Sqrt(new Complex(1, 0) -
                                      Complex.Pow(material.N / material_2.N * Complex.Sin(material.Angle), 2)),
                    t0 = material_2.N*Complex.Cos(incedentAngle)
                ;
            temp[0][1] = temp[1][0] = (t1 - t2)/(t1 + t2);
            temp = (t1 + t2) / (2 * t0) * temp;
            return temp;
        }

        private Matrix interfaceMatrixS(Material material, Material material_2)
        {
            Matrix temp = new Matrix(2, 2);
            temp[0][0] = new Complex(1, 0);
            temp[1][1] = new Complex(1, 0);
            Complex t11 = material_2.N * Complex.Cos(material.Angle)
                    ,
                    t22 = material.N*
                          Complex.Sqrt(new Complex(1, 0) -
                                       Complex.Pow(material.N / material_2.N * Complex.Sin(material.Angle), 2))

                ;
            temp[0][1] = temp[1][0] = (t11 - t22)/(t11 + t22);
            temp = (t11 + t22) / (2 * t11) *temp;
            return temp;
        }

    }

}
