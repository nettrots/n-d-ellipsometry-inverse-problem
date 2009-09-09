
using System;
using System.Collections.Generic;
using ComplexMath;
using SbB.Diploma;


namespace InvertEllipsometryClass
{
  
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
            return scaterringMatrixP[1][0] / scaterringMatrixP[0][0] * scaterringMatrixS[0                                   ][0] / scaterringMatrixS[1][0];

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
                     

                       lMatricies[i] = (stack[i].calcL(stack[i-1].Angle, stack[i - 1].N, stack[i].D, lambda));
                       interfaceMatrixPs[i]=(interfaceMatrixP(stack[i], stack[i + 1]));
                       interfaceMatrixSs[i]=(interfaceMatrixS(stack[i], stack[i + 1]));
                   }
                   else
                   {
                       lMatricies.Add(stack[i].calcL(stack[i - 1].Angle, stack[i-1].N, stack[i].D, lambda));
                       interfaceMatrixPs.Add(interfaceMatrixP(stack[i], stack[i + 1]));
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
            Complex t1 = material_2.N * Complex.Cos(material.Angle)
                    ,
                   
                    t2 = material.N*
                         Complex.Sqrt(new Complex(1, 0) -
                                      Complex.Pow(material.N / material_2.N * Complex.Sin(material.Angle), 2)),
                    t0 = material.N * Complex.Cos(material.Angle)
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
            Complex t11 = material.N * Complex.Cos(material.Angle)
                    ,
                    t22 = material_2.N*
                          Complex.Sqrt(new Complex(1, 0) -
                                       Complex.Pow(material.N / material_2.N * Complex.Sin(material.Angle), 2))

                ;
            temp[0][1] = temp[1][0] = (t11 - t22)/(t11 + t22);
            temp = (t11 + t22) / (2 * t11) *temp;
            return temp;
        }

    }

}
