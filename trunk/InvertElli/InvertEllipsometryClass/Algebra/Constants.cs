namespace SbB.Diploma
{
    public abstract class Constants
    {
        public const double EPS = 1e-20;
    }

    // Перелічення положення точки відносно прямої
    public enum VertexPos
    {
        LEFT,       //зліва
        RIGHT,      //справа
        BEYOND,     //попереду
        BEHIND,     //позаду
        BETWEEN,    //між
        ORIGIN,     //початок 
        DESTINATION //кінець
    }

    public enum BoundaryType
    {
        KINEMATIC,
        STATIC,
        MORTAR,
        NONMORTAR
    }

    public enum Norma
    {
        Linear = 1,
        Euclidean,
        Maximum
    }

    public struct PairInt
    {
        public int m;
        public int n;

        public PairInt(int m, int n)
        {
            this.m = m;
            this.n = n;
        }

        public override string ToString()
        {
            return string.Format("[ {0} x {1} ]", m, n);
        }
    }

    public delegate double Phi(int i, double x, double y);

    public enum Functions
    {
        U,
        V,
        Exx,
        Eyy,
        Exy,
        Sxx,
        Syy,
        Sxy
    }

    public delegate double FunctionXY(double x, double y);

    public delegate void EmptyDelegate();
}