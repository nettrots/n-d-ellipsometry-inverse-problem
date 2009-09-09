
using System;

namespace alglib
{
    public class principalaxis
    {
        /*
        This members must be defined by you:
        static double f(double[] x)
        */


        public static double fx = 0;
        public static double ldt = 0;
        public static double dmin = 0;
        public static int nl = 0;
        public static int nf = 0;
        public static double[,] v = new double[0,0];
        public static double[] q0 = new double[0];
        public static double[] q1 = new double[0];
        public static double qa = 0;
        public static double qb = 0;
        public static double qc = 0;
        public static double qd0 = 0;
        public static double qd1 = 0;
        public static double qf1 = 0;


        /*************************************************************************
        Поиск минимума функции N аргументов (N>=2) методом главных направлений.

        Исходный код получен путем перевода оригинального кода из языка Fortran.
        Автор оригинального кода и алгоритма в целом: Р.П.Брент.

        Входные параметры:
            N   -   Размерность задачи. Не меньше 2.
            X   -   массив с нумерацией элементов от 1 до N. Начальная точка.
            T0  -   Требуемая точность. Алгоритм пытается найти такую точку, которая
                    отличается от локального минимума не более, чем на
                    T0+SquareRoot(MachineEpsilon)*NORM(Minimum).
            H0  -   Максимально допустимая длина шага. Должно по порядку
                    соответствовать максимально возможному расстоянию от начальной
                    точки до минимума. Если задано слишком большое или слишком маленькое
                    значение, начальная скорость сходимости может быть слишком низкой.

        Выходные параметры:
            X   -   найденное приближение к минимуму

        Результат:
            Значение функции в найденной точке
        *************************************************************************/
        public static double principalaxisminimize(int n,
            ref double[] x,
            double t0,
            double h0)
        {
            double result = 0;
            bool illc = new bool();
            int kl = 0;
            int kt = 0;
            int ktm = 0;
            int i = 0;
            int j = 0;
            int k = 0;
            int k2 = 0;
            int km1 = 0;
            int klmk = 0;
            int ii = 0;
            int im1 = 0;
            double s = 0;
            double sl = 0;
            double dn = 0;
            double f1 = 0;
            double lds = 0;
            double t = 0;
            double h = 0;
            double sf = 0;
            double df = 0;
            double m2 = 0;
            double m4 = 0;
            double small = 0;
            double vsmall = 0;
            double large = 0;
            double vlarge = 0;
            double scbd = 0;
            double ldfac = 0;
            double t2 = 0;
            double dni = 0;
            double value = 0;
            double[] d = new double[0];
            double[] y = new double[0];
            double[] z = new double[0];
            double temparrayelement = 0;

            d = new double[n+1];
            y = new double[n+1];
            z = new double[n+1];
            q0 = new double[n+1];
            q1 = new double[n+1];
            v = new double[n+1, n+1];
            small = AP.Math.Sqr(AP.Math.MachineEpsilon);
            vsmall = small*small;
            large = 1/small;
            vlarge = 1/vsmall;
            m2 = Math.Sqrt(AP.Math.MachineEpsilon);
            m4 = Math.Sqrt(m2);
            scbd = 1;
            illc = false;
            ktm = 1;
            ldfac = 0.01;
            if( illc )
            {
                ldfac = 0.1;
            }
            kt = 0;
            nl = 0;
            nf = 1;
            fx = f(x);
            qf1 = fx;
            t = small+Math.Abs(t0);
            t2 = t;
            dmin = small;
            h = h0;
            if( h<100*t )
            {
                h = 100*t;
            }
            ldt = h;
            for(i=1; i<=n; i++)
            {
                for(j=1; j<=n; j++)
                {
                    v[i,j] = 0;
                }
                v[i,i] = 1;
            }
            d[1] = 0;
            qd0 = 0;
            for(i=1; i<=n; i++)
            {
                q0[i] = x[i];
                q1[i] = x[i];
            }
            while( true )
            {
                sf = d[1];
                d[1] = 0;
                s = 0;
                value = fx;
                temparrayelement = d[1];
                directionalminimize(n, 1, 2, ref temparrayelement, ref s, ref value, false, ref x, ref t, AP.Math.MachineEpsilon, h);
                d[1] = temparrayelement;
                if( s<=0 )
                {
                    for(i=1; i<=n; i++)
                    {
                        v[i,1] = -v[i,1];
                    }
                }
                if( sf<=0.9*d[1] | 0.9*sf>=d[1] )
                {
                    for(i=2; i<=n; i++)
                    {
                        d[i] = 0;
                    }
                }
                for(k=2; k<=n; k++)
                {
                    for(i=1; i<=n; i++)
                    {
                        y[i] = x[i];
                    }
                    sf = fx;
                    if( kt>0 )
                    {
                        illc = true;
                    }
                    while( true )
                    {
                        kl = k;
                        df = 0;
                        if( illc )
                        {
                            for(i=1; i<=n; i++)
                            {
                                s = (0.1*ldt+t2*Math.Pow(10, kt))*(AP.Math.RandomReal()-0.5);
                                z[i] = s;
                                for(j=1; j<=n; j++)
                                {
                                    x[j] = x[j]+s*v[j,i];
                                }
                            }
                            fx = f(x);
                            nf = nf+1;
                        }
                        for(k2=k; k2<=n; k2++)
                        {
                            sl = fx;
                            s = 0;
                            value = fx;
                            temparrayelement = d[k2];
                            directionalminimize(n, k2, 2, ref temparrayelement, ref s, ref value, false, ref x, ref t, AP.Math.MachineEpsilon, h);
                            d[k2] = temparrayelement;
                            if( illc )
                            {
                                s = d[k2]*AP.Math.Sqr(s+z[k2]);
                            }
                            else
                            {
                                s = sl-fx;
                            }
                            if( df>s )
                            {
                                continue;
                            }
                            df = s;
                            kl = k2;
                        }
                        if( illc | df>=Math.Abs(100*AP.Math.MachineEpsilon*fx) )
                        {
                            break;
                        }
                        illc = true;
                    }
                    km1 = k-1;
                    for(k2=1; k2<=km1; k2++)
                    {
                        s = 0;
                        value = fx;
                        temparrayelement = d[k2];
                        directionalminimize(n, k2, 2, ref temparrayelement, ref s, ref value, false, ref x, ref t, AP.Math.MachineEpsilon, h);
                        d[k2] = temparrayelement;
                    }
                    f1 = fx;
                    fx = sf;
                    lds = 0;
                    for(i=1; i<=n; i++)
                    {
                        sl = x[i];
                        x[i] = y[i];
                        sl = sl-y[i];
                        y[i] = sl;
                        lds = lds+sl*sl;
                    }
                    lds = Math.Sqrt(lds);
                    if( lds>small )
                    {
                        klmk = kl-k;
                        if( klmk>=1 )
                        {
                            for(ii=1; ii<=klmk; ii++)
                            {
                                i = kl-ii;
                                for(j=1; j<=n; j++)
                                {
                                    v[j,i+1] = v[j,i];
                                }
                                d[i+1] = d[i];
                            }
                        }
                        d[k] = 0;
                        for(i=1; i<=n; i++)
                        {
                            v[i,k] = y[i]/lds;
                        }
                        value = f1;
                        temparrayelement = d[k];
                        directionalminimize(n, k, 4, ref temparrayelement, ref lds, ref value, true, ref x, ref t, AP.Math.MachineEpsilon, h);
                        d[k] = temparrayelement;
                        if( lds<=0 )
                        {
                            lds = -lds;
                            for(i=1; i<=n; i++)
                            {
                                v[i,k] = -v[i,k];
                            }
                        }
                    }
                    ldt = ldfac*ldt;
                    if( ldt<lds )
                    {
                        ldt = lds;
                    }
                    t2 = 0;
                    for(i=1; i<=n; i++)
                    {
                        t2 = t2+AP.Math.Sqr(x[i]);
                    }
                    t2 = m2*Math.Sqrt(t2)+t;
                    if( ldt>0.5*t2 )
                    {
                        kt = -1;
                    }
                    kt = kt+1;
                    if( kt>ktm )
                    {
                        result = fx;
                        return result;
                    }
                }
                quadsearch(n, ref x, ref t, AP.Math.MachineEpsilon, h);
                dn = 0;
                for(i=1; i<=n; i++)
                {
                    d[i] = 1/Math.Sqrt(d[i]);
                    if( dn<d[i] )
                    {
                        dn = d[i];
                    }
                }
                for(j=1; j<=n; j++)
                {
                    s = d[j]/dn;
                    for(i=1; i<=n; i++)
                    {
                        v[i,j] = s*v[i,j];
                    }
                }
                if( scbd>1 )
                {
                    s = vlarge;
                    for(i=1; i<=n; i++)
                    {
                        sl = 0;
                        for(j=1; j<=n; j++)
                        {
                            sl = sl+v[i,j]*v[i,j];
                        }
                        z[i] = Math.Sqrt(sl);
                        if( z[i]<m4 )
                        {
                            z[i] = m4;
                        }
                        if( s>z[i] )
                        {
                            s = z[i];
                        }
                    }
                    for(i=1; i<=n; i++)
                    {
                        sl = s/z[i];
                        z[i] = 1/sl;
                        if( z[i]>scbd )
                        {
                            sl = 1/scbd;
                            z[i] = scbd;
                        }
                        for(j=1; j<=n; j++)
                        {
                            v[i,j] = sl*v[i,j];
                        }
                    }
                }
                for(i=2; i<=n; i++)
                {
                    im1 = i-1;
                    for(j=1; j<=im1; j++)
                    {
                        s = v[i,j];
                        v[i,j] = v[j,i];
                        v[j,i] = s;
                    }
                }
                svdfit(n, vsmall, ref v, ref d);
                if( scbd>1 )
                {
                    for(i=1; i<=n; i++)
                    {
                        s = z[i];
                        for(j=1; j<=n; j++)
                        {
                            v[i,j] = s*v[i,j];
                        }
                    }
                    for(i=1; i<=n; i++)
                    {
                        s = 0;
                        for(j=1; j<=n; j++)
                        {
                            s = s+AP.Math.Sqr(v[j,i]);
                        }
                        s = Math.Sqrt(s);
                        d[i] = s*d[i];
                        s = 1/s;
                        for(j=1; j<=n; j++)
                        {
                            v[j,i] = s*v[j,i];
                        }
                    }
                }
                for(i=1; i<=n; i++)
                {
                    dni = dn*d[i];
                    if( dni>large )
                    {
                        d[i] = vsmall;
                    }
                    else
                    {
                        if( dni<small )
                        {
                            d[i] = vlarge;
                        }
                        else
                        {
                            d[i] = 1/(dni*dni);
                        }
                    }
                }
                svdsort(n, ref d, ref v);
                dmin = d[n];
                if( dmin<small )
                {
                    dmin = small;
                }
                illc = false;
                if( m2*d[1]>dmin )
                {
                    illc = true;
                }
            }
            return result;
        }


        /*************************************************************************
        Служебная подпрограмма SVD-разложения.

        Улучшенная версия, адаптированная под данную задачу.
        *************************************************************************/
        private static void svdfit(int n,
            double tol,
            ref double[,] ab,
            ref double[] q)
        {
            int i = 0;
            int j = 0;
            int k = 0;
            int ii = 0;
            int kk = 0;
            int ll2 = 0;
            int l = 0;
            int kt = 0;
            int lp1 = 0;
            int l2 = 0;
            double af = 0;
            double c = 0;
            double g = 0;
            double s = 0;
            double h = 0;
            double x = 0;
            double y = 0;
            double z = 0;
            double temp = 0;
            double[] e = new double[0];
            double eps = 0;
            int flag1 = 0;

            e = new double[n+1];
            if( n==1 )
            {
                q[1] = ab[1,1];
                ab[1,1] = 1;
                return;
            }
            eps = AP.Math.MachineEpsilon;
            g = 0;
            x = 0;
            for(i=1; i<=n; i++)
            {
                e[i] = g;
                s = 0;
                l = i+1;
                for(j=i; j<=n; j++)
                {
                    s = s+AP.Math.Sqr(ab[j,i]);
                }
                g = 0;
                if( s>=tol )
                {
                    af = ab[i,i];
                    g = Math.Sqrt(s);
                    if( af>=0 )
                    {
                        g = -g;
                    }
                    h = af*g-s;
                    ab[i,i] = af-g;
                    if( l<=n )
                    {
                        for(j=l; j<=n; j++)
                        {
                            af = 0;
                            for(k=i; k<=n; k++)
                            {
                                af = af+ab[k,i]*ab[k,j];
                            }
                            af = af/h;
                            for(k=i; k<=n; k++)
                            {
                                ab[k,j] = ab[k,j]+af*ab[k,i];
                            }
                        }
                    }
                }
                q[i] = g;
                s = 0;
                if( i!=n )
                {
                    for(j=l; j<=n; j++)
                    {
                        s = s+ab[i,j]*ab[i,j];
                    }
                }
                g = 0;
                if( s>=tol )
                {
                    if( i!=n )
                    {
                        af = ab[i,i+1];
                    }
                    g = Math.Sqrt(s);
                    if( af>=0 )
                    {
                        g = -g;
                    }
                    h = af*g-s;
                    if( i!=n )
                    {
                        ab[i,i+1] = af-g;
                        for(j=l; j<=n; j++)
                        {
                            e[j] = ab[i,j]/h;
                        }
                        for(j=l; j<=n; j++)
                        {
                            s = 0;
                            for(k=l; k<=n; k++)
                            {
                                s = s+ab[j,k]*ab[i,k];
                            }
                            for(k=l; k<=n; k++)
                            {
                                ab[j,k] = ab[j,k]+s*e[k];
                            }
                        }
                    }
                }
                y = Math.Abs(q[i])+Math.Abs(e[i]);
                if( y>x )
                {
                    x = y;
                }
            }
            ab[n,n] = 1;
            g = e[n];
            l = n;
            for(ii=2; ii<=n; ii++)
            {
                i = n-ii+1;
                if( g!=0 )
                {
                    h = ab[i,i+1]*g;
                    for(j=l; j<=n; j++)
                    {
                        ab[j,i] = ab[i,j]/h;
                    }
                    for(j=l; j<=n; j++)
                    {
                        s = 0;
                        for(k=l; k<=n; k++)
                        {
                            s = s+ab[i,k]*ab[k,j];
                        }
                        for(k=l; k<=n; k++)
                        {
                            ab[k,j] = ab[k,j]+s*ab[k,i];
                        }
                    }
                }
                for(j=l; j<=n; j++)
                {
                    ab[i,j] = 0;
                    ab[j,i] = 0;
                }
                ab[i,i] = 1;
                g = e[i];
                l = i;
            }
            eps = eps*x;
            for(kk=1; kk<=n; kk++)
            {
                k = n-kk+1;
                kt = 0;
                while( true )
                {
                    kt = kt+1;
                    if( kt>30 )
                    {
                        e[k] = 0;
                    }
                    flag1 = 110;
                    for(ll2=1; ll2<=k; ll2++)
                    {
                        l2 = k-ll2+1;
                        l = l2;
                        if( Math.Abs(e[l])<eps )
                        {
                            flag1 = 120;
                            break;
                        }
                        if( l==1 )
                        {
                            continue;
                        }
                        if( Math.Abs(q[l-1])<=eps )
                        {
                            break;
                        }
                    }
                    if( flag1!=120 )
                    {
                        c = 0;
                        s = 1;
                        for(i=l; i<=k; i++)
                        {
                            af = s*e[i];
                            e[i] = c*e[i];
                            if( Math.Abs(af)<=eps )
                            {
                                break;
                            }
                            g = q[i];
                            if( Math.Abs(af)<Math.Abs(g) )
                            {
                                h = Math.Abs(g)*Math.Sqrt(1+AP.Math.Sqr(af/g));
                            }
                            else
                            {
                                if( af==0 )
                                {
                                    h = 0;
                                }
                                else
                                {
                                    h = Math.Abs(af)*Math.Sqrt(1+AP.Math.Sqr(g/af));
                                }
                            }
                            q[i] = h;
                            if( h==0 )
                            {
                                g = 1;
                                h = 1;
                            }
                            c = g/h;
                            s = -(af/h);
                        }
                    }
                    z = q[k];
                    if( l==k )
                    {
                        break;
                    }
                    x = q[l];
                    y = q[k-1];
                    g = e[k-1];
                    h = e[k];
                    af = ((y-z)*(y+z)+(g-h)*(g+h))/(2*h*y);
                    g = Math.Sqrt(af*af+1);
                    temp = af-g;
                    if( af>=0 )
                    {
                        temp = af+g;
                    }
                    af = ((x-z)*(x+z)+h*(y/temp-h))/x;
                    c = 1;
                    s = 1;
                    lp1 = l+1;
                    if( lp1<=k )
                    {
                        for(i=lp1; i<=k; i++)
                        {
                            g = e[i];
                            y = q[i];
                            h = s*g;
                            g = g*c;
                            if( Math.Abs(af)<Math.Abs(h) )
                            {
                                z = Math.Abs(h)*Math.Sqrt(1+AP.Math.Sqr(af/h));
                            }
                            else
                            {
                                if( af==0 )
                                {
                                    z = 0;
                                }
                                else
                                {
                                    z = Math.Abs(af)*Math.Sqrt(1+AP.Math.Sqr(h/af));
                                }
                            }
                            e[i-1] = z;
                            if( z==0 )
                            {
                                af = 1;
                                z = 1;
                            }
                            c = af/z;
                            s = h/z;
                            af = x*c+g*s;
                            g = -(x*s)+g*c;
                            h = y*s;
                            y = y*c;
                            for(j=1; j<=n; j++)
                            {
                                x = ab[j,i-1];
                                z = ab[j,i];
                                ab[j,i-1] = x*c+z*s;
                                ab[j,i] = -(x*s)+z*c;
                            }
                            if( Math.Abs(af)<Math.Abs(h) )
                            {
                                z = Math.Abs(h)*Math.Sqrt(1+AP.Math.Sqr(af/h));
                            }
                            else
                            {
                                if( af==0 )
                                {
                                    z = 0;
                                }
                                else
                                {
                                    z = Math.Abs(af)*Math.Sqrt(1+AP.Math.Sqr(h/af));
                                }
                            }
                            q[i-1] = z;
                            if( z==0 )
                            {
                                af = 1;
                                z = 1;
                            }
                            c = af/z;
                            s = h/z;
                            af = c*g+s*y;
                            x = -(s*g)+c*y;
                        }
                    }
                    e[l] = 0;
                    e[k] = af;
                    q[k] = x;
                }
                if( z>=0 )
                {
                    continue;
                }
                q[k] = -z;
                for(j=1; j<=n; j++)
                {
                    ab[j,k] = -ab[j,k];
                }
            }
        }


        /*************************************************************************
        Подпрограмма минимизации функции вдоль заданного направления или вдоль
        кривой второго порядка. Далее - выдержка из оригинального комментария.

        ...THE SUBROUTINE MIN MINIMIZES F FROM X IN THE DIRECTION V(*,J) UNLESS
           J IS LESS THAN 1, WHEN A QUADRATIC SEARCH IS MADE IN THE PLANE
           DEFINED BY Q0,Q1,X.
           D2 IS EITHER ZERO OR AN APPROXIMATION TO HALF F".
           ON ENTRY, X1 IS AN ESTIMATE OF THE DISTANCE FROM X TO THE MINIMUM
           ALONG V(*,J) (OR, IF J=0, A CURVE).  ON RETURN, X1 IS THE DISTANCE
           FOUND.
           IF FK=.TRUE., THEN F1 IS FLIN(X1).  OTHERWISE X1 AND F1 ARE IGNORED
           ON ENTRY UNLESS FINAL FX IS GREATER THAN F1.
           NITS CONTROLS THE NUMBER OF TIMES AN ATTEMPT WILL BE MADE TO HALVE
           THE INTERVAL.
        *************************************************************************/
        private static void directionalminimize(int n,
            int j,
            int nits,
            ref double d2,
            ref double x1,
            ref double f1,
            bool fk,
            ref double[] x,
            ref double t,
            double machep,
            double h)
        {
            int i = 0;
            int k = 0;
            double sf1 = 0;
            double sx1 = 0;
            double xm = 0;
            double fm = 0;
            double f0 = 0;
            double s = 0;
            double temp = 0;
            double t2 = 0;
            double f2 = 0;
            double x2 = 0;
            double d1 = 0;
            bool dz = new bool();
            double m2 = 0;
            double m4 = 0;
            double small = 0;
            int cycleflag = 0;

            small = AP.Math.Sqr(machep);
            m2 = Math.Sqrt(machep);
            m4 = Math.Sqrt(m2);
            sf1 = f1;
            sx1 = x1;
            k = 0;
            xm = 0;
            fm = fx;
            f0 = fx;
            dz = d2<machep;
            s = 0;
            for(i=1; i<=n; i++)
            {
                s = s+AP.Math.Sqr(x[i]);
            }
            s = Math.Sqrt(s);
            temp = d2;
            if( dz )
            {
                temp = dmin;
            }
            t2 = m4*Math.Sqrt(Math.Abs(fx)/temp+s*ldt)+m2*ldt;
            s = m4*s+t;
            if( dz & t2>s )
            {
                t2 = s;
            }
            if( t2<small )
            {
                t = small;
            }
            if( t2>0.01*h )
            {
                t2 = 0.01*h;
            }
            if( fk & f1<=fm )
            {
                xm = x1;
                fm = f1;
            }
            if( !fk | Math.Abs(x1)<t2 )
            {
                temp = 1;
                if( x1<0 )
                {
                    temp = -1;
                }
                x1 = temp*t2;
                f1 = flin(n, j, x1, ref x);
            }
            if( f1<=fm )
            {
                xm = x1;
                fm = f1;
            }
            cycleflag = 4;
            while( true )
            {
                if( cycleflag==4 )
                {
                    if( dz )
                    {
                        x2 = -x1;
                        if( f0>=f1 )
                        {
                            x2 = 2*x1;
                        }
                        f2 = flin(n, j, x2, ref x);
                        if( f2<=fm )
                        {
                            xm = x2;
                            fm = f2;
                        }
                        d2 = (x2*(f1-f0)-x1*(f2-f0))/(x1*x2*(x1-x2));
                    }
                    d1 = (f1-f0)/x1-x1*d2;
                    dz = true;
                    if( d2>small )
                    {
                        x2 = -(0.5*d1)/d2;
                    }
                    else
                    {
                        x2 = h;
                        if( d1>=0 )
                        {
                            x2 = -x2;
                        }
                    }
                    if( Math.Abs(x2)>h )
                    {
                        if( x2<=0 )
                        {
                            x2 = -h;
                        }
                        else
                        {
                            x2 = h;
                        }
                    }
                }
                f2 = flin(n, j, x2, ref x);
                if( k>=nits | f2<=f0 )
                {
                    break;
                }
                k = k+1;
                if( f0<f1 & x1*x2>0 )
                {
                    cycleflag = 4;
                    continue;
                }
                x2 = 0.5*x2;
                cycleflag = 11;
            }
            nl = nl+1;
            if( f2<=fm )
            {
                fm = f2;
            }
            else
            {
                x2 = xm;
            }
            if( Math.Abs(x2*(x2-x1))<=small )
            {
                if( k>0 )
                {
                    d2 = 0;
                }
            }
            else
            {
                d2 = (x2*(f1-f0)-x1*(fm-f0))/(x1*x2*(x1-x2));
            }
            if( d2<=small )
            {
                d2 = small;
            }
            x1 = x2;
            fx = fm;
            if( sf1<fx )
            {
                fx = sf1;
                x1 = sx1;
            }
            if( j==0 )
            {
                return;
            }
            for(i=1; i<=n; i++)
            {
                x[i] = x[i]+x1*v[i,j];
            }
        }


        /*************************************************************************
        Функция одной переменной, минимизируемая подпрограммой DirectionalMinimize.


        ...FLIN IS THE FUNCTION OF ONE REAL VARIABLE L THAT IS MINIMIZED
           BY THE SUBROUTINE DirectionalMinimize...
        *************************************************************************/
        private static double flin(int n,
            int j,
            double l,
            ref double[] x)
        {
            double result = 0;
            double[] t = new double[0];
            int i = 0;

            t = new double[n+1];
            if( j!=0 )
            {
                for(i=1; i<=n; i++)
                {
                    t[i] = x[i]+l*v[i,j];
                }
            }
            else
            {
                qa = l*(l-qd1)/(qd0*(qd0+qd1));
                qb = (l+qd0)*(qd1-l)/(qd0*qd1);
                qc = l*(l+qd0)/(qd1*(qd0+qd1));
                for(i=1; i<=n; i++)
                {
                    t[i] = qa*q0[i]+qb*x[i]+qc*q1[i];
                }
            }
            nf = nf+1;
            result = f(t);
            return result;
        }


        /*************************************************************************
        Подпрограмма сортировки столбцов матрицы главных направлений в соответствии
        с их сингулярными значениями
        *************************************************************************/
        private static void svdsort(int n,
            ref double[] vecd,
            ref double[,] matv)
        {
            double s = 0;
            int i = 0;
            int j = 0;
            int k = 0;
            int nm1 = 0;

            if( n==1 )
            {
                return;
            }
            nm1 = n-1;
            for(i=1; i<=nm1; i++)
            {
                k = i;
                s = vecd[i];
                for(j=i+1; j<=n; j++)
                {
                    if( vecd[j]<=s )
                    {
                        continue;
                    }
                    k = j;
                    s = vecd[j];
                }
                if( k<=i )
                {
                    continue;
                }
                vecd[k] = vecd[i];
                vecd[i] = s;
                for(j=1; j<=n; j++)
                {
                    s = matv[j,i];
                    matv[j,i] = matv[j,k];
                    matv[j,k] = s;
                }
            }
        }


        /*************************************************************************
        ...QUAD LOOKS FOR THE MINIMUM OF F ALONG A CURVE DEFINED BY Q0,Q1,X...
        *************************************************************************/
        private static void quadsearch(int n,
            ref double[] x,
            ref double t,
            double machep,
            double h)
        {
            int i = 0;
            double s = 0;
            double value = 0;
            double l = 0;

            s = fx;
            fx = qf1;
            qf1 = s;
            qd1 = 0;
            for(i=1; i<=n; i++)
            {
                s = x[i];
                l = q1[i];
                x[i] = l;
                q1[i] = s;
                qd1 = qd1+AP.Math.Sqr(s-l);
            }
            qd1 = Math.Sqrt(qd1);
            l = qd1;
            s = 0;
            if( qd0<=0 | qd1<=0 | nl<3*n*n )
            {
                fx = qf1;
                qa = 0;
                qb = qa;
                qc = 1;
            }
            else
            {
                value = qf1;
                directionalminimize(n, 0, 2, ref s, ref l, ref value, true, ref x, ref t, machep, h);
                qa = l*(l-qd1)/(qd0*(qd0+qd1));
                qb = (l+qd0)*(qd1-l)/(qd0*qd1);
                qc = l*(l+qd0)/(qd1*(qd0+qd1));
            }
            qd0 = qd1;
            for(i=1; i<=n; i++)
            {
                s = q0[i];
                q0[i] = x[i];
                x[i] = qa*s+qb*x[i]+qc*q1[i];
            }
        }
    }
}
