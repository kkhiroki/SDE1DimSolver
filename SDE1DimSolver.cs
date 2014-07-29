using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace SDESolver
{
   
    class Program
    {
        static void Main(string[] args)
        {
            double lam = 2.0;
            double mu = 1.0;
            double T = 1.0;
            int R = 20;
            int count = 256;

            Func<double, double> f = x =>
            {
                return lam*x;
            };

            Func<double,double> g = x =>
            {
              return mu*x;      
            };

      //      SDE1DimSolver sol = new SDE1DimSolver(f, g, 1.0, T, R, count);
            SDESolverMH sol = new SDESolverMH(f, g, 1.0, T, R, count);
            double[] res = sol.solveSDE();
            double[] ex = sol.getExactSol();
            Util.plotSaveVector(res,T,"SDENum","\"blue\"",false);
            Util.plotSaveVector(ex, T, "SDEExact","\"red\"",true);
            Console.Read();
        }
    }

    public static class Util
    {
        static Random Rn = new Random(System.Environment.TickCount);
        static double Alpha, Beta, BoxMuller1, BoxMuller2;
        static bool Flag = false;
        static System.Diagnostics.Process gnuplot;
        static string command;

        public static double NormalDistribution(double mu, double sigma)
        {
            if (Flag)
            {
                Alpha = Rn.NextDouble();
                Beta = Rn.NextDouble() * Math.PI * 2;
                BoxMuller1 = Math.Sqrt(-2 * Math.Log(Alpha));
                BoxMuller2 = Math.Sin(Beta);
            }
            else
                BoxMuller2 = Math.Cos(Beta);
            Flag = !Flag;
            return sigma * (BoxMuller1 * BoxMuller2) + mu;
        }

        public static void plotSaveVector(double[] vec,double T,string fName, string color, bool isContinue)
        {

            //Save
            string date = DateTime.Now.ToString("yyyyMMdd");
            string createDirName = "\\" + date;
            string fileName = fName + ".csv";

        
                // データ保存
            System.IO.StreamWriter sw = new System.IO.StreamWriter(@fileName);
            string dateTime = DateTime.Now.ToString("yyyy/MM/dd HH:mm:ss");
            sw.WriteLine("#" + dateTime);
            sw.WriteLine("#count,X");

            double dt = T / (double)vec.Length;
            for (int i = 0; i < vec.Length ; i++)
            {
                sw.Write(i*dt + "," + vec[i] + "\n");
            }
            sw.Close();

            //Plot
            //外部アプリケーションを呼び出すにはProcessクラスを使えばよい
            if ((gnuplot == null))
            {
                gnuplot = new System.Diagnostics.Process();
                gnuplot.StartInfo.FileName = "pgnuplot.exe";
                gnuplot.StartInfo.UseShellExecute = false;
                gnuplot.StartInfo.RedirectStandardInput = true;
                gnuplot.StartInfo.RedirectStandardOutput = true;
                gnuplot.Start();
            }
            
            if (isContinue)
            {
                command +=  "'" + fileName + "' using '%lf,%lf' w l title '"+ fName +"' lc rgb " + color;
                gnuplot.StandardInput.WriteLine(command);
            }
            else
            {
                command = ("plot '" + fileName + "' using '%lf,%lf' w lp title '"+ fName + "' lc rgb " + color + ",");
            }
            
        }
    } 

    class SDE1DimSolver
    {
        //Solving SDE dX = f dt + g dW

        protected Func<double, double> f;
        protected Func<double, double> g;
        protected double iv;//initial value
        protected double T;//final time
        protected int R;
        protected int count;
        protected double[] exact;

        public SDE1DimSolver(Func<double,double> _f, Func<double ,double> _g, double _iv, double _T, int _R, int _count)
        {
            this.f = _f;
            this.g = _g;
            this.iv = _iv;
            this.T = _T;
            this.R = _R;
            this.count = _count;
            this.exact = new double[this.R*this.count+1];
        }

        public virtual double[] solveSDE()
        {
            double[] resVec = new double[this.count+1];
            double[] WVec = new double[this.R * this.count + 1];
            double dt = T / count;
            double dtt = dt / R;


            resVec[0] = this.iv;
            this.exact[0] = this.iv;
            double W = 0.0;
            double lam = this.f(1.0);
            double mu = this.g(1.0);
            double dW = 0.0;

            WVec[0] = 0.0;
            for (int i = 0; i < this.R * this.count; i++)
            {

                dW = Util.NormalDistribution(0, Math.Sqrt(dtt));
                W += dW;
                WVec[i + 1] = dW;
                exact[i + 1] = Math.Exp((lam - 0.5 * mu * mu) * (i * dtt) + mu * W);
            }

            double ddW;
            for (int i = 0; i < this.count; i++)
            {
                ddW = 0.0;
                    double Xi = resVec[i];
                    for (int j = 0; j < R; j++) 
                    {
                        ddW += WVec[this.R*i + j];
                    }
                resVec[i + 1] = Xi + this.f(Xi) * dt + this.g(Xi) * ddW;  
            }
            return resVec;
        }

        public double[] getExactSol()
        {
            return this.exact;
        }

    }

    class SDESolverMH : SDE1DimSolver
    {
        public SDESolverMH(Func<double,double> _f, Func<double ,double> _g, double _iv, double _T, int _R, int _count)
            :base(_f,_g,_iv,_T,_R,_count)
        {
        }
        public override double[] solveSDE()
        {
            double[] resVec = new double[this.count + 1];
            double[] WVec = new double[this.R * this.count + 1];
            double dt = T / this.count;
            double dtt = dt / R;



            resVec[0] = this.iv;
            this.exact[0] = this.iv;
            double W = 0.0;
            double lam = this.f(1.0);
            double mu = this.g(1.0);
            double dW = 0.0;

            WVec[0] = 0.0;
            for (int i = 0; i < this.R * this.count; i++)
            {

                dW = Util.NormalDistribution(0, Math.Sqrt(dtt));
                W += dW;
                WVec[i + 1] = dW;
                exact[i + 1] = Math.Exp((lam - 0.5 * mu * mu) * (i * dtt) + mu * W);
            }

            double ddW;
            double dx = 0.001;
            for (int i = 0; i < this.count; i++)
            {
                ddW = 0.0;
                double Xi = resVec[i];
                for (int j = 0; j < R; j++)
                {
                    ddW += WVec[this.R * i + j];
                }
                double LTerm = 0.5 * this.g(Xi) * this.dg(Xi,dx) * (ddW * ddW - dt);
                resVec[i + 1] = Xi + this.f(Xi) * dt + this.g(Xi) * ddW + LTerm;
            }
            return resVec;
        }

        double df (double x,double dx)
        {
            return ( this.f(x + dx) - this.f(x) ) / dx;

        }
        double dg(double x, double dx)
        {
            return (this.g(x + dx) - this.g(x)) / dx;
        }

    }
}
