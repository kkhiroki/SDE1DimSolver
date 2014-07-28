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
            int count = 500;

            Func<double, double> f = x =>
            {
                return lam*x;
            };

            Func<double,double> g = x =>
            {
              return mu*x;      
            };

            SDE1DimSolver sol = new SDE1DimSolver(f, g, 1.0, T, count);
            double[] res = sol.solveSDE();
            Util.plotSaveVector(res,T,"SDESolve");
        }
    }

    public static class Util
    {
        static Random Rn = new Random(System.Environment.TickCount);
        static double Alpha, Beta, BoxMuller1, BoxMuller2;
        static bool Flag = false;

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

        public static void plotSaveVector(double[] vec,double T,string fName)
        {

            //Save
            string date = DateTime.Now.ToString("yyyyMMdd");
            string createDirName = fName + "\\" + date;
            string fileName = createDirName + "\\SDE.csv";

            if (System.IO.Directory.Exists(createDirName))
            {
                Console.WriteLine("ファイルが存在します");
            }
            else
            {
                // フォルダ作成
                    System.IO.Directory.CreateDirectory(@createDirName);                
            }
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
            System.Diagnostics.Process gnuplot = new System.Diagnostics.Process();

            gnuplot.StartInfo.FileName = "pgnuplot.exe";
            gnuplot.StartInfo.UseShellExecute = false;
            gnuplot.StartInfo.RedirectStandardInput = true;
            gnuplot.StartInfo.RedirectStandardOutput = true;

            gnuplot.Start();
            gnuplot.StandardInput.WriteLine("cd '" + createDirName + "'");
            gnuplot.StandardInput.WriteLine("plot 'SDE.csv' using '%lf,%lf' w l title 'SDE1Dim'");
        }
    } 

    class SDE1DimSolver
    {
        //Solving SDE dX = f dt + g dW

        Func<double, double> f;
        Func<double, double> g;
        double iv;//initial value
        double T;//final time
        int count;

        public SDE1DimSolver(Func<double,double> _f, Func<double ,double> _g, double _iv, double _T, int _count)
        {
            this.f = _f;
            this.g = _g;
            this.iv = _iv;
            this.T = _T;
            this.count = _count;
        }

        public double[] solveSDE()
        {
            double[] resVec = new double[this.count+1];
            double dt = T / count;


            resVec[0] = this.iv;
            for (int i = 0; i < this.count; i++)
            { 
                double Xi = resVec[i];
                double dW = Util.NormalDistribution(0,Math.Sqrt(dt));
                resVec[i + 1] = Xi + this.f(Xi) * dt + this.g(Xi) * dW; 
            }
            return resVec;
        }
    }
}
