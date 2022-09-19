using Microsoft.VisualBasic;
using System;
using System.Globalization;
using System.IO;
using System.Linq;
using System.Linq.Expressions;
using System.Diagnostics;

namespace Decoder
{

    public class Load
    {
        //количество точек
        //по оси x
        private int N { get; }
        //по оси y
        private int M { get; }

        private int type;
        public double LoadValue;
        public double Frequency;


        public bool FlagTimeProblem { get; }
        public double[,] F { set; get; }
        public void FullLoad()
        {
            ///Квадрат нагрузки в центре равен Q остальное 0
            if (type == 2)
                for (int x = 0; x < N; x++)
                    for (int y = 0; y < M; y++)
                        if ((1.0 / 2 - 1.0 / 8) * (N - 1) <= x && (1.0 / 2 + 1.0 / 8) * (M - 1) >= x && (1.0 / 2 - 1.0 / 8) * (M - 1) <= y && (1.0 / 2 + 1.0 / 8) * (M - 1) >= y)
                            F[x, y] = LoadValue;
                        else
                            F[x, y] = 0;

            ///Квадрат нагрузки в четвертинке равен Q остальное 0
            if (type == 3)
                for (int x = 0; x < N; x++)
                    for (int y = 0; y < M; y++)
                        if ((1.0 / 4 - 1.0 / 8) * (N - 1) <= x && (1.0 / 4 + 1.0 / 8) * (N - 1) >= x && (1.0 / 4 - 1.0 / 8) * (M - 1) <= y && (1.0 / 4 + 1.0 / 8) * (M - 1) >= y)
                            F[x, y] = LoadValue;
                        else
                            F[x, y] = 0;
            ///Равномерно распределённая нагрузка
            if (type == 1)
                for (int x = 0; x < N; x++)
                    for (int y = 0; y < M; y++)
                        F[x, y] = LoadValue;
        }
        public void FullLoadSint(double t)
        {
            if (!FlagTimeProblem)
            {
                Console.WriteLine("Error, in constructor not detected freqeuncy");
                FullLoad();
                return;
            }
            ///Квадрат нагрузки в центре равен Q остальное 0
            if (type == 2)
                for (int x = 0; x < N; x++)
                    for (int y = 0; y < M; y++)
                        if ((1.0 / 2 - 1.0 / 8) * (N - 1) <= x && (1.0 / 2 + 1.0 / 8) * (M - 1) >= x && (1.0 / 2 - 1.0 / 8) * (M - 1) <= y && (1.0 / 2 + 1.0 / 8) * (M - 1) >= y)
                            F[x, y] = LoadValue * Math.Cos(Frequency * t);
                        else
                            F[x, y] = 0;

            ///Квадрат нагрузки в четвертинке равен Q остальное 0
            if (type == 3)
                for (int x = 0; x < N; x++)
                    for (int y = 0; y < M; y++)
                        if ((1.0 / 4 - 1.0 / 8) * (N - 1) <= x && (1.0 / 4 + 1.0 / 8) * (N - 1) >= x && (1.0 / 4 - 1.0 / 8) * (M - 1) <= y && (1.0 / 4 + 1.0 / 8) * (M - 1) >= y)
                            F[x, y] = LoadValue * Math.Cos(Frequency * t);
                        else
                            F[x, y] = 0;
            ///Равномерно распределённая нагрузка
            if (type == 1)
                for (int x = 0; x < N; x++)
                    for (int y = 0; y < M; y++)
                        F[x, y] = LoadValue * Math.Cos(Frequency * t);
        }

        public Load(int N, int M)
        {
            this.N = N;
            this.M = M;
            this.LoadValue = 0;
            this.Frequency = 0;
            this.FlagTimeProblem = false;
            this.type = 1;
            F = new double[N, M];
            FullLoad();
        }
        public Load(int N, int M, int type, double LoadValue)
        {
            this.N = N;
            this.M = M;
            this.LoadValue = LoadValue;
            this.Frequency = 0;
            this.FlagTimeProblem = false;
            if (type == 1 || type == 2 || type == 3)
                this.type = type;
            else
            {
                this.type = 1;
                Console.WriteLine("Error, No true type Load!!!");
            }
            F = new double[N, M];
            FullLoad();
        }
        public Load(int N, int M, int type, double LoadValue, double Frequency)
        {
            this.N = N;
            this.M = M;
            this.LoadValue = LoadValue;
            this.FlagTimeProblem = true;
            if (type == 1 || type == 2 || type == 3)
                this.type = type;
            else
            {
                this.type = 1;
                Console.WriteLine("Error, No true type Load!!!");
            }
            this.Frequency = Frequency;
            F = new double[N, M];
            FullLoad();
        }

    }
}
