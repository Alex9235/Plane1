using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace Decoder
{
    public class PlaneTimoshenco : Plane
    {
        // массив используемых интегралов
        public double[,] Fun { get; }// Интегралл введённой функции связанный со смещением слоёв 
        public double[,] gamx { get; }// Фунция смещения слоёв относительно оси х
        public double[,] gamy { get; }// Фунция смещения слоёв относительно оси y
        // Масссивы для подсчёта переменных констант при системе дифференциального уравнения
        public double[,] DM1 { get; }
        public double[,] DM2 { get; }
        public double[,] RM1 { get; }
        public double[,] RM2 { get; }
        public double[,] SM1 { get; }
        public double[,] SM2 { get; }
        public double[,] IM { get; }// Массив интеграллов с темпретурой 
        public PlaneTimoshenco(int N, int M, int P, double n, double m, double p) : base(N, M, P, n, m, p)//
        {
            this.gamx = new double[N + 2, M + 2];
            this.gamy = new double[N + 2, M + 2];
            this.DM1 = new double[N + 2, M + 2];
            this.RM1 = new double[N + 2, M + 2];
            this.SM1 = new double[N + 2, M + 2];
            this.DM2 = new double[N + 2, M + 2];
            this.RM2 = new double[N + 2, M + 2];
            this.SM2 = new double[N + 2, M + 2];
            this.IM = new double[N + 2, M + 2];
            this.Fun = new double[N, M];
        }
        //Вычисление интенсивности деформаций
        //Вводимая условная функция
        private double fun(double z)
        {
            return 6 * ((1.0 / 4) - z * z);
        }
        private void FullFun()
        {
            double S;

            for (int x = 0; x < N; x++)
                for (int y = 0; y < M; y++)
                {
                    S = 0;
                    //Формула симсона
                    for (int z = 1; z < P; z = z + 2)
                    {
                        double Z_1 = (z - (P - 1) / 2 - 1) * DZ[x, y];
                        double Z0 = (z - (P - 1) / 2) * DZ[x, y];
                        double Z1 = (z - (P - 1) / 2 + 1) * DZ[x, y];
                        S = S + (fun(Z_1) * fun(Z_1) * 2 * (1 + NU[x, y, z - 1]) / (E[x, y, z - 1] + 2 * l * l) +
                            4 * fun(Z0) * fun(Z0) * 2 * (1 + NU[x, y, z]) / (E[x, y, z] + 2 * l * l) +
                            fun(Z1) * fun(Z1) * 2 * (1 + NU[x, y, z + 1]) / (E[x, y, z + 1] + 2 * l * l)) * DZ[x, y] / 3;
                    }
                    Fun[x, y] = 1 / S;
                }
        }

        protected override double ei(int x, int y, int z)
        {
            double Z0 = (z - (P - 1) / 2) * DZ[x, y];
            double DGamxX = (gamx[x + 2, y + 1] - gamx[x, y + 1]) / (2 * dx);
            double DGamyY = (gamy[x + 1, y + 2] - gamy[x + 1, y]) / (2 * dy);
            double DGamxY = (gamx[x + 1, y + 2] - gamx[x + 1, y]) / (2 * dy);
            double DGamyX = (gamy[x + 2, y + 1] - gamy[x, y + 1]) / (2 * dx);
            double dwx = (W[x + 2, y + 1] - W[x, y + 1]) / (2 * dx);
            double dwy = (W[x + 1, y + 2] - W[x + 1, y]) / (2 * dy);
            double d2wx = (W[x, y + 1] - 2 * W[x + 1, y + 1] + W[x + 2, y + 1]) / (dx * dx);
            double d2wy = (W[x + 1, y] - 2 * W[x + 1, y + 1] + W[x + 1, y + 2]) / (dy * dy);
            double d2wxy = (W[x + 2, y + 2] - W[x + 2, y] - W[x, y + 2] + W[x, y]) / (4 * dx * dy);
            double lam0 = lam * lam;
            double lam01 = lam1 * lam1;
            double Tet = 1.0 + NU[x, y, z] / Math.Pow(1 - NU[x, y, z], 2.0);
            if (Z0 != 0)
                if (FlagMultiModularIntensiveDeformations)
                    return (double)(2.0 * Z0 *
                        Math.Sqrt(Tet * Math.Pow(DGamxX - Z0 * Z0 * (DGamxX + d2wx) / 3 + lam0 * (DGamyY - Z0 * Z0 * (DGamyY + d2wy) / 3), 2.0) -
                        3 * lam0 * (DGamxX - Z0 * Z0 * (DGamxX + d2wx) / 3) * (DGamyY - Z0 * Z0 * (DGamyY + d2wy) / 3) +
                        3 * (lam0 * Math.Pow(DGamxY + DGamyX - Z0 * Z0 * (DGamxY + DGamyX + d2wxy) / 3, 2.0) +
                        G0 * G0 * lam0 * lam01 * Math.Pow((gamy[x + 1, y + 1] + dwy) * Fun[x, y] * fun(Z0) / Z0, 2.0) +
                        G0 * G0 * lam01 * Math.Pow((gamx[x + 1, y + 1] + dwx) * Fun[x, y] * fun(Z0) / Z0, 2.0)) / 2) / 3);
                else
                    return (double)(2.0 * Math.Abs(Z0) *
                        Math.Sqrt(Tet * Math.Pow(DGamxX - Z0 * Z0 * (DGamxX + d2wx) / 3 + lam0 * (DGamyY - Z0 * Z0 * (DGamyY + d2wy) / 3), 2.0) -
                        3 * lam0 * (DGamxX - Z0 * Z0 * (DGamxX + d2wx) / 3) * (DGamyY - Z0 * Z0 * (DGamyY + d2wy) / 3) +
                        3 * (lam0 * Math.Pow(DGamxY + DGamyX - Z0 * Z0 * (DGamxY + DGamyX + d2wxy) / 3, 2.0) +
                        G0 * G0 * lam0 * lam01 * Math.Pow((gamy[x + 1, y + 1] + dwy) * Fun[x, y] * fun(Z0) / Z0, 2.0) +
                        G0 * G0 * lam01 * Math.Pow((gamx[x + 1, y + 1] + dwx) * Fun[x, y] * fun(Z0) / Z0, 2.0)) / 2) / 3);
            //) / 2) / 3);
            else
                return 0;
        }
        // Заполнение W Начальными значениями

        private void LoadW()
        {
            for (int x = 0; x < W.GetLength(0); x++)
                for (int y = 0; y < W.GetLength(1); y++)
                    W[x, y] = 0;
        }
        private void LoadGamx()
        {
            for (int x = 0; x < W.GetLength(0); x++)
                for (int y = 0; y < W.GetLength(1); y++)
                    gamx[x, y] = 0;
        }
        private void LoadGamy()
        {
            for (int x = 0; x < W.GetLength(0); x++)
                for (int y = 0; y < W.GetLength(1); y++)
                    gamy[x, y] = 0;
        }
        // расчёт интегралов, метод сведения 3-х мерной задачи к двумерной 
        private void LoadD()
        {
            double S = 0;
            double Z_1 = 0;
            double Z0 = 0;
            double Z1 = 0;

            for (int x = 0; x < N; x++)
                for (int y = 0; y < M; y++)
                {
                    S = 0;
                    //Формула симсона
                    for (int z = 1; z < P; z = z + 2)
                    {
                        Z_1 = (z - (P - 1) / 2 - 1) * DZ[x, y];
                        Z0 = (z - (P - 1) / 2) * DZ[x, y];
                        Z1 = (z - (P - 1) / 2 + 1) * DZ[x, y];
                        S = S + (E[x, y, z - 1] * ((1 - Z_1 * Z_1 / 3) * Z_1 * Z_1 / (1 - NU[x, y, z - 1] * NU[x, y, z - 1]) + l * l / (1 + NU[x, y, z - 1]))
                                     + 4 * E[x, y, z] * ((1 - Z0 * Z0 / 3) * Z0 * Z0 / (1 - NU[x, y, z] * NU[x, y, z]) + l * l / (1 + NU[x, y, z]))
                                     + E[x, y, z + 1] * ((1 - Z1 * Z1 / 3) * Z1 * Z1 / (1 - NU[x, y, z + 1] * NU[x, y, z + 1]) + l * l / (1 + NU[x, y, z + 1]))) * DZ[x, y] / 3;
                    }
                    DM1[x + 1, y + 1] = S;
                    S = 0;
                    //Формула симсона
                    for (int z = 1; z < P; z = z + 2)
                    {
                        Z_1 = (z - (P - 1) / 2 - 1) * DZ[x, y];
                        Z0 = (z - (P - 1) / 2) * DZ[x, y];
                        Z1 = (z - (P - 1) / 2 + 1) * DZ[x, y];
                        S = S + (E[x, y, z - 1] * (Z_1 * Z_1 * Z_1 * Z_1 / (1 - NU[x, y, z - 1] * NU[x, y, z - 1]) + l * l / (1 + NU[x, y, z - 1]))
                                     + 4 * E[x, y, z] * (Z0 * Z0 * Z0 * Z0 / (1 - NU[x, y, z] * NU[x, y, z]) + l * l / (1 + NU[x, y, z]))
                                     + E[x, y, z + 1] * (Z1 * Z1 * Z1 * Z1 / (1 - NU[x, y, z + 1] * NU[x, y, z + 1]) + l * l / (1 + NU[x, y, z + 1]))) * DZ[x, y] / 3;
                    }
                    DM2[x + 1, y + 1] = S;
                }
            //Установление за граничных значений по линейной зависимости
            for (int y = 0; y < M; y++)
            {
                DM1[0, y + 1] = 2 * DM1[1, y + 1] - DM1[2, y + 1];
                DM1[N + 1, y + 1] = 2 * DM1[N, y + 1] - DM1[N - 1, y + 1];
                DM2[0, y + 1] = 2 * DM2[1, y + 1] - DM2[2, y + 1];
                DM2[N + 1, y + 1] = 2 * DM2[N, y + 1] - DM2[N - 1, y + 1];
            }
            for (int x = 0; x < N; x++)
            {
                DM1[x + 1, 0] = 2 * DM1[x + 1, 1] - DM1[x + 1, 2];
                DM1[x + 1, M + 1] = 2 * DM1[x + 1, M] - DM1[x + 1, M - 1];
                DM2[x + 1, 0] = 2 * DM2[x + 1, 1] - DM2[x + 1, 2];
                DM2[x + 1, M + 1] = 2 * DM2[x + 1, M] - DM2[x + 1, M - 1];
            }
            //угловые точки
            DM1[0, 0] = DM1[0, 1] - 0.5 * DM1[0, 2] + DM1[1, 0] - 0.5 * DM1[2, 0];
            DM1[0, M + 1] = DM1[0, M] - 0.5 * DM1[0, M - 1] + DM1[1, M + 1] - 0.5 * DM1[2, M + 1];
            DM1[N + 1, 0] = DM1[N + 1, 1] - 0.5 * DM1[N + 1, 2] + DM1[N, 0] - 0.5 * DM1[N - 1, 0];
            DM1[N + 1, M + 1] = DM1[N + 1, M] - 0.5 * DM1[N + 1, M - 1] + DM1[N, M + 1] - 0.5 * DM1[N - 1, M + 1];
            //угловые точки
            DM2[0, 0] = DM2[0, 1] - 0.5 * DM2[0, 2] + DM2[1, 0] - 0.5 * DM2[2, 0];
            DM2[0, M + 1] = DM2[0, M] - 0.5 * DM2[0, M - 1] + DM2[1, M + 1] - 0.5 * DM2[2, M + 1];
            DM2[N + 1, 0] = DM2[N + 1, 1] - 0.5 * DM2[N + 1, 2] + DM2[N, 0] - 0.5 * DM2[N - 1, 0];
            DM2[N + 1, M + 1] = DM2[N + 1, M] - 0.5 * DM2[N + 1, M - 1] + DM2[N, M + 1] - 0.5 * DM2[N - 1, M + 1];


            //ShowMassiv2(DM);
        }
        private void LoadR()
        {
            double S = 0;
            double Z_1 = 0;
            double Z0 = 0;
            double Z1 = 0;

            for (int x = 0; x < N; x++)
                for (int y = 0; y < M; y++)
                {
                    S = 0;
                    //Формула симсона
                    for (int z = 1; z < P; z = z + 2)
                    {
                        Z_1 = (z - (P - 1) / 2 - 1) * DZ[x, y];
                        Z0 = (z - (P - 1) / 2) * DZ[x, y];
                        Z1 = (z - (P - 1) / 2 + 1) * DZ[x, y];
                        S = S + (E[x, y, z - 1] * NU[x, y, z - 1] * ((1 - Z_1 * Z_1 / 3) * Z_1 * Z_1 / (1 - NU[x, y, z - 1] * NU[x, y, z - 1]) + l * l / (1 + NU[x, y, z - 1]))
                                     + 4 * E[x, y, z] * NU[x, y, z] * ((1 - Z0 * Z0 / 3) * Z0 * Z0 / (1 - NU[x, y, z] * NU[x, y, z]) + l * l / (1 + NU[x, y, z]))
                                     + E[x, y, z + 1] * NU[x, y, z + 1] * ((1 - Z1 * Z1 / 3) * Z1 * Z1 / (1 - NU[x, y, z + 1] * NU[x, y, z + 1]) + l * l / (1 + NU[x, y, z + 1]))) * DZ[x, y] / 3;
                    }
                    RM1[x + 1, y + 1] = S;
                    S = 0;
                    //Формула симсона
                    for (int z = 1; z < P; z = z + 2)
                    {
                        Z_1 = (z - (P - 1) / 2 - 1) * DZ[x, y];
                        Z0 = (z - (P - 1) / 2) * DZ[x, y];
                        Z1 = (z - (P - 1) / 2 + 1) * DZ[x, y];
                        S = S + (E[x, y, z - 1] * NU[x, y, z - 1] * (Z_1 * Z_1 * Z_1 * Z_1 / (1 - NU[x, y, z - 1] * NU[x, y, z - 1]) + l * l / (1 + NU[x, y, z - 1]))
                                     + 4 * E[x, y, z] * NU[x, y, z] * (Z0 * Z0 * Z0 * Z0 / (1 - NU[x, y, z] * NU[x, y, z]) + l * l / (1 + NU[x, y, z]))
                                     + E[x, y, z + 1] * NU[x, y, z + 1] * (Z1 * Z1 * Z1 * Z1 / (1 - NU[x, y, z + 1] * NU[x, y, z + 1]) + l * l / (1 + NU[x, y, z + 1]))) * DZ[x, y] / 3;
                    }
                    RM2[x + 1, y + 1] = S;
                }
            //Установление за граничных значений по линейной зависимости
            for (int y = 0; y < M; y++)
            {
                RM1[0, y + 1] = 2 * RM1[1, y + 1] - RM1[2, y + 1];
                RM1[N + 1, y + 1] = 2 * RM1[N, y + 1] - RM1[N - 1, y + 1];
                RM2[0, y + 1] = 2 * RM2[1, y + 1] - RM2[2, y + 1];
                RM2[N + 1, y + 1] = 2 * RM2[N, y + 1] - RM2[N - 1, y + 1];
            }
            for (int x = 0; x < N; x++)
            {
                RM1[x + 1, 0] = 2 * RM1[x + 1, 1] - RM1[x + 1, 2];
                RM1[x + 1, M + 1] = 2 * RM1[x + 1, M] - RM1[x + 1, M - 1];
                RM2[x + 1, 0] = 2 * RM2[x + 1, 1] - RM2[x + 1, 2];
                RM2[x + 1, M + 1] = 2 * RM2[x + 1, M] - RM2[x + 1, M - 1];
            }
            //угловые точки
            RM1[0, 0] = RM1[0, 1] - 0.5 * RM1[0, 2] + RM1[1, 0] - 0.5 * RM1[2, 0];
            RM1[0, M + 1] = RM1[0, M] - 0.5 * RM1[0, M - 1] + RM1[1, M + 1] - 0.5 * RM1[2, M + 1];
            RM1[N + 1, 0] = RM1[N + 1, 1] - 0.5 * RM1[N + 1, 2] + RM1[N, 0] - 0.5 * RM1[N - 1, 0];
            RM1[N + 1, M + 1] = RM1[N + 1, M] - 0.5 * RM1[N + 1, M - 1] + RM1[N, M + 1] - 0.5 * RM1[N - 1, M + 1];
            //угловые точки
            RM2[0, 0] = RM2[0, 1] - 0.5 * RM2[0, 2] + RM2[1, 0] - 0.5 * RM2[2, 0];
            RM2[0, M + 1] = RM2[0, M] - 0.5 * RM2[0, M - 1] + RM2[1, M + 1] - 0.5 * RM2[2, M + 1];
            RM2[N + 1, 0] = RM2[N + 1, 1] - 0.5 * RM2[N + 1, 2] + RM2[N, 0] - 0.5 * RM2[N - 1, 0];
            RM2[N + 1, M + 1] = RM2[N + 1, M] - 0.5 * RM2[N + 1, M - 1] + RM2[N, M + 1] - 0.5 * RM2[N - 1, M + 1];


        }
        private void LoadS()
        {
            double S = 0;
            double Z_1 = 0;
            double Z0 = 0;
            double Z1 = 0;

            for (int x = 0; x < N; x++)
                for (int y = 0; y < M; y++)
                {
                    S = 0;
                    //Формула симсона
                    for (int z = 1; z < P; z = z + 2)
                    {
                        Z_1 = (z - (P - 1) / 2 - 1) * DZ[x, y];
                        Z0 = (z - (P - 1) / 2) * DZ[x, y];
                        Z1 = (z - (P - 1) / 2 + 1) * DZ[x, y];
                        S = S + (E[x, y, z - 1] * ((1 - Z_1 * Z_1 / 3) * Z_1 * Z_1 + 2 * l * l) / (2 * (1 + NU[x, y, z - 1]))
                                + 4 * E[x, y, z] * ((1 - Z0 * Z0 / 3) * Z0 * Z0 + 2 * l * l) / (2 * (1 + NU[x, y, z]))
                                + E[x, y, z + 1] * ((1 - Z1 * Z1 / 3) * Z1 * Z1 + 2 * l * l) / (2 * (1 + NU[x, y, z]))) * DZ[x, y] / 3;
                    }
                    SM1[x + 1, y + 1] = S;
                    S = 0;
                    //Формула симсона
                    for (int z = 1; z < P; z = z + 2)
                    {
                        Z_1 = (z - (P - 1) / 2 - 1) * DZ[x, y];
                        Z0 = (z - (P - 1) / 2) * DZ[x, y];
                        Z1 = (z - (P - 1) / 2 + 1) * DZ[x, y];
                        S = S + (E[x, y, z - 1] * (Z_1 * Z_1 * Z_1 * Z_1 + 2 * l * l) / (2 * (1 + NU[x, y, z - 1]))
                                + 4 * E[x, y, z] * (Z0 * Z0 * Z0 * Z0 + 2 * l * l) / (2 * (1 + NU[x, y, z]))
                                + E[x, y, z + 1] * (Z1 * Z1 * Z1 * Z1 + 2 * l * l) / (2 * (1 + NU[x, y, z]))) * DZ[x, y] / 3;
                    }
                    SM2[x + 1, y + 1] = S;
                }
            //Установление за граничных значений по линейной зависимости
            for (int y = 0; y < M; y++)
            {
                SM1[0, y + 1] = 2 * SM1[1, y + 1] - SM1[2, y + 1];
                SM1[N + 1, y + 1] = 2 * SM1[N, y + 1] - SM1[N - 1, y + 1];
                SM2[0, y + 1] = 2 * SM2[1, y + 1] - SM2[2, y + 1];
                SM2[N + 1, y + 1] = 2 * SM2[N, y + 1] - SM2[N - 1, y + 1];
            }
            for (int x = 0; x < N; x++)
            {
                SM1[x + 1, 0] = 2 * SM1[x + 1, 1] - SM1[x + 1, 2];
                SM1[x + 1, M + 1] = 2 * SM1[x + 1, M] - SM1[x + 1, M - 1];
                SM2[x + 1, 0] = 2 * SM2[x + 1, 1] - SM2[x + 1, 2];
                SM2[x + 1, M + 1] = 2 * SM2[x + 1, M] - SM2[x + 1, M - 1];
            }
            //угловые точки
            SM1[0, 0] = SM1[0, 1] - 0.5 * SM1[0, 2] + SM1[1, 0] - 0.5 * SM1[2, 0];
            SM1[0, M + 1] = SM1[0, M] - 0.5 * SM1[0, M - 1] + SM1[1, M + 1] - 0.5 * SM1[2, M + 1];
            SM1[N + 1, 0] = SM1[N + 1, 1] - 0.5 * SM1[N + 1, 2] + SM1[N, 0] - 0.5 * SM1[N - 1, 0];
            SM1[N + 1, M + 1] = SM1[N + 1, M] - 0.5 * SM1[N + 1, M - 1] + SM1[N, M + 1] - 0.5 * SM1[N - 1, M + 1];
            //угловые точки
            SM2[0, 0] = SM2[0, 1] - 0.5 * SM2[0, 2] + SM2[1, 0] - 0.5 * SM2[2, 0];
            SM2[0, M + 1] = SM2[0, M] - 0.5 * SM2[0, M - 1] + SM2[1, M + 1] - 0.5 * SM2[2, M + 1];
            SM2[N + 1, 0] = SM2[N + 1, 1] - 0.5 * SM2[N + 1, 2] + SM2[N, 0] - 0.5 * SM2[N - 1, 0];
            SM2[N + 1, M + 1] = SM2[N + 1, M] - 0.5 * SM2[N + 1, M - 1] + SM2[N, M + 1] - 0.5 * SM2[N - 1, M + 1];


        }
        private void LoadI()
        {

            double S = 0;
            double Z_1 = 0;
            double Z0 = 0;
            double Z1 = 0;

            for (int x = 0; x < N; x++)
                for (int y = 0; y < M; y++)
                {
                    S = 0;
                    //Формула симсона
                    for (int z = 1; z < P; z = z + 2)
                    {
                        Z_1 = (z - (P - 1) / 2 - 1) * DZ[x, y];
                        Z0 = (z - (P - 1) / 2) * DZ[x, y];
                        Z1 = (z - (P - 1) / 2 + 1) * DZ[x, y];
                        S = S + (E[x, y, z - 1] * T[x, y, z - 1] * Z_1 / (1 - NU[x, y, z - 1]) + l * l / (1 + NU[x, y, z - 1])
                                     + 4 * E[x, y, z] * T[x, y, z] * Z0 / (1 - NU[x, y, z]) + l * l / (1 + NU[x, y, z])
                                     + E[x, y, z + 1] * T[x, y, z + 1] * Z1 / (1 - NU[x, y, z + 1]) + l * l / (1 + NU[x, y, z + 1])) * DZ[x, y] / 3;
                    }
                    IM[x + 1, y + 1] = S;
                }
            //Установление за граничных значений по линейной зависимости
            for (int y = 0; y < M; y++)
            {
                IM[0, y + 1] = 2 * IM[1, y + 1] - IM[2, y + 1];
                IM[N + 1, y + 1] = 2 * IM[N, y + 1] - IM[N - 1, y + 1];
            }
            for (int x = 0; x < N; x++)
            {
                IM[x + 1, 0] = 2 * IM[x + 1, 1] - IM[x + 1, 2];
                IM[x + 1, M + 1] = 2 * IM[x + 1, M] - IM[x + 1, M - 1];
            }
            //угловые точки
            IM[0, 0] = IM[0, 1] - 0.5 * IM[0, 2] + IM[1, 0] - 0.5 * IM[2, 0];
            IM[0, M + 1] = IM[0, M] - 0.5 * IM[0, M - 1] + IM[1, M + 1] - 0.5 * IM[2, M + 1];
            IM[N + 1, 0] = IM[N + 1, 1] - 0.5 * IM[N + 1, 2] + IM[N, 0] - 0.5 * IM[N - 1, 0];
            IM[N + 1, M + 1] = IM[N + 1, M] - 0.5 * IM[N + 1, M - 1] + IM[N, M + 1] - 0.5 * IM[N - 1, M + 1];
            //ShowMassiv2(IM);
        }


        // загрузка начальных значений
        public void RELoad()
        {

            //обнуление прогибов
            LoadW();
            LoadGamx();
            LoadGamy();
            // для начала забиваем коэффициент пуассона 
            for (int x = 0; x < N; x++)
                for (int y = 0; y < M; y++)
                    for (int z = 0; z < P; z++)
                        NU[x, y, z] = nu;
            //Обнуление деформций
            FULZero(Eii);
            //Установление модулей юнга 
            LoadE();
            if (FlagPorysotyProblem)
            {
                PoritostbLoadNU();
                PoritostbLoadE();
            }
            if (FlagGradientProblem)
            {
                GradientLoadE();
                GradientLoadNU();
            }

            //подсчёт моментов
            if (FlagTemperatureProblem)
                LoadI();

            LoadD();
            LoadR();
            LoadS();
            FullFun();

        }
        // загрузка физ параметров относительно установленного прогиба.
        private void Load()
        {
            //вычисление деформции
            LoadEii();
            //ShowMassiv3CentrR3(Eii, 6);
            //вычисление новых коэффициентов пуассона
            LoadNU();
            //вычисление новых модулей Юнга
            LoadE();
            if (FlagPorysotyProblem)
            {
                PoritostbLoadNU();
                PoritostbLoadE();
            }
            if (FlagGradientProblem)
            {
                GradientLoadE();
                GradientLoadNU();
            }
            //ShowMassiv2(E);
            //подсчёт моментов
            if (FlagTemperatureProblem)
                LoadI();
            LoadD();
            LoadR();
            LoadS();
            FullFun();
        }

        private double

            DM1dx, DM1dy, DM1d2x, DM1d2y,
            DM2dx, DM2dy, DM2d2x, DM2d2y,
            RM1d2x, RM1d2y, RM1dx, RM1dy,
            RM2d2x, RM2d2y, RM2dx, RM2dy,
            SM1dx, SM1dy, SM1dxdy,
            SM2dx, SM2dy, SM2dxdy,
            IMdx, IMd2x, IMdy, IMd2y;


        private int CountIterationInMethods = 0;// Счётчик итерационного метода
        //public int CountIterationABS = 0;
        public int MaxCountIterationSpendOfMetods = 0; // Максимальное количество итераций затриченных методом
        public int MaximumCountIterationOfMethodAllDecision = 0; // Суммарное число всех вложанных друг в друга итераций
        private double BorderExitIterationMethods = 0; //граница выхода метода вариационных итераций
        //private double BorderExitMethodABS = 0; //граница выхода метода Аграновского Баглая Смирнова
        private double BorderExitPhisicalIteration = 0; //граница выхода подсчёта упрого-пластических дифформаций
        public int CountIterationABS = 0;

        //для проверки выхода через максимум для итерационных методов
        private double MaxWForMethodIteration1;
        private double MaxWForMethodIteration2;
        // Делегаты для смены метода расчёта уравнения и исаользования в его физической нелинейности
        delegate void DelegatForMethods(Load F, int TypeBorder1, int TypeBorder2, int TypeBorder3, int TypeBorder4);
        public override void StaticDecisionMethodVariationIteration(Load F, int BorderExitIterationMethods,
            int BorderExitMethodABS, int BorderExitPhisicalIteration, int TypeBorder1, int TypeBorder2, int TypeBorder3, int TypeBorder4)
        {
            RELoad();
            this.BorderExitIterationMethods = BorderExitIterationMethods;
            //this.BorderExitMethodABS = BorderExitMethodABS;
            this.BorderExitPhisicalIteration = BorderExitPhisicalIteration;
            MaxCountIterationSpendOfMetods = 0;
            MaximumCountIterationOfMethodAllDecision = 0;
            if (BorderExitMethodABS == 0)
                PhisicalNolenearyPart(F, TypeBorder1, TypeBorder2, TypeBorder3, TypeBorder4, false, MethodVariationIteration);
            else
                PhisicalNolenearyPart(F, TypeBorder1, TypeBorder2, TypeBorder3, TypeBorder4, true, MethodVariationIteration);
            //Заполнения значиний прогибов для печати
            for (int i = 0; i < N; i++)
                for (int j = 0; j < M; j++)
                    PrintW[i, j] = W[i + 1, j + 1];

        }
        ///Расчёт  пластинки методом вариационных итераций............................

        //Параметры использованнае в методе Вариационных итераций
        // Интегальные функции играющие роль в виде не постоянных коэффициентов связанных 
        //с решение cистемы линейного дифференциального уранвения уравнения  4 порядка

        private double[,,] A; //массив переменных коэфициентов находящисхся при производных
        double[,] Xn;// по х
        double[,] Yn;// по х

        // Интегральные функции при линейном дифференциальном уравнени 
        // Переменные для производных
        private double
            XY0d2x, XY0dx, XY0d2y, XY0dy, XY0d4x, XY0d3x, XY0d4y, XY0d3y,
            XY1d2x, XY1dx, XY1d2y, XY1dy, XY1d4x, XY1d3x, XY1d4y, XY1d3y,
            XY2d2x, XY2dx, XY2d2y, XY2dy, XY2d4x, XY2d3x, XY2d4y, XY2d3y;
        private double
            XY0d, XY0d2, XY0d3, XY1d, XY1d2, XY1d3, XY2d, XY2d2, XY2d3;
        private double Aij(double[,] XY, int x, int y, bool xy, double[,] F, int Ai, int Aj)
        {

            switch ((Ai, Aj))
            {
                case (0, 0):
                    if (xy) return (double)(Math.Pow(XY[0, x + 1], 2.0) * DM2[x, y] * Math.Pow(lam, 4.0) / 3);
                    else return (double)Math.Pow(XY[0, y + 1], 2.0) * DM2[x, y] / 3;
                    break;
                case (0, 1):
                    DM2dx = (DM2[x + 1, y] - DM2[x - 1, y]) / (2 * dx);
                    DM2dy = (DM2[x, y + 1] - DM2[x, y - 1]) / (2 * dy);
                    if (xy) return 2 * (double)(Math.Pow(XY[0, x + 1], 2.0) * DM2dy * Math.Pow(lam, 4.0)) / 3;
                    else return 2 * (double)Math.Pow(XY[0, y + 1], 2.0) * DM2dx / 3;

                case (0, 2):
                    DM2d2x = (DM2[x + 1, y] - 2 * DM2[x, y] + DM2[x - 1, y]) / (dx * dx);
                    DM2d2y = (DM2[x, y + 1] - 2 * DM2[x, y] + DM2[x, y - 1]) / (dy * dy);
                    RM2d2x = (RM2[x + 1, y] - 2 * RM2[x, y] + RM2[x - 1, y]) / (dx * dx);
                    RM2d2y = (RM2[x, y + 1] - 2 * RM2[x, y] + RM2[x, y - 1]) / (dy * dy);
                    RM2dx = (RM2[x + 1, y] - RM2[x - 1, y]) / (2 * dx);
                    RM2dy = (RM2[x, y + 1] - RM2[x, y - 1]) / (2 * dy);
                    SM2dx = (SM2[x + 1, y] - SM2[x - 1, y]) / (2 * dx);
                    SM2dy = (SM2[x, y + 1] - SM2[x, y - 1]) / (2 * dy);
                    if (xy)
                    {
                        XY0d2x = (XY[0, x + 2] - 2 * XY[0, x + 1] + XY[0, x]) / (dx * dx);
                        XY0dx = (XY[0, x + 2] - XY[0, x]) / (2 * dx);
                        return XY[0, x + 1] * ((double)Math.Pow(lam, 4.0) * XY[0, x + 1] * (DM2d2y + RM2d2x / (lam * lam)) +
                            2 * lam * lam * (XY0d2x * (RM2[x, y] + SM2[x, y]) + XY0dx * (RM2dx + SM2dx))) / 3;
                    }
                    else
                    {
                        XY0d2y = (XY[0, y + 2] - 2 * XY[0, y + 1] + XY[0, y]) / (dy * dy);
                        XY0dy = (XY[0, y + 2] - XY[0, y]) / (2 * dy);
                        return XY[0, y + 1] * (XY[0, y + 1] * (DM2d2x + lam * lam * RM2d2y) +
                            2 * lam * lam * (XY0d2y * (RM2[x, y] + SM2[x, y]) + XY0dy * (RM2dy + SM2dy))) / 3;
                    }

                case (0, 3):
                    RM2dx = (RM2[x + 1, y] - RM2[x - 1, y]) / (2 * dx);
                    RM2dy = (RM2[x, y + 1] - RM2[x, y - 1]) / (2 * dy);
                    SM2dx = (SM2[x + 1, y] - SM2[x - 1, y]) / (2 * dx);
                    SM2dy = (SM2[x, y + 1] - SM2[x, y - 1]) / (2 * dy);
                    SM2dxdy = (SM2[x + 1, y + 1] - SM2[x + 1, y - 1] - SM2[x - 1, y + 1] + SM2[x - 1, y - 1]) / (4 * dx * dy);
                    if (xy)
                    {
                        XY0d2x = (XY[0, x + 2] - 2 * XY[0, x + 1] + XY[0, x]) / (dx * dx);
                        XY0dx = (XY[0, x + 2] - XY[0, x]) / (2 * dx);
                        return 2 * lam * lam * XY[0, x + 1] * (XY0d2x * (RM2dy + SM2dy) + XY0dx * SM2dxdy) / 3;
                    }
                    else
                    {
                        XY0d2y = (XY[0, y + 2] - 2 * XY[0, y + 1] + XY[0, y]) / (dy * dy);
                        XY0dy = (XY[0, y + 2] - XY[0, y]) / (2 * dy);
                        return 2 * lam * lam * XY[0, y + 1] * (XY0d2y * (RM2dx + SM2dx) + XY0dy * SM2dxdy) / 3;
                    }

                case (0, 4):
                    DM2d2x = (DM2[x + 1, y] - 2 * DM2[x, y] + DM2[x - 1, y]) / (dx * dx);
                    DM2d2y = (DM2[x, y + 1] - 2 * DM2[x, y] + DM2[x, y - 1]) / (dy * dy);
                    DM2dx = (DM2[x + 1, y] - DM2[x - 1, y]) / (2 * dx);
                    DM2dy = (DM2[x, y + 1] - DM2[x, y - 1]) / (2 * dy);
                    RM2d2x = (RM2[x + 1, y] - 2 * RM2[x, y] + RM2[x - 1, y]) / (dx * dx);
                    RM2d2y = (RM2[x, y + 1] - 2 * RM2[x, y] + RM2[x, y - 1]) / (dy * dy);

                    if (xy)
                    {
                        XY0d4x = (XY[0, x + 3] - 4 * XY[0, x + 2] + 6 * XY[0, x + 1] - 4 * XY[0, x] + XY[0, x - 1]) / (dx * dx * dx * dx);
                        XY0d3x = (XY[0, x + 3] - 2 * XY[0, x + 2] + 2 * XY[0, x] - XY[0, x - 1]) / (2 * dx * dx * dx);
                        XY0d2x = (XY[0, x + 2] - 2 * XY[0, x + 1] + XY[0, x]) / (dx * dx);
                        return XY[0, x + 1] * (XY0d4x * DM2[x, y] + 2 * XY0d3x * DM2dx + XY0d2x * (DM2d2x + RM2d2y * lam * lam)) / 3;
                    }
                    else
                    {
                        XY0d4y = (XY[0, y + 3] - 4 * XY[0, y + 2] + 6 * XY[0, y + 1] - 4 * XY[0, y] + XY[0, y - 1]) / (dy * dy * dy * dy);
                        XY0d3y = (XY[0, y + 3] - 2 * XY[0, y + 2] + 2 * XY[0, y] - XY[0, y - 1]) / (2 * dy * dy * dy);
                        XY0d2y = (XY[0, y + 2] - 2 * XY[0, y + 1] + XY[0, y]) / (dy * dy);
                        return (double)Math.Pow(lam, 4.0) * XY[0, y + 1] * (XY0d4y * DM2[x, y] + 2 * XY0d3y * DM2dy +
                            XY0d2y * (DM2d2y + RM2d2x / (lam * lam))) / 3;
                    }

                case (0, 5):
                    if (xy) return 0;
                    else
                    {
                        return XY[0, y + 1] * XY[1, y + 1] * DM1[x, y];
                    }

                case (0, 6):
                    DM1dx = (DM1[x + 1, y] - DM1[x - 1, y]) / (2 * dx);
                    SM1dx = (SM1[x + 1, y] - SM1[x - 1, y]) / (2 * dx);
                    if (xy)
                    {
                        XY1dx = (XY[1, x + 2] - XY[1, x]) / (2 * dx);
                        return lam * lam * XY[0, x + 1] * (XY1dx * (RM1[x, y] + SM1[x, y]) + XY[1, x + 1] * SM1dx);
                    }
                    else
                    {
                        return 2 * XY[0, y + 1] * XY[1, y + 1] * DM1dx;
                    }

                case (0, 7):
                    DM1d2x = (DM1[x + 1, y] - 2 * DM1[x, y] + DM1[x - 1, y]) / (dx * dx);
                    RM1d2y = (RM1[x, y + 1] - 2 * RM1[x, y] + RM1[x, y - 1]) / (dy * dy);
                    RM1dy = (RM1[x, y + 1] - RM1[x, y - 1]) / (2 * dy);
                    SM1dy = (SM1[x, y + 1] - SM1[x, y - 1]) / (2 * dy);
                    SM1dxdy = (SM1[x + 1, y + 1] - SM1[x + 1, y - 1] - SM1[x - 1, y + 1] + SM1[x - 1, y - 1]) / (4 * dx * dy);
                    if (xy)
                    {
                        XY1dx = (XY[1, x + 2] - XY[1, x]) / (2 * dx);
                        return lam * lam * XY[0, x + 1] * (XY1dx * (2 * RM1dy + SM1dy) + XY[1, x + 1] * SM1dxdy);
                    }
                    else
                    {
                        XY1d2y = (XY[1, y + 2] - 2 * XY[1, y + 1] + XY[1, y]) / (dy * dy);
                        XY1dy = (XY[1, y + 2] - XY[1, y]) / (2 * dy);
                        return XY[0, y + 1] * (XY[1, y + 1] * (DM1d2x + lam * lam * RM1d2y) + lam * lam * (XY1d2y * (RM1[x, y] + SM1[x, y]) +
                            XY1dy * (2 * RM1dy + SM1dy)));
                    }

                case (0, 8):
                    DM1d2x = (DM1[x + 1, y] - 2 * DM1[x, y] + DM1[x - 1, y]) / (dx * dx);
                    DM1dx = (DM1[x + 1, y] - DM1[x - 1, y]) / (2 * dx);
                    RM1d2y = (RM1[x, y + 1] - 2 * RM1[x, y] + RM1[x, y - 1]) / (dy * dy);
                    SM1dx = (SM1[x + 1, y] - SM1[x - 1, y]) / (2 * dx);
                    SM1dxdy = (SM1[x + 1, y + 1] - SM1[x + 1, y - 1] - SM1[x - 1, y + 1] + SM1[x - 1, y - 1]) / (4 * dx * dy);
                    if (xy)
                    {
                        XY1d3x = (XY[1, x + 3] - 2 * XY[1, x + 2] + 2 * XY[1, x] - XY[1, x - 1]) / (2 * dx * dx * dx);
                        XY1d2x = (XY[1, x + 2] - 2 * XY[1, x + 1] + XY[1, x]) / (dx * dx);
                        XY1dx = (XY[1, x + 2] - XY[1, x]) / (2 * dx);
                        return XY[0, x + 1] * (XY1d3x * DM1[x, y] + 2 * XY1d2x * DM1dx + XY1dx * (DM1d2x + lam * lam * RM1d2y));
                    }
                    else
                    {
                        XY1d2y = (XY[1, y + 2] - 2 * XY[1, y + 1] + XY[1, y]) / (dy * dy);
                        XY1dy = (XY[1, y + 2] - XY[1, y]) / (2 * dy);
                        return XY[0, y + 1] * lam * lam * (XY1d2y * SM1dx + XY1dy * SM1dxdy);
                    }

                case (0, 9):
                    if (xy)
                    {
                        return (double)Math.Pow(lam, 6.0) * XY[2, x + 1] * XY[0, x + 1] * DM1[x, y];
                    }
                    else
                    {
                        return 0;
                    }

                case (0, 10):
                    DM1dy = (DM1[x, y + 1] - DM1[x, y - 1]) / (2 * dy);
                    SM1dy = (SM1[x, y + 1] - SM1[x, y - 1]) / (2 * dy);
                    if (xy)
                    {
                        return XY[0, x + 1] * (double)Math.Pow(lam, 6.0) * 2 * XY[2, x + 1] * DM1dy;
                    }
                    else
                    {
                        XY2dy = (XY[2, y + 2] - XY[2, y]) / (2 * dy);
                        return lam * lam * XY[0, y + 1] * (XY2dy * (RM1[x, y] + SM1[x, y]) + XY[2, y + 1] * SM1dy);
                    }

                case (0, 11):
                    DM1d2y = (DM1[x, y + 1] - 2 * DM1[x, y] + DM1[x, y - 1]) / (dy * dy);
                    RM1d2x = (RM1[x + 1, y] - 2 * RM1[x, y] + RM1[x - 1, y]) / (dx * dx);
                    RM1dx = (RM1[x + 1, y] - RM1[x - 1, y]) / (2 * dx);
                    SM1dx = (SM1[x + 1, y] - SM1[x - 1, y]) / (2 * dx);
                    SM1dxdy = (SM1[x + 1, y + 1] - SM1[x + 1, y - 1] - SM1[x - 1, y + 1] + SM1[x - 1, y - 1]) / (4 * dx * dy);
                    if (xy)
                    {
                        XY2d2x = (XY[2, x + 2] - 2 * XY[2, x + 1] + XY[2, x]) / (dx * dx);
                        XY2dx = (XY[2, x + 2] - XY[2, x]) / (2 * dx);
                        return lam * lam * XY[0, x + 1] * (lam * lam * lam * lam * XY[2, x + 1] * (DM1d2y + lam * lam * RM1d2x) + XY2d2x * (RM1[x, y] + SM1[x, y]) +
                            XY2dx * (2 * RM1dx + SM1dx));
                    }
                    else
                    {
                        XY2dy = (XY[2, y + 2] - XY[2, y]) / (dy * dy);
                        return lam * lam * XY[0, y + 1] * (XY2dy * (2 * RM1dx + SM1dx) + XY[2, y + 1] * SM1dxdy);
                    }

                case (0, 12):
                    DM1d2y = (DM1[x, y + 1] - 2 * DM1[x, y] + DM1[x, y - 1]) / (dy * dy);
                    DM1dy = (DM1[x, y + 1] - DM1[x, y - 1]) / (2 * dy);
                    RM1d2x = (RM1[x + 1, y] - 2 * RM1[x, y] + RM1[x - 1, y]) / (dx * dx);
                    SM1dy = (SM1[x, y + 1] - SM1[x, y - 1]) / (2 * dy);
                    SM1dxdy = (SM1[x + 1, y + 1] - SM1[x + 1, y - 1] - SM1[x - 1, y + 1] + SM1[x - 1, y - 1]) / (4 * dx * dy);
                    if (xy)
                    {
                        XY2d2x = (XY[2, x + 2] - 2 * XY[2, x + 1] + XY[2, x]) / (dx * dx);
                        XY2dx = (XY[2, x + 2] - XY[1, x]) / (2 * dx);
                        return XY[0, y + 1] * lam * lam * (XY2d2x * SM1dy + XY2dx * SM1dxdy);
                    }
                    else
                    {
                        XY2d3y = (XY[2, y + 3] - 2 * XY[2, y + 2] + 2 * XY[2, y] - XY[2, y - 1]) / (2 * dy * dy * dy);
                        XY2d2y = (XY[2, y + 2] - 2 * XY[2, y + 1] + XY[2, y]) / (dy * dy);
                        XY2dy = (XY[2, y + 2] - XY[2, y]) / (2 * dy);
                        return Math.Pow(lam, 6.0) * XY[0, y + 1] * (XY2d3y * DM1[x, y] + 2 * XY2d2y * DM1dy + XY2dy * (DM1d2y + RM1d2x / (lam * lam)));
                    }

                case (0, 13):
                    if (FlagTemperatureProblem)
                    {
                        IMd2x = (IM[x + 1, y] - 2 * IM[x, y] + IM[x - 1, y]) / (dx * dx);
                        IMd2y = (IM[x, y + 1] - 2 * IM[x, y] + IM[x, y - 1]) / (dy * dy);
                    }
                    if (xy)
                    {
                        return (F[x - 1, y - 1] - IMd2x - lam * lam * IMd2y) * XY[0, x + 1];
                    }
                    else
                    {
                        return (F[x - 1, y - 1] - IMd2x - lam * lam * IMd2y) * XY[0, y + 1];
                    }

                case (1, 0):
                    if (xy)
                    {
                        return 0;
                    }
                    else
                    {
                        return XY[1, y + 1] * XY[0, y + 1] * DM2[x, y] / 3;
                    }

                case (1, 1):
                    DM2dx = (DM2[x + 1, y] - DM2[x - 1, y]) / (2 * dx);
                    RM2dx = (RM2[x + 1, y] - RM2[x - 1, y]) / (2 * dx);
                    if (xy)
                    {
                        XY0dx = (XY[0, x + 2] - XY[0, x]) / (2 * dx);
                        return XY[1, x + 1] * lam * lam * (XY0dx * (SM2[x, y] + RM2[x, y]) + XY[0, x + 1] * RM2dx) / 3;
                    }
                    else
                    {
                        return XY[1, y + 1] * XY[0, y + 1] * DM2dx / 3;
                    }

                case (1, 2):
                    SM2dy = (SM2[x, y + 1] - SM2[x, y - 1]) / (2 * dy);
                    if (xy)
                    {
                        XY0dx = (XY[0, x + 2] - XY[0, x]) / (2 * dx);
                        return XY[1, x + 1] * lam * lam * XY0dx * SM2dy / 3;
                    }
                    else
                    {
                        XY0d2y = (XY[0, y + 2] - 2 * XY[0, y + 1] + XY[0, y]) / (dy * dy);
                        XY0dy = (XY[0, y + 2] - XY[0, y]) / (2 * dy);
                        return XY[1, y + 1] * (lam * lam * (XY0d2y * (SM2[x, y] + RM2[x, y]) + XY0dy * SM2dy) / 3 +
                            lam1 * lam1 * Fun[x - 1, y - 1] * XY[0, y + 1]);
                    }

                case (1, 3):
                    RM2dx = (RM2[x + 1, y] - RM2[x - 1, y]) / (2 * dx);
                    DM2dx = (DM2[x + 1, y] - DM2[x - 1, y]) / (2 * dx);
                    if (xy)
                    {
                        XY0d3x = (XY[0, x + 3] - 2 * XY[0, x + 2] + 2 * XY[0, x] - XY[0, x - 1]) / (2 * dx * dx * dx);
                        XY0d2x = (XY[0, x + 2] - 2 * XY[0, x + 1] + XY[0, x]) / (dx * dx);
                        XY0dx = (XY[0, x + 2] - XY[0, x]) / (2 * dx);
                        return XY[1, x + 1] * ((XY0d3x * DM2[x, y] + XY0d2x * DM2dx) / 3 +
                            lam1 * lam1 * Fun[x - 1, y - 1] * XY0dx);
                    }
                    else
                    {
                        XY0d2y = (XY[0, y + 2] - 2 * XY[0, y + 1] + XY[0, y]) / (dy * dy);
                        return XY[1, y + 1] * lam * lam * XY0d2y * RM2dx / 3;
                    }

                case (1, 4):
                    if (xy)
                    {
                        return XY[1, x + 1] * lam * lam * XY[1, x + 1] * SM1[x, y];
                    }
                    else
                    {
                        return XY[1, y + 1] * XY[1, y + 1] * DM1[x, y];
                    }

                case (1, 5):
                    DM1dx = (DM1[x + 1, y] - DM1[x - 1, y]) / (2 * dx);
                    SM1dy = (SM1[x, y + 1] - SM1[x, y - 1]) / (2 * dy);
                    if (xy)
                    {
                        return XY[1, x + 1] * lam * lam * XY[1, x + 1] * SM1dy;
                    }
                    else
                    {
                        return XY[1, y + 1] * XY[1, y + 1] * DM1dx;
                    }

                case (1, 6):
                    DM1dx = (DM1[x + 1, y] - DM1[x - 1, y]) / (2 * dx);
                    SM1dy = (SM1[x, y + 1] - SM1[x, y - 1]) / (2 * dy);
                    if (xy)
                    {
                        XY1d2x = (XY[1, x + 2] - 2 * XY[1, x + 1] + XY[1, x]) / (dx * dx);
                        XY1dx = (XY[1, x + 2] - XY[1, x]) / (2 * dx);
                        return XY[1, x + 1] * (XY1d2x * DM1[x, y] + XY1dx * DM1dx
                            - lam1 * lam1 * Fun[x - 1, y - 1] * XY[1, x + 1]);
                    }
                    else
                    {
                        XY1d2y = (XY[1, y + 2] - 2 * XY[1, y + 1] + XY[1, y]) / (dy * dy);
                        XY1dy = (XY[1, y + 2] - XY[1, y]) / (2 * dy);
                        return XY[1, y + 1] * (lam * lam * (XY1d2y * SM1[x, y] + XY1dy * SM1dy)
                            - lam1 * lam1 * Fun[x - 1, y - 1] * XY[1, y + 1]);
                    }

                case (1, 7):
                    SM1dy = (SM1[x, y + 1] - SM1[x, y - 1]) / (2 * dy);
                    RM1dx = (RM1[x + 1, y] - RM1[x - 1, y]) / (2 * dx);
                    if (xy)
                    {
                        XY2dx = (XY[2, x + 2] - XY[2, x]) / (2 * dx);
                        return XY[1, x + 1] * lam * lam * (XY2dx * (RM1[x, y] + SM1[x, y]) + XY[2, x + 1] * RM1dx);

                    }
                    else
                    {
                        XY2dy = (XY[2, y + 2] - XY[2, y]) / (2 * dy);
                        return XY[1, y + 1] * lam * lam * (XY2dy * (RM1[x, y] + SM1[x, y]) + XY[2, y + 1] * SM1dy);

                    }
                case (1, 8):
                    SM1dy = (SM1[x, y + 1] - SM1[x, y - 1]) / (2 * dy);
                    RM1dx = (RM1[x + 1, y] - RM1[x - 1, y]) / (2 * dx);
                    if (xy)
                    {
                        XY2dx = (XY[2, x + 2] - XY[2, x]) / (2 * dx);
                        return XY[1, x + 1] * lam * lam * XY2dx * SM1dy;

                    }
                    else
                    {
                        XY2dy = (XY[2, y + 2] - XY[2, y]) / (2 * dy);
                        return XY[1, y + 1] * lam * lam * XY2dy * RM1dx;
                    }
                case (1, 9):
                    IMdx = (IM[x + 1, y] - IM[x - 1, y]) / (2 * dx);
                    if (xy)
                    {
                        return XY[1, x + 1] * IMdx;
                    }
                    else
                    {
                        return XY[1, y + 1] * IMdx;
                    }
                case (2, 0):
                    if (xy)
                    {
                        return XY[2, x + 1] * XY[0, x + 1] * DM2[x, y] / 3;
                    }
                    else
                    {
                        return 0;
                    }

                case (2, 1):
                    DM2dy = (DM2[x, y + 1] - DM2[x, y - 1]) / (2 * dy);
                    RM2dy = (RM2[x, y + 1] - RM2[x, y - 1]) / (2 * dy);
                    if (xy)
                    {
                        return XY[2, x + 1] * XY[0, x + 1] * DM2dy / 3;
                    }
                    else
                    {
                        XY0dy = (XY[0, y + 2] - XY[0, y]) / (2 * dy);
                        return XY[2, y + 1] / (lam * lam) * (XY0dy * (SM2[x, y] + RM2[x, y]) + XY[0, y + 1] * RM2dy) / 3;
                    }

                case (2, 2):
                    SM2dx = (SM2[x + 1, y] - SM2[x - 1, y]) / (2 * dx);
                    if (xy)
                    {
                        XY0d2x = (XY[0, x + 2] - 2 * XY[0, x + 1] + XY[0, x]) / (dx * dx);
                        XY0dx = (XY[0, x + 2] - XY[0, x]) / (2 * dx);
                        return XY[2, x + 1] * ((XY0d2x * (SM2[x, y] + RM2[x, y]) + XY0dx * SM2dx) / (3 * lam * lam) +
                            lam2 * lam2 * Fun[x - 1, y - 1] * XY[0, x + 1]);
                    }
                    else
                    {
                        XY0dy = (XY[0, y + 2] - XY[0, y]) / (2 * dy);
                        return XY[2, y + 1] / (lam * lam) * XY0dy * SM2dx / 3;
                    }

                case (2, 3):
                    RM2dy = (RM2[x, y + 1] - RM2[x, y - 1]) / (2 * dy);
                    DM2dy = (DM2[x, y + 1] - DM2[x, y - 1]) / (2 * dy);
                    if (xy)
                    {
                        XY0d2x = (XY[0, x + 2] - 2 * XY[0, x + 1] + XY[0, x]) / (dx * dx);
                        return XY[2, x + 1] / (lam * lam) * XY0d2x * RM2dy / 3;
                    }
                    else
                    {
                        XY0d3y = (XY[0, y + 3] - 2 * XY[0, y + 2] + 2 * XY[0, y] - XY[0, y - 1]) / (2 * dy * dy * dy);
                        XY0d2y = (XY[0, y + 2] - 2 * XY[0, y + 1] + XY[0, y]) / (dy * dy);
                        XY0dy = (XY[0, y + 2] - XY[0, y]) / (2 * dy);
                        return XY[2, y + 1] * ((XY0d3y * DM2[x, y] + XY0d2y * DM2dy) / 3 +
                            lam2 * lam2 * Fun[x - 1, y - 1] * XY0dy);
                    }

                case (2, 4):
                    if (xy)
                    {
                        return XY[2, x + 1] * XY[2, x + 1] * DM1[x, y];
                    }
                    else
                    {
                        return XY[2, y + 1] / (lam * lam) * XY[2, y + 1] * SM1[x, y];
                    }

                case (2, 5):
                    DM1dy = (DM1[x, y + 1] - DM1[x, y - 1]) / (2 * dy);
                    SM1dx = (SM1[x + 1, y] - SM1[x - 1, y]) / (2 * dx);
                    if (xy)
                    {
                        return XY[2, x + 1] * XY[2, x + 1] * DM1dy;
                    }
                    else
                    {
                        return XY[2, x + 1] / (lam * lam) * XY[2, x + 1] * SM1dx;

                    }

                case (2, 6):
                    DM1dy = (DM1[x, y + 1] - DM1[x, y - 1]) / (2 * dy);
                    SM1dx = (SM1[x + 1, y] - SM1[x - 1, y]) / (2 * dx);
                    if (xy)
                    {
                        XY2d2x = (XY[2, x + 2] - 2 * XY[2, x + 1] + XY[2, x]) / (dx * dx);
                        XY2dx = (XY[2, x + 2] - XY[2, x]) / (2 * dx);
                        return XY[2, x + 1] * ((XY2d2x * SM1[x, y] + XY2dx * SM1dx) / (lam * lam)
                            - lam2 * lam2 * Fun[x - 1, y - 1] * XY[2, x + 1]);
                    }
                    else
                    {
                        XY2d2y = (XY[2, y + 2] - 2 * XY[2, y + 1] + XY[2, y]) / (dy * dy);
                        XY2dy = (XY[2, y + 2] - XY[2, y]) / (2 * dy);
                        return XY[2, y + 1] * (XY2d2y * DM1[x, y] + XY2dy * DM1dy
                            - lam2 * lam2 * Fun[x - 1, y - 1] * XY[2, y + 1]);
                    }

                case (2, 7):
                    SM1dx = (SM1[x + 1, y] - SM1[x - 1, y]) / (2 * dx);
                    RM1dy = (RM1[x, y + 1] - RM1[x, y - 1]) / (2 * dy);
                    if (xy)
                    {
                        XY2dx = (XY[1, x + 2] - XY[1, x]) / (2 * dx);
                        return XY[2, x + 1] / (lam * lam) * (XY2dx * (RM1[x, y] + SM1[x, y]) + XY[1, x + 1] * SM1dx);

                    }
                    else
                    {
                        XY2dy = (XY[1, y + 2] - XY[1, y]) / (2 * dy);
                        return XY[2, y + 1] / (lam * lam) * (XY2dy * (RM1[x, y] + SM1[x, y]) + XY[1, y + 1] * RM1dy);

                    }
                case (2, 8):
                    SM1dx = (SM1[x + 1, y] - SM1[x - 1, y]) / (2 * dx);
                    RM1dy = (RM1[x, y + 1] - RM1[x, y - 1]) / (2 * dy);
                    if (xy)
                    {
                        XY2dx = (XY[1, x + 2] - XY[1, x]) / (2 * dx);
                        return XY[2, x + 1] / (lam * lam) * XY2dx * RM1dy;

                    }
                    else
                    {
                        XY2dy = (XY[1, y + 2] - XY[1, y]) / (2 * dy);
                        return XY[2, y + 1] / (lam * lam) * XY2dy * SM1dx;
                    }
                case (2, 9):
                    IMdy = (IM[x, y + 1] - IM[x, y - 1]) / (2 * dy);
                    if (xy)
                    {
                        return XY[2, x + 1] * IMdy;
                    }
                    else
                    {
                        return XY[2, y + 1] * IMdy;
                    }
                default: return 0;

            }
        }
        private void Pn(double[,,] AXnYn, double[,] XnYn, bool xy, double[,] F)
        {
            // расчёт функций стоящих при роизводных ЛДУ, 
            //в момент нахождения функциии по одной переменной

            // Размерность Массивов N и M
            int iLength = N;
            int jLength = M;
            double f1, f2, f3;
            // Зависит от того по какой переменно идёт находжения функций Pn
            if (xy) { iLength = M; jLength = N; }
            for (int Ai = 0; Ai < A.GetLength(0); Ai++)
            {
                for (int Aj = 0; Aj < A.GetLength(1); Aj++)
                {

                    for (int i = 0; i < iLength - 2; i++)
                    {
                        // обнуление массива
                        AXnYn[Ai, Aj, i] = 0;
                        for (int j = 1; j < jLength; j = j + 2)
                        {
                            //Интегралы считаются на основе формулы симсона
                            if (xy)
                            {
                                // XnYn-это один из набора функций X1(x),X2(x),X3(x) или Y1(y),Y2(y),Y3(y)
                                // Вариация относительно Y1(y),Y2(y),Y3(y)
                                // интегрование по х
                                if (j == jLength - 1)
                                {
                                    //формула трапеции на случай если количество точек разбиение чётное
                                    f1 = Aij(XnYn, j, i + 1, xy, F, Ai, Aj);
                                    f2 = Aij(XnYn, j + 1, i + 1, xy, F, Ai, Aj);
                                    AXnYn[Ai, Aj, i] = A[Ai, Aj, i] + (f1 + f2) * dx / 2;
                                }
                                else
                                {
                                    f1 = Aij(XnYn, j, i + 1, xy, F, Ai, Aj);
                                    f2 = Aij(XnYn, j + 1, i + 1, xy, F, Ai, Aj);
                                    f3 = Aij(XnYn, j + 2, i + 1, xy, F, Ai, Aj);
                                    AXnYn[Ai, Aj, i] = A[Ai, Aj, i] + (f1 + 4 * f2 + f3) * dx / 3;
                                }
                            }
                            else
                            {
                                // Вариация относительно A(x) 
                                // интегрование по у
                                if (j == jLength - 1)
                                {
                                    //формула трапеции на случай если количество точек разбиение чётное
                                    f1 = Aij(XnYn, i + 1, j, xy, F, Ai, Aj);
                                    f2 = Aij(XnYn, i + 1, j + 1, xy, F, Ai, Aj);
                                    AXnYn[Ai, Aj, i] = A[Ai, Aj, i] + (f1 + f2) * dx / 2;
                                }
                                else
                                {
                                    f1 = Aij(XnYn, i + 1, j, xy, F, Ai, Aj);
                                    f2 = Aij(XnYn, i + 1, j + 1, xy, F, Ai, Aj);
                                    f3 = Aij(XnYn, i + 1, j + 2, xy, F, Ai, Aj);
                                    AXnYn[Ai, Aj, i] = A[Ai, Aj, i] + (f1 + 4 * f2 + f3) * dx / 3;
                                }
                            }
                        }
                    }
                }

            }

        }
        private double[,] SolutionSystemsProblemMetodGaussaForVariationIteration(double[,] XY, double d, double[,,] A, int GrUsl1, int GrUsl2, bool xy, double eps)
        {
            //Грачиные условия типа известных значений функции с двух концов 
            //p1(x)y''''(x)+p2(x)y'''(x)+p3(x)y''(x)+p4(x)y'(x)+p5(x)A(x)=f(x)
            //GrUsl=0 - (по умолчанию )жёсткая заделка на граница первая производная равна нулю
            //GrUsl=1 - шарнирное операние на границе вторая производная равна нулю
            //d - шаг разбиения
            //Размерность решения зависит от подаваемых функций. если размерно P функций N, то выходной будет N+4, 
            //так как 4 производная и функции P на границах не учитываются
            double MaxDifferents = 0;

            //Определяем размерность
            int Ra = A.GetLength(2);
            //задаём массив решение 
            //0 1(г) 2 .... Ra+1 Ra+2(г) Ra+3 Cчёт от (2) до (Ra+1) 
            //создаём копию
            double[,] FXnYn = new double[3, XY.GetLength(1)];
            for (int i = 0; i < FXnYn.GetLength(0); i++)
                for (int j = 0; j < FXnYn.GetLength(1); j++)
                    FXnYn[i, j] = XY[i, j];
            //установка по какой переменной идёт счёт 
            //определяем длину
            double Lxy;
            if (xy) Lxy = m;
            else Lxy = n;

            double[] CopyXY = new double[XY.GetLength(1)];
            // матрица для метода Гаусса
            double[,] Gauss = new double[Ra, Ra + 1];

            while (true)
            {

                //Заполнение матрицы Гаусса
                //Создаётся 5 диагональная матрица
                //Решение относительно W 
                for (int i = 0; i < Ra; i++)
                {
                    XY1d = (XY[1, i + 1 + 3] - XY[1, i - 1 + 3]) / (2 * d);
                    XY1d2 = (XY[1, i + 1 + 3] - 2 * XY[1, i + 3] + XY[1, i - 1 + 3]) / (d * d);
                    XY1d3 = (XY[1, i + 2 + 3] - 2 * XY[1, i + 1 + 3] + 2 * XY[1, i - 1 + 3] - XY[1, i - 2 + 3]) / (2 * d * d * d);
                    XY2d = (XY[2, i + 1 + 3] - XY[2, i - 1 + 3]) / (2 * d);
                    XY2d2 = (XY[2, i + 1 + 3] - 2 * XY[2, i + 3] + XY[2, i - 1 + 3]) / (d * d);
                    XY2d3 = (XY[2, i + 2 + 3] - 2 * XY[2, i + 1 + 3] + 2 * XY[2, i - 1 + 3] - XY[2, i - 2 + 3]) / (2 * d * d * d);
                    if (i > 1 && i < Ra - 2)
                    {
                        Gauss[i, i + 2] = A[0, 0, i] / (d * d * d * d) + A[0, 1, i] / (2 * d * d * d);
                        Gauss[i, i + 1] = A[0, 2, i] / (d * d) + A[0, 3, i] / (2 * d) - 4 * A[0, 0, i] / (d * d * d * d)
                            - A[0, 1, i] / (d * d * d);
                        Gauss[i, i] = 6 * A[0, 0, i] / (d * d * d * d) - 2 * A[0, 2, i] / (d * d) + A[0, 4, i];
                        Gauss[i, i - 1] = A[0, 2, i] / (d * d) - A[0, 3, i] / (2 * d) - 4 * A[0, 0, i] / (d * d * d * d)
                            + A[0, 1, i] / (d * d * d);
                        Gauss[i, i - 2] = A[0, 0, i] / (d * d * d * d) - A[0, 1, i] / (2 * d * d * d);
                    }
                    Gauss[i, Ra] = A[0, 5, i] * XY1d3 + A[0, 6, i] * XY1d2 + A[0, 7, i] * XY1d + A[0, 8, i] * XY[1, i + 3]
                        + A[0, 9, i] * XY2d3 + A[0, 10, i] * XY2d2 + A[0, 11, i] * XY2d + A[0, 12, i] * XY[2, i + 3] + A[0, 13, i];
                    //Gauss[i, Ra] = A[0, 13, i];

                    if (i == 0)
                    {
                        if (GrUsl1 == 0)
                        {
                            Gauss[i, i + 2] = A[0, 0, i] / (d * d * d * d) + A[0, 1, i] / (2 * d * d * d);
                            Gauss[i, i + 1] = A[0, 2, i] / (d * d) + A[0, 3, i] / (2 * d) - 4 * A[0, 0, i] / (d * d * d * d)
                                - A[0, 1, i] / (d * d * d);
                            Gauss[i, i] = 6 * A[0, 0, i] / (d * d * d * d) - 2 * A[0, 2, i] / (d * d) + A[0, 4, i]
                                + A[0, 0, i] / (d * d * d * d) - A[0, 1, i] / (2 * d * d * d);
                        }
                        if (GrUsl1 != 0)
                        {
                            Gauss[i, i + 2] = A[0, 0, i] / (d * d * d * d) + A[0, 1, i] / (2 * d * d * d);
                            Gauss[i, i + 1] = A[0, 2, i] / (d * d) + A[0, 3, i] / (2 * d) - 4 * A[0, 0, i] / (d * d * d * d)
                                - A[0, 1, i] / (d * d * d);
                            Gauss[i, i] = 6 * A[0, 0, i] / (d * d * d * d) - 2 * A[0, 2, i] / (d * d) + A[0, 4, i]
                                - A[0, 0, i] / (d * d * d * d) + A[0, 1, i] / (2 * d * d * d);
                        }
                    }
                    if (i == 1)
                    {
                        Gauss[i, i + 2] = A[0, 0, i] / (d * d * d * d) + A[0, 1, i] / (2 * d * d * d);
                        Gauss[i, i + 1] = A[0, 2, i] / (d * d) + A[0, 3, i] / (2 * d) - 4 * A[0, 0, i] / (d * d * d * d)
                            - A[0, 1, i] / (d * d * d);
                        Gauss[i, i] = 6 * A[0, 0, i] / (d * d * d * d) - 2 * A[0, 2, i] / (d * d) + A[0, 4, i];
                        Gauss[i, i - 1] = A[0, 2, i] / (d * d) - A[0, 3, i] / (2 * d) - 4 * A[0, 0, i] / (d * d * d * d)
                            + A[0, 1, i] / (d * d * d);
                    }
                    if (i == Ra - 2)
                    {
                        Gauss[i, i + 1] = A[0, 2, i] / (d * d) + A[0, 3, i] / (2 * d) - 4 * A[0, 0, i] / (d * d * d * d)
                            - A[0, 1, i] / (d * d * d);
                        Gauss[i, i] = 6 * A[0, 0, i] / (d * d * d * d) - 2 * A[0, 2, i] / (d * d) + A[0, 4, i];
                        Gauss[i, i - 1] = A[0, 2, i] / (d * d) - A[0, 3, i] / (2 * d) - 4 * A[0, 0, i] / (d * d * d * d)
                            + A[0, 1, i] / (d * d * d);
                        Gauss[i, i - 2] = A[0, 0, i] / (d * d * d * d) - A[0, 1, i] / (2 * d * d * d);
                    }
                    if (i == Ra - 1)
                    {
                        if (GrUsl2 == 0)
                        {
                            Gauss[i, i] = 6 * A[0, 0, i] / (d * d * d * d) - 2 * A[0, 2, i] / (d * d) + A[0, 4, i] + A[0, 0, i] / (d * d * d * d) + A[0, 1, i] / (2 * d * d * d);
                            Gauss[i, i - 1] = A[0, 2, i] / (d * d) - A[0, 3, i] / (2 * d) - 4 * A[0, 0, i] / (d * d * d * d)
                                + A[0, 1, i] / (d * d * d);
                            Gauss[i, i - 2] = A[0, 0, i] / (d * d * d * d) - A[0, 1, i] / (2 * d * d * d);
                        }
                        if (GrUsl2 != 0)
                        {
                            Gauss[i, i] = 6 * A[0, 0, i] / (d * d * d * d) - 2 * A[0, 2, i] / (d * d) + A[0, 4, i] - A[0, 0, i] / (d * d * d * d) - A[0, 1, i] / (2 * d * d * d);
                            Gauss[i, i - 1] = A[0, 2, i] / (d * d) - A[0, 3, i] / (2 * d) - 4 * A[0, 0, i] / (d * d * d * d)
                                + A[0, 1, i] / (d * d * d);
                            Gauss[i, i - 2] = A[0, 0, i] / (d * d * d * d) - A[0, 1, i] / (2 * d * d * d);
                        }
                    }
                }
                //Нахождение прогибов методом гаусса-жордана
                MethodGordanGauss(Gauss);

                //копирование прогибов из главной диагонали матрицы

                for (int i = 0; i < Ra; i++)
                    XY[0, i + 3] = (double)(Gauss[i, Ra]);
                //установка граничных значениц для расчёта производных высших порядков
                for (int i = 0; i < CopyXY.Length; i++)
                    CopyXY[i] = XY[0, i];
                BorderZadelSharnkMVIForW(CopyXY, d, Lxy, GrUsl1, GrUsl2, 0);
                for (int i = 0; i < CopyXY.Length; i++)
                    XY[0, i] = CopyXY[i];



                for (int i = 0; i < Ra; i++)
                {
                    XY0d = (XY[0, i + 1 + 3] - XY[0, i - 1 + 3]) / (2 * d);
                    XY0d2 = (XY[0, i + 1 + 3] - 2 * XY[0, i + 3] + XY[0, i - 1 + 3]) / (d * d);
                    XY0d3 = (XY[0, i + 2 + 3] - 2 * XY[0, i + 1 + 3] + 2 * XY[0, i - 1 + 3] - XY[0, i - 2 + 3]) / (2 * d * d * d);
                    XY2d = (XY[2, i + 1 + 3] - XY[2, i - 1 + 3]) / (2 * d);

                    if (i > 0 && i < Ra - 1)
                    {
                        Gauss[i, i + 1] = A[1, 4, i] / (d * d) + A[1, 5, i] / (2 * d);
                        Gauss[i, i] = A[1, 6, i] - 2 * A[1, 4, i] / (d * d);
                        Gauss[i, i - 1] = A[1, 4, i] / (d * d) - A[1, 5, i] / (2 * d);

                    }
                    Gauss[i, Ra] = A[1, 0, i] * XY0d3 + A[1, 1, i] * XY0d2 + A[1, 2, i] * XY0d + A[1, 3, i] * XY[0, i + 3]
                        - A[1, 7, i] * XY2d - A[1, 8, i] * XY[2, i + 3] + A[1, 9, i];
                    //Gauss[i, Ra] =  A[1, 9, i];
                    if (i == 0)
                    {
                        if (GrUsl1 == 0)
                        {
                            Gauss[i, i + 1] = A[1, 4, i] / (d * d) + A[1, 5, i] / (2 * d);
                            Gauss[i, i] = A[1, 6, i] - 2 * A[1, 4, i] / (d * d);
                        }
                        if (GrUsl1 != 0)
                        {
                            Gauss[i, i + 1] = A[1, 4, i] / (d * d) + A[1, 5, i] / (2 * d);
                            Gauss[i, i] = A[1, 6, i] - 2 * A[1, 4, i] / (d * d);
                            Gauss[i, Ra] = Gauss[i, Ra] - (A[1, 4, i] / (d * d) - A[1, 5, i] / (2 * d)) * XY[1, i + 2];
                        }
                    }
                    if (i == Ra - 1)
                    {
                        if (GrUsl2 == 0)
                        {
                            Gauss[i, i] = A[1, 6, i] - 2 * A[1, 4, i] / (d * d);
                            Gauss[i, i - 1] = A[1, 4, i] / (d * d) - A[1, 5, i] / (2 * d);

                        }
                        if (GrUsl2 != 0)
                        {
                            Gauss[i, i] = A[1, 6, i] - 2 * A[1, 4, i] / (d * d);
                            Gauss[i, i - 1] = A[1, 4, i] / (d * d) - A[1, 5, i] / (2 * d);
                            Gauss[i, Ra] = Gauss[i, Ra] - (A[1, 4, i] / (d * d) + A[1, 5, i] / (2 * d)) * XY[1, i + 2];
                        }
                    }
                }
                //Нахождение гамма x методом гаусса-жордана
                MethodGordanGauss(Gauss);


                for (int i = 0; i < Ra; i++)
                    XY[1, i + 3] = (double)(Gauss[i, Ra]);

                for (int i = 0; i < CopyXY.Length; i++)
                    CopyXY[i] = XY[1, i];
                BorderZadelSharnkMVIForW(CopyXY, d, Lxy, GrUsl1, GrUsl2, 1);
                for (int i = 0; i < CopyXY.Length; i++)
                    XY[1, i] = CopyXY[i];


                for (int i = 0; i < Ra; i++)
                {
                    XY0d = (XY[0, i + 1 + 3] - XY[0, i - 1 + 3]) / (2 * d);
                    XY0d2 = (XY[0, i + 1 + 3] - 2 * XY[0, i + 3] + XY[0, i - 1 + 3]) / (d * d);
                    XY0d3 = (XY[0, i + 2 + 3] - 2 * XY[0, i + 1 + 3] + 2 * XY[0, i - 1 + 3] - XY[0, i - 2 + 3]) / (2 * d * d * d);
                    XY1d = (XY[1, i + 1 + 3] - XY[1, i - 1 + 3]) / (2 * d);

                    if (i > 0 && i < Ra - 1)
                    {
                        Gauss[i, i + 1] = A[2, 4, i] / (d * d) + A[2, 5, i] / (2 * d);
                        Gauss[i, i] = A[2, 6, i] - 2 * A[2, 4, i] / (d * d);
                        Gauss[i, i - 1] = A[2, 4, i] / (d * d) - A[2, 5, i] / (2 * d);

                    }
                    Gauss[i, Ra] = A[2, 0, i] * XY0d3 + A[2, 1, i] * XY0d2 + A[2, 2, i] * XY0d + A[2, 3, i] * XY[0, i + 3]
                        - A[2, 7, i] * XY1d - A[2, 8, i] * XY[1, i + 3] + A[2, 9, i];
                    if (i == 0)
                    {
                        if (GrUsl1 == 0)
                        {
                            Gauss[i, i + 1] = A[2, 4, i] / (d * d) + A[2, 5, i] / (2 * d);
                            Gauss[i, i] = A[2, 6, i] - 2 * A[2, 4, i] / (d * d);
                        }
                        if (GrUsl1 != 0)
                        {
                            Gauss[i, i + 1] = A[2, 4, i] / (d * d) + A[2, 5, i] / (2 * d);
                            Gauss[i, i] = A[2, 6, i] - 2 * A[2, 4, i] / (d * d);
                            Gauss[i, Ra] = Gauss[i, Ra] - (A[2, 4, i] / (d * d) - A[2, 5, i] / (2 * d)) * XY[2, i + 2];
                        }
                    }
                    if (i == Ra - 1)
                    {
                        if (GrUsl2 == 0)
                        {
                            Gauss[i, i] = A[2, 6, i] - 2 * A[2, 4, i] / (d * d);
                            Gauss[i, i - 1] = A[2, 4, i] / (d * d) - A[2, 5, i] / (2 * d);

                        }
                        if (GrUsl2 != 0)
                        {
                            Gauss[i, i] = A[2, 6, i] - 2 * A[2, 4, i] / (d * d);
                            Gauss[i, i - 1] = A[2, 4, i] / (d * d) - A[2, 5, i] / (2 * d);
                            Gauss[i, Ra] = Gauss[i, Ra] - (A[2, 4, i] / (d * d) + A[2, 5, i] / (2 * d)) * XY[2, i + 2];
                        }
                    }
                }
                //Нахождение гамма y методом гаусса-жордана
                MethodGordanGauss(Gauss);


                for (int i = 0; i < Ra; i++)
                    XY[2, i + 3] = (double)(Gauss[i, Ra]);

                for (int i = 0; i < CopyXY.Length; i++)
                    CopyXY[i] = XY[2, i];
                BorderZadelSharnkMVIForW(CopyXY, d, Lxy, GrUsl1, GrUsl2, 1);
                for (int i = 0; i < CopyXY.Length; i++)
                    XY[2, i] = CopyXY[i];

                //for (int i = 0; i < FXnYn.GetLength(0); i++)
                //    for (int j = 0; j < FXnYn.GetLength(1); j++)
                //    {
                //        if (MaxDifferents < Math.Abs(XY[i, j] - FXnYn[i, j])) MaxDifferents = Math.Abs(XY[i, j] - FXnYn[i, j]);
                //        FXnYn[i, j] = XY[i, j];
                //    }
                //ShowMassiv2(FXnYn);
                //if (MaxDifferents < eps) break;
                break;
            }

            return XY;
        }

        //Метод изменения граничных условий при счёте МВИ
        private void BorderZadelSharnkMVIForW(double[] AB, double d, double l, int GrUsl1, int GrUsl2, int type)
        {
            //Граничные условия жёсткой заделки
            //Слева 
            if (type == 0)
            {
                if (GrUsl1 == 0)
                {
                    AB[2] = 0;
                    AB[1] = AB[3];
                }
                else
                {
                    AB[2] = 0;
                    AB[1] = -AB[3];
                }
                //Справа
                if (GrUsl2 == 0)
                {
                    AB[AB.Length - 3] = 0;
                    AB[AB.Length - 2] = AB[AB.Length - 4];
                }
                else
                {
                    AB[AB.Length - 3] = 0;
                    AB[AB.Length - 2] = -AB[AB.Length - 4];
                }
            }
            else
            {
                if (GrUsl1 == 0)
                {
                    AB[2] = 0;
                    AB[1] = -AB[3];
                }
                else
                {
                    AB[1] = AB[3];
                }
                //Справа
                if (GrUsl2 == 0)
                {
                    AB[AB.Length - 3] = 0;
                    AB[AB.Length - 2] = -AB[AB.Length - 4];
                }
                else
                {
                    AB[AB.Length - 2] = AB[AB.Length - 4];
                }
            }
            // Далее забиваются крайние значения на основе кубической интерполяции для нахождения
            //для типа 1 необходима интерполяция граничной точки и экстрополяция за граничной точки 
            // для типа 1 учитывается степень свободы.

            // производной 4 порядка на границе
            double[,] GausA = new double[4, 5];


            if (type != 0)
            {
                if (GrUsl1 != 0)
                {
                    for (int i = 0; i < 4; i++)
                        for (int j = 0; j < 5; j++)
                            GausA[i, j] = 0;
                    // первый слой 
                    GausA[0, 0] = 1;
                    GausA[0, 1] = (double)Math.Pow(-d, 1.0);
                    GausA[0, 2] = (double)Math.Pow(-d, 2.0);
                    GausA[0, 3] = (double)Math.Pow(-d, 3.0);
                    GausA[0, 4] = AB[1];

                    for (int i = 1; i < 4; i++)
                    {
                        //Заполнение массива А
                        GausA[i, 0] = 1;
                        GausA[i, 1] = (double)Math.Pow(i * d, 1.0);
                        GausA[i, 2] = (double)Math.Pow(i * d, 2.0);
                        GausA[i, 3] = (double)Math.Pow(i * d, 3.0);
                        GausA[i, 4] = AB[i + 2];
                    }
                    MethodGordanGauss(GausA);

                    //double[,] A = new double[3, 4];
                    AB[2] = (double)(GausA[0, 4]);
                    AB[0] = (double)(GausA[3, 4] * (double)Math.Pow(-2 * d, 3.0) + GausA[2, 4] * (double)Math.Pow(-2 * d, 2.0) + GausA[1, 4] * (-2 * d) + GausA[0, 4]);
                }
                else
                {
                    for (int i = 0; i < 4; i++)
                        for (int j = 0; j < 5; j++)
                            GausA[i, j] = 0;
                    // первый слой 
                    GausA[0, 0] = 1;
                    // второй слой
                    GausA[1, 0] = 1;
                    GausA[1, 1] = (double)Math.Pow(-d, 1.0);
                    GausA[1, 2] = (double)Math.Pow(-d, 2.0);
                    GausA[1, 3] = (double)Math.Pow(-d, 3.0);
                    GausA[1, 4] = AB[1];

                    for (int i = 2; i < 4; i++)
                    {
                        //Заполнение массива А
                        GausA[i, 0] = 1;
                        GausA[i, 1] = (double)Math.Pow((i - 1) * d, 1.0);
                        GausA[i, 2] = (double)Math.Pow((i - 1) * d, 2.0);
                        GausA[i, 3] = (double)Math.Pow((i - 1) * d, 3.0);
                        GausA[i, 4] = AB[i + 1];
                    }
                    MethodGordanGauss(GausA);
                    AB[0] = (double)(GausA[3, 4] * (double)Math.Pow(-2 * d, 3.0) + GausA[2, 4] * (double)Math.Pow(-2 * d, 2.0) + GausA[1, 4] * (-2 * d) + GausA[0, 4]);
                }
                if (GrUsl2 != 0)
                {
                    for (int i = 0; i < 3; i++)
                    {
                        //Заполнение массива А
                        GausA[i, 0] = 1;
                        GausA[i, 1] = (double)Math.Pow(l + (i - 4) * d, 1.0);
                        GausA[i, 2] = (double)Math.Pow(l + (i - 4) * d, 2.0);
                        GausA[i, 3] = (double)Math.Pow(l + (i - 4) * d, 3.0);
                        GausA[i, 4] = AB[AB.Length - 7 + i];
                    }
                    GausA[3, 0] = 1;
                    GausA[3, 1] = (double)Math.Pow(l + d, 1.0);
                    GausA[3, 2] = (double)Math.Pow(l + d, 2.0);
                    GausA[3, 3] = (double)Math.Pow(l + d, 3.0);
                    GausA[3, 4] = AB[AB.Length - 2];
                    MethodGordanGauss(GausA);
                    AB[AB.Length - 1] = (double)(GausA[3, 4] * (double)Math.Pow(l + 2 * d, 3.0) + GausA[2, 4]
                        * (double)Math.Pow(l + 2 * d, 2.0) + GausA[1, 4] * (l + 2 * d) + GausA[0, 4]);
                    AB[AB.Length - 3] = (double)(GausA[3, 4] * (double)Math.Pow(l, 3.0) + GausA[2, 4]
                        * (double)Math.Pow(l, 2.0) + GausA[1, 4] * (l) + GausA[0, 4]);
                }
                else
                {
                    for (int i = 0; i < 4; i++)
                    {
                        //Заполнение массива А
                        GausA[i, 0] = 1;
                        GausA[i, 1] = (double)Math.Pow(l + (i - 2) * d, 1.0);
                        GausA[i, 2] = (double)Math.Pow(l + (i - 2) * d, 2.0);
                        GausA[i, 3] = (double)Math.Pow(l + (i - 2) * d, 3.0);
                        GausA[i, 4] = AB[AB.Length - 5 + i];
                    }
                    MethodGordanGauss(GausA);
                    AB[AB.Length - 1] = (double)(GausA[3, 4] * (double)Math.Pow(l + 2 * d, 3.0) + GausA[2, 4]
                        * (double)Math.Pow(l + 2 * d, 2.0) + GausA[1, 4] * (l + 2 * d) + GausA[0, 4]);
                }
            }
            else
            {
                for (int i = 0; i < 4; i++)
                    for (int j = 0; j < 5; j++)
                        GausA[i, j] = 0;
                // первый слой 
                GausA[0, 0] = 1;
                // второй слой
                GausA[1, 0] = 1;
                GausA[1, 1] = (double)Math.Pow(-d, 1.0);
                GausA[1, 2] = (double)Math.Pow(-d, 2.0);
                GausA[1, 3] = (double)Math.Pow(-d, 3.0);
                GausA[1, 4] = AB[1];

                for (int i = 2; i < 4; i++)
                {
                    //Заполнение массива А
                    GausA[i, 0] = 1;
                    GausA[i, 1] = (double)Math.Pow((i - 1) * d, 1.0);
                    GausA[i, 2] = (double)Math.Pow((i - 1) * d, 2.0);
                    GausA[i, 3] = (double)Math.Pow((i - 1) * d, 3.0);
                    GausA[i, 4] = AB[i + 1];
                }
                MethodGordanGauss(GausA);
                AB[0] = (double)(GausA[3, 4] * (double)Math.Pow(-2 * d, 3.0) + GausA[2, 4] * (double)Math.Pow(-2 * d, 2.0) + GausA[1, 4] * (-2 * d) + GausA[0, 4]);

                for (int i = 0; i < 4; i++)
                {
                    //Заполнение массива А
                    GausA[i, 0] = 1;
                    GausA[i, 1] = (double)Math.Pow(l + (i - 2) * d, 1.0);
                    GausA[i, 2] = (double)Math.Pow(l + (i - 2) * d, 2.0);
                    GausA[i, 3] = (double)Math.Pow(l + (i - 2) * d, 3.0);
                    GausA[i, 4] = AB[AB.Length - 5 + i];
                }
                MethodGordanGauss(GausA);
                AB[AB.Length - 1] = (double)(GausA[3, 4] * (double)Math.Pow(l + 2 * d, 3.0) + GausA[2, 4]
                    * (double)Math.Pow(l + 2 * d, 2.0) + GausA[1, 4] * (l + 2 * d) + GausA[0, 4]);
            }
        }
        private void MethodVariationIteration(Load F, int TypeBorder1, int TypeBorder2, int TypeBorder3, int TypeBorder4)
        {
            //Основные 2 Функции по каздой переменной
            //кол.т    0   1 .....  N-2  N-1(г)
            //A = 0 1 2(г) 3 ....    N  N+1(г) N+2 N+3
            //p = - - -(г) 0 1 2..  N-3  -(г)   -   -
            //DM = - 0 1(г) 2 ....   N-1  N(г)  N+1  - 
            Xn = new double[3, N + 4];// по x
            Yn = new double[3, M + 4];// по y
                                      // Задаём начальную функцию по y 
            for (int y = 0; y < M + 4; y++)
            {
                Xn[0, y] = (double)(Math.Sin(dy * (y - 2)));
                Xn[1, y] = (double)(Math.Sin(dy * (y - 2)));
                Xn[2, y] = (double)(Math.Sin(dy * (y - 2)));
            }
            for (int x = 0; x < M + 4; x++)
            {
                Yn[1, x] = 0;
                Yn[2, x] = 0;
            }
            CountIterationInMethods = 0;
            MaxWForMethodIteration1 = 0;
            MaxWForMethodIteration2 = 0;
            // Метод вариационных итераций
            // Будем искать функцию относитльно x
            for (int k = 1; k <= BorderExitIterationMethods; k++)
            {
                // 
                MaxWForMethodIteration1 = 0;
                MaxWForMethodIteration2 = 0;

                // задаём размерность массивам по переменной х, 3 - уравнения, 14 - максимальное количество переменных коэффициентов при уравнении из системы
                A = new double[3, 14, M - 2];

                // находим переменные функции стоящие при производных в ДУ.
                Pn(A, Xn, false, F.F);
                //ShowMassiv(p6);
                // находим решение уравнения методом гаусса жордана
                Yn = SolutionSystemsProblemMetodGaussaForVariationIteration(Yn, dx, A, TypeBorder4, TypeBorder2, false, 0.0000001);
                // относительно граничного условия и кубической интерполяцией задаём заграничные значения.
                // задаём размерность массивам по переменной х, 3 - уравнения, 14 - максимальное количество переменных коэффициентов при уравнении из системы
                A = new double[3, 14, N - 2];
                //ShowMassiv2(Yn);

                // находим переменные функции стоящие при производных в ДУ.
                Pn(A, Yn, true, F.F);
                //ShowMassiv(p6);
                // находим решение уравнения методом гаусса жордана
                Xn = SolutionSystemsProblemMetodGaussaForVariationIteration(Xn, dy, A, TypeBorder4, TypeBorder2, true, 0.0000001);
                // относительно граничного условия и кубической интерполяцией задаём заграничные значения.
                //ShowMassiv2(Xn);

                CountIterationInMethods++;
                for (int i = 0; i < N + 2; i++)
                    for (int j = 0; j < M + 2; j++)
                    {
                        if (Math.Abs(W[i, j] - Xn[0, i + 1] * Yn[0, j + 1]) > MaxWForMethodIteration1) MaxWForMethodIteration1 = Math.Abs(W[i, j] - Xn[0, i + 1] * Yn[0, j + 1]);// нахождение максмального отклонения найденного решения.
                        W[i, j] = Xn[0, i + 1] * Yn[0, j + 1];
                        //if (Math.Abs(gamx[i, j] - Xn[1, i + 1] * Yn[1, j + 1]) > MaxWForMethodIteration1) MaxWForMethodIteration1 = Math.Abs(gamx[i, j] - Xn[1, i + 1] * Yn[1, j + 1]);// нахождение максмального отклонения найденного решения.
                        gamx[i, j] = Xn[1, i + 1] * Yn[1, j + 1];
                        //if (Math.Abs(gamy[i, j] - Xn[2, i + 1] * Yn[2, j + 1]) > MaxWForMethodIteration1) MaxWForMethodIteration1 = Math.Abs(gamy[i, j] - Xn[2, i + 1] * Yn[2, j + 1]);// нахождение максмального отклонения найденного решения.
                        gamy[i, j] = Xn[2, i + 1] * Yn[2, j + 1];
                        if (Math.Abs(W[i, j]) > MaxWForMethodIteration2) MaxWForMethodIteration2 = Math.Abs(W[i, j]); // Нахождение максимального прогиба
                    }
                if (Math.Abs(MaxWForMethodIteration1) < 0.0000001) break;
                //ShowMassiv2(W);
                //ShowMassiv2(gamx);
                //ShowMassiv2(gamy);
            }
            if (MaximumCountIterationOfMethodAllDecision < CountIterationInMethods) MaximumCountIterationOfMethodAllDecision = CountIterationInMethods;
            //ShowMassiv2(W);
            //ShowMassiv2(gamx);
            //ShowMassiv2(gamy);
            MaxCountIterationSpendOfMetods += CountIterationInMethods;
            mW = MaxWForMethodIteration2;
        }

        /// ........................................................................
        private void MethodABS(Load F, int TypeBorder1, int TypeBorder2, int TypeBorder3, int TypeBorder4, DelegatForMethods Method)
        {
        }
        /// Метод пересчёта переменных параметров упрогости
        public double CountIterationPhisicalNoneleneary { get; protected set; } //Количесвто итераций затраченных на подсчёт переменных параметров упругости
        public double MaxCountMethodOfIteration { get; protected set; } // максимальное значение итераций затрачиваемых итерационным методом
        private void PhisicalNolenearyPart(Load F, int TypeBorder1, int TypeBorder2, int TypeBorder3, int TypeBorder4, bool FlagABS, DelegatForMethods Method)
        {
            double[,] Wlast = new double[N + 2, M + 2];
            double epsilon = 0;
            for (int i = 0; i < N + 2; i++)
                for (int j = 0; j < M + 2; j++)
                    Wlast[i, j] = 0;
            //CountIterationABS = 0;
            CountIterationPhisicalNoneleneary = 0;
            for (int fizi = 1; fizi <= BorderExitPhisicalIteration; fizi++)
            {

                MaxWForMethodIteration1 = 0;
                //Итерационный метод
                if (FlagABS)
                {
                    MethodABS(F, TypeBorder1, TypeBorder2, TypeBorder3, TypeBorder4, Method);
                    epsilon = 0.00001;
                }
                else
                {
                    Method(F, TypeBorder1, TypeBorder2, TypeBorder3, TypeBorder4);
                    epsilon = 0.000001;
                }
                if (!FlagFisicalNonLineary) break;
                //ShowMassiv2(W);
                //Максимальное колличество итераци затрачиваемых итерационным методом;
                //Счётчик подсчёта переменных параметров упругости
                CountIterationPhisicalNoneleneary++;
                //Console.WriteLine(CountIterationPhisicalNoneleneary);
                for (int i = 0; i < N + 2; i++)
                    for (int j = 0; j < M + 2; j++)
                    {
                        if (Math.Abs(W[i, j] - Wlast[i, j]) > MaxWForMethodIteration1) MaxWForMethodIteration1 = Math.Abs(W[i, j] - Wlast[i, j]);// нахождение максмального отклонения найденного решения.
                        Wlast[i, j] = W[i, j];
                    }

                if (Math.Abs(MaxWForMethodIteration1) < epsilon || Math.Abs(MaxWForMethodIteration1) > 100) break;

                //ShowMassiv2(W);
                //Console.WriteLine(mW+" "+ CountIterationPhisicalNoneleneary);
                // пересчёт всех значений с новыми прогибами
                Load();
                //Console.WriteLine(fizi);

            }
        }
    }
}
