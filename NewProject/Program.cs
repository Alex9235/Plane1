using Microsoft.VisualBasic;
using System;
using System.Globalization;
using System.IO;
using System.Linq;
using System.Linq.Expressions;
using System.Diagnostics;
using System.Threading;

namespace Decoder
{
    class Program
    {
        static void WriteMassivInFile(double[,] Mass, string j,bool invers)
        {
            System.Threading.Thread.CurrentThread.CurrentCulture = new System.Globalization.CultureInfo("en-US");// смена точки на запятую
            FileStream file1 = new FileStream("File_" + j.ToString() + ".txt", FileMode.Create); //создаем файловый поток
            StreamWriter writer = new StreamWriter(file1);//создаем «потоковый писатель» и связываем его с файловым потоком 

            //Вывод массива данных в файл
            if (invers)
                for (int i = 0; i < Mass.GetLength(0); i++)
                {
                    for (int k = 0; k < Mass.GetLength(1); k++)
                        writer.Write(String.Format("{0:0.00000000 }", Mass[i, k]));
                    writer.WriteLine();
                }
            else
                for (int k = 0; k < Mass.GetLength(1); k++)
                {
                    for (int i = 0; i < Mass.GetLength(0); i++)
                        writer.Write(String.Format("{0:0.00000000 }", Mass[i, k]));
                    writer.WriteLine();
                }
            writer.Close(); //закрываем поток. Не закрыв поток, в файл ничего не запишется так как первичная запись идёт в ОЗУ
        }
        static void WriteMassivInFile(double[,,] Mass, int index, string j)
        {
            System.Threading.Thread.CurrentThread.CurrentCulture = new System.Globalization.CultureInfo("en-US");// смена точки на запятую
            FileStream file1 = new FileStream("File_" + j.ToString() + ".txt", FileMode.Create); //создаем файловый поток
            StreamWriter writer = new StreamWriter(file1);//создаем «потоковый писатель» и связываем его с файловым потоком 

            //Вывод массива данных в файл

            for (int i = 0; i < Mass.GetLength(0); i++)
            {
                for (int k = 0; k < Mass.GetLength(1); k++)
                {
                    writer.Write(String.Format("{0:0.00000000 }", Mass[i, k, index]));
                }
                writer.WriteLine();
            }

            writer.Close(); //закрываем поток. Не закрыв поток, в файл ничего не запишется так как первичная запись идёт в ОЗУ
        }
        static void WriteNumberInFile(double N, string j)
        {
            System.Threading.Thread.CurrentThread.CurrentCulture = new System.Globalization.CultureInfo("en-US");// смена точки на запятую
            FileStream file1 = new FileStream("File_" + j.ToString() + ".txt", FileMode.Append); //создаем файловый поток
            StreamWriter writer = new StreamWriter(file1);//создаем «потоковый писатель» и связываем его с файловым потоком 
            writer.WriteLine(String.Format("{0:0.00000000 }", N));
            writer.Close(); //закрываем поток. Не закрыв поток, в файл ничего не запишется так как первичная запись идёт в ОЗУ
        }
        static void WriteNumberInFile(double N, double W, string j)
        {
            System.Threading.Thread.CurrentThread.CurrentCulture = new System.Globalization.CultureInfo("en-US");// смена точки на запятую
            FileStream file1 = new FileStream("File_" + j.ToString() + ".txt", FileMode.Append); //создаем файловый поток
            StreamWriter writer = new StreamWriter(file1);//создаем «потоковый писатель» и связываем его с файловым потоком 
            writer.WriteLine(String.Format("{0:0.00000000} {1:0.00000000} ", N,W));
            writer.Close(); //закрываем поток. Не закрыв поток, в файл ничего не запишется так как первичная запись идёт в ОЗУ
        }
        static void WriteMassivInFile(double[] Mass, string j)
        {
            System.Threading.Thread.CurrentThread.CurrentCulture = new System.Globalization.CultureInfo("en-US");// смена точки на запятую
            FileStream file1 = new FileStream("File_" + j.ToString() + ".txt", FileMode.Create); //создаем файловый поток
            StreamWriter writer = new StreamWriter(file1);//создаем «потоковый писатель» и связываем его с файловым потоком 

            //Вывод массива данных в файл
            for (int i = 0; i < Mass.Length; i++)
            {
                writer.Write(String.Format("{0:0.00000000 }", Mass[i]));
                writer.WriteLine();
            }

            writer.Close(); //закрываем поток. Не закрыв поток, в файл ничего не запишется так как первичная запись идёт в ОЗУ
        }
        static void WriteMasssiveListInFile(List<double>[] MassiveLists)
        {
            System.Threading.Thread.CurrentThread.CurrentCulture = new System.Globalization.CultureInfo("en-US");// смена точки на запятую
            FileStream file1 = new FileStream("File_Fuul" + ".txt", FileMode.OpenOrCreate); //создаем файловый поток
            StreamWriter writer = new StreamWriter(file1);//создаем «потоковый писатель» и связываем его с файловым потоком 
            int MaxlengthLists = 0;
            foreach (var Mass in MassiveLists)
                if (MaxlengthLists < Mass.Count) MaxlengthLists = Mass.Count;
            bool FlagForNegativeW = false;
            for (int j = 0; j < MaxlengthLists; j++)
            {
                //Вывод массива данных в файл
                for (int i = 0; i < MassiveLists.Length; i++)
                {
                    if(j<MassiveLists[i].Count)
                        if(MassiveLists[i][j]<0)
                            writer.Write(String.Format("{0:0.00000000} ", MassiveLists[i][j]));
                        else
                            writer.Write(String.Format(" {0:0.00000000} ", MassiveLists[i][j]));
                    else
                        writer.Write("            ");
 
                }
                writer.WriteLine();
            }
            writer.Close(); //закрываем поток. Не закрыв поток, в файл ничего не запишется так как первичная запись идёт в ОЗУ
        }
        static void RasspredDiform(double[,,] Diff, double esi, string s)
        {
            double[,] Exit = new double[Diff.GetLength(0), Diff.GetLength(1)];
            for (int i = 0; i < Diff.GetLength(0); i++)
            {
                for (int j = 0; j < Diff.GetLength(1); j++)
                {
                    for (int k = 0; k <= (Diff.GetLength(2)-1)/2; k++)
                    {
                        if (Diff[i, j, k] <= esi)
                        {
                            Exit[i, j] = 1.0 - k / ((Diff.GetLength(2) - 1) / 2.0);
                            break;
                        }
                    }
                }
            }
            WriteMassivInFile(Exit, s,true);
        }

        static void Main(string[] args)
        {
            
            // наилучшие 30 30 6
            // начальные данные 
            int N = 30;
            int M = 30;
            int P = 10;//(2*P-1)   
            double n = 1;
            double m = 1;
            double p = 0.05; //0.05   0.039935;
            double nu = 0;// коэффициент Пуассона
            double sigmas = 0; //предел текучести 
            double es = 0; //деформация текучести она получается за счёт предела текучести и модуля сдвига es=sigmas/(3*G0) 
            double G0 = 0;//коэффициент модуля сдвига для сжатия
            double QTemp = 0; //Температурная нагрузка
            double alf = 0; //коэффициент линейного расширения
            double alf1 = 0; //коэффициент линейного расширения для второго материала в гридиентном распределении
            double l = 0;//коэффициент масшабности для нано структур
            // для разномодульной задачи коэффициенты для сжатия
            double nu2 = 0;// коэффициент Пуассона
            double sigmas1 = 0;
            double es1 = 0;
            double G1 = 0;//коэффициент модуля сдвига для растяжения
            // для задачи с линеныйм упрочнением
            double G01 = 0;
            // переменные для граничных условий
            // ориентация границ по часовой стрелки!!!!!
            int TypeBorder1 = 0;
            int TypeBorder2 = 0;
            int TypeBorder3 = 0;
            int TypeBorder4 = 0;
            //массивы для вывода эпюров прогиба
            double[] PrintW2X = new double[N];
            double[] PrintW2Y = new double[M];
            double[] PrintWX = new double[N];
            double[] PrintWY = new double[M];
            //массив для двойного интегралла
            double DoubleInt = 0;
            //Ессли задача связанна со временем
            double t0 = 0;
            double t1 = 0;
            double dissipation = 0;
            //градиент с материалом на ходящийся на нижней грани пластины
            double nu1 = 0.24; //коэффициент Пуассонта материала
            double E1 = 348430; //модуль юнга материала
            double RaspredGradient = 0; // коэффициент определяющий тип распределение градиента.
            double Poristostb = 0; // коэффициент пористости
            int ModelPoristostb = 0; // модель пористости (4) (1) x-P (2)
            int ModelFunctionСurvilinearPlane = 0; // модель криволиненой поверхности пластинки
            double C = 3; //величина влажности 
            double Cgamma = 0;


            // определяющие параметры используемы при счёте
            int model = 26;
            // тип граничных условий
            int typeGrandIF = 1;
            //Тип нагрузки
            int typeQ = 1;
            // Метод 
            int TypeMethod = 1;
            // Слой печатания
            int sloi = 5;
            // Тип темпераутры 
            int typeT = 0;
            //количество коэф бубнова галёркин
            int a = 9;
            // данные по керамике
            //Режим счёта с одной нагрузкой и множеством
            bool flagNagruzki = true;
            int U1 = 0; // начальное знавение нагрузки 
            int U2 = 400; // коненчное значение нагрузки 
            double Q = 20; //шаг нагрузки
            // https://www.lib.tpu.ru/fulltext/v/Bulletin_TPU/2006/v309/i2/06.pdf
            int COUNTbet = 4;
            int countb = 0;
            double[,] MassivePointsFisicalNonlineary = new double[4, COUNTbet];
            List<double>[] MasiiveOfList = new List<double>[COUNTbet];
            
            for (int bet = 0; bet < 4; bet++) 
            {

                MasiiveOfList [countb] = new List<double>();
                switch (bet)
                {
                    case 0:
                        RaspredGradient = 0;
                        Poristostb = 0;
                        ModelPoristostb = 0;
                        typeGrandIF = 1;
                        typeT = 0;
                        break;
                    case 1:
                        RaspredGradient = 0.1;
                        Poristostb = 0.8;
                        ModelPoristostb = 1;
                        typeGrandIF = 1;
                        typeT = 0;
                        break;
                    case 2:
                        RaspredGradient = 0.1;
                        Poristostb = 0.8;
                        ModelPoristostb = 2;
                        typeGrandIF = 1;
                        typeT = 0;

                        break;
                    case 3:
                        RaspredGradient = 0.1;
                        Poristostb = 0.8;
                        ModelPoristostb = 3;
                        typeGrandIF = 1;
                        typeT = 0;

                        break;
                    case 4:
                        RaspredGradient = 0;
                        Poristostb = 0;
                        ModelPoristostb = 0;
                        typeGrandIF = 2;
                        typeT = 0;

                        break;
                    case 5:
                        RaspredGradient = 0.1;
                        Poristostb = 0.8;
                        ModelPoristostb = 1;
                        typeGrandIF = 2;
                        typeT = 0;

                        break;
                    case 6:
                        RaspredGradient = 0.1;
                        Poristostb = 0.8;
                        ModelPoristostb = 2;
                        typeGrandIF = 2;
                        typeT = 0;

                        break;
                    case 7:
                        RaspredGradient = 0.1;
                        Poristostb = 0.8;
                        ModelPoristostb = 3;
                        typeGrandIF = 2;
                        typeT = 0;


                        break;
                    case 8:
                        RaspredGradient = 0;
                        Poristostb = 0;
                        ModelPoristostb = 0;
                        typeGrandIF = 3;
                        typeT = 0;


                        break;
                    case 9:
                        RaspredGradient = 0.1;
                        Poristostb = 0.8;
                        ModelPoristostb = 1;
                        typeGrandIF = 3;
                        typeT = 0;


                        break;
                    case 10:
                        RaspredGradient = 0.1;
                        Poristostb = 0.8;
                        ModelPoristostb = 2;
                        typeGrandIF = 3;
                        typeT = 0;


                        break;
                    case 11:
                        RaspredGradient = 0.1;
                        Poristostb = 0.8;
                        ModelPoristostb = 3;
                        typeGrandIF = 3;
                        typeT = 0;


                        break;
                    case 12:
                        RaspredGradient = 0;
                        Poristostb = 0;
                        ModelPoristostb = 0;
                        typeGrandIF = 1;
                        typeT = 1;


                        break;
                    case 13:
                        RaspredGradient = 0.1;
                        Poristostb = 0.8;
                        ModelPoristostb = 1;
                        typeGrandIF = 1;
                        typeT = 1;


                        break;
                    case 14:
                        RaspredGradient = 0.1;
                        Poristostb = 0.8;
                        ModelPoristostb = 2;
                        typeGrandIF = 1;
                        typeT = 1;


                        break;
                    case 15:
                        RaspredGradient = 0.1;
                        Poristostb = 0.8;
                        ModelPoristostb = 3;
                        typeGrandIF = 1;
                        typeT = 1;


                        break;
                    case 16:
                        RaspredGradient = 0;
                        Poristostb = 0;
                        ModelPoristostb = 0;
                        typeGrandIF = 2;
                        typeT = 1;


                        break;
                    case 17:
                        RaspredGradient = 0.1;
                        Poristostb = 0.8;
                        ModelPoristostb = 1;
                        typeGrandIF = 2;
                        typeT = 1;


                        break;
                    case 18:
                        RaspredGradient = 0.1;
                        Poristostb = 0.8;
                        ModelPoristostb = 2;
                        typeGrandIF = 2;
                        typeT = 1;


                        break;
                    case 19:
                        RaspredGradient = 0.1;
                        Poristostb = 0.8;
                        ModelPoristostb = 3;
                        typeGrandIF = 2;
                        typeT = 1;


                        break;
                    case 20:
                        RaspredGradient = 0;
                        Poristostb = 0;
                        ModelPoristostb = 0;
                        typeGrandIF = 3;
                        typeT = 1;


                        break;
                    case 21:
                        RaspredGradient = 0.1;
                        Poristostb = 0.8;
                        ModelPoristostb = 1;
                        typeGrandIF = 3;
                        typeT = 1;


                        break;
                    case 22:
                        RaspredGradient = 0.1;
                        Poristostb = 0.8;
                        ModelPoristostb = 2;
                        typeGrandIF = 3;
                        typeT = 1;


                        break;
                    case 23:
                        RaspredGradient = 0.1;
                        Poristostb = 0.8;
                        ModelPoristostb = 3;
                        typeGrandIF = 3;
                        typeT = 1;


                        break;
                    case 24:
                        RaspredGradient = 0;
                        Poristostb = 0;
                        ModelPoristostb = 0;
                        typeGrandIF = 1;
                        typeT = 2;


                        break;
                    case 25:
                        RaspredGradient = 0.1;
                        Poristostb = 0.8;
                        ModelPoristostb = 1;
                        typeGrandIF = 1;
                        typeT = 2;


                        break;
                    case 26:
                        RaspredGradient = 0.1;
                        Poristostb = 0.8;
                        ModelPoristostb = 2;
                        typeGrandIF = 1;
                        typeT = 2;


                        break;
                    case 27:
                        RaspredGradient = 0.1;
                        Poristostb = 0.8;
                        ModelPoristostb = 3;
                        typeGrandIF = 1;
                        typeT = 2;


                        break;
                    case 28:
                        RaspredGradient = 0;
                        Poristostb = 0;
                        ModelPoristostb = 0;
                        typeGrandIF = 2;
                        typeT = 2;


                        break;
                    case 29:
                        RaspredGradient = 0.1;
                        Poristostb = 0.8;
                        ModelPoristostb = 1;
                        typeGrandIF = 2;
                        typeT = 2;


                        break;
                    case 30:
                        RaspredGradient = 0.1;
                        Poristostb = 0.8;
                        ModelPoristostb = 2;
                        typeGrandIF = 2;
                        typeT = 2;


                        break;
                    case 31:
                        RaspredGradient = 0.1;
                        Poristostb = 0.8;
                        ModelPoristostb = 3;
                        typeGrandIF = 2;
                        typeT = 2;


                        break;
                    case 32:
                        RaspredGradient = 0;
                        Poristostb = 0;
                        ModelPoristostb = 0;
                        typeGrandIF = 3;
                        typeT = 2;


                        break;
                    case 33:
                        RaspredGradient = 0.1;
                        Poristostb = 0.8;
                        ModelPoristostb = 1;
                        typeGrandIF = 3;
                        typeT = 2;


                        break;
                    case 34:
                        RaspredGradient = 0.1;
                        Poristostb = 0.8;
                        ModelPoristostb = 2;
                        typeGrandIF = 3;
                        typeT = 2;


                        break;
                    case 35:
                        RaspredGradient = 0.1;
                        Poristostb = 0.8;
                        ModelPoristostb = 3;
                        typeGrandIF = 3;
                        typeT = 2;


                        break;
                    case 36:
                        RaspredGradient = 0;
                        Poristostb = 0;
                        ModelPoristostb = 0;
                        typeGrandIF = 1;
                        typeT = 3;


                        break;
                    case 37:
                        RaspredGradient = 0.1;
                        Poristostb = 0.8;
                        ModelPoristostb = 1;
                        typeGrandIF = 1;
                        typeT = 3;


                        break;
                    case 38:
                        RaspredGradient = 0.1;
                        Poristostb = 0.8;
                        ModelPoristostb = 2;
                        typeGrandIF = 1;
                        typeT = 3;


                        break;
                    case 39:
                        RaspredGradient = 0.1;
                        Poristostb = 0.8;
                        ModelPoristostb = 3;
                        typeGrandIF = 1;
                        typeT = 3;


                        break;
                    case 40:
                        RaspredGradient = 0;
                        Poristostb = 0;
                        ModelPoristostb = 0;
                        typeGrandIF = 2;
                        typeT = 3;


                        break;
                    case 41:
                        RaspredGradient = 0.1;
                        Poristostb = 0.8;
                        ModelPoristostb = 1;
                        typeGrandIF = 2;
                        typeT = 3;


                        break;
                    case 42:
                        RaspredGradient = 0.1;
                        Poristostb = 0.8;
                        ModelPoristostb = 2;
                        typeGrandIF = 2;
                        typeT = 3;


                        break;
                    case 43:
                        RaspredGradient = 0.1;
                        Poristostb = 0.8;
                        ModelPoristostb = 3;
                        typeGrandIF = 2;
                        typeT = 3;


                        break;
                    case 44:
                        RaspredGradient = 0;
                        Poristostb = 0;
                        ModelPoristostb = 0;
                        typeGrandIF = 3;
                        typeT = 3;


                        break;
                    case 45:
                        RaspredGradient = 0.1;
                        Poristostb = 0.8;
                        ModelPoristostb = 1;
                        typeGrandIF = 3;
                        typeT = 3;


                        break;
                    case 46:
                        RaspredGradient = 0.1;
                        Poristostb = 0.8;
                        ModelPoristostb = 2;
                        typeGrandIF = 3;
                        typeT = 3;


                        break;
                    case 47:
                        RaspredGradient = 0.1;
                        Poristostb = 0.8;
                        ModelPoristostb = 3;
                        typeGrandIF = 3;
                        typeT = 3;


                        break;
                    case 48:
                        RaspredGradient = 0;
                        Poristostb = 0;
                        ModelPoristostb = 0;
                        typeGrandIF = 1;
                        typeT = 4;


                        break;
                    case 49:
                        RaspredGradient = 0.1;
                        Poristostb = 0.8;
                        ModelPoristostb = 1;
                        typeGrandIF = 1;
                        typeT = 4;


                        break;
                    case 50:
                        RaspredGradient = 0.1;
                        Poristostb = 0.8;
                        ModelPoristostb = 2;
                        typeGrandIF = 1;
                        typeT = 4;


                        break;
                    case 51:
                        RaspredGradient = 0.1;
                        Poristostb = 0.8;
                        ModelPoristostb = 3;
                        typeGrandIF = 1;
                        typeT = 4;


                        break;
                    case 52:
                        RaspredGradient = 0;
                        Poristostb = 0;
                        ModelPoristostb = 0;
                        typeGrandIF = 2;
                        typeT = 4;


                        break;
                    case 53:
                        RaspredGradient = 0.1;
                        Poristostb = 0.8;
                        ModelPoristostb = 1;
                        typeGrandIF = 2;
                        typeT = 4;


                        break;
                    case 54:
                        RaspredGradient = 0.1;
                        Poristostb = 0.8;
                        ModelPoristostb = 2;
                        typeGrandIF = 2;
                        typeT = 4;


                        break;
                    case 55:
                        RaspredGradient = 0.1;
                        Poristostb = 0.8;
                        ModelPoristostb = 3;
                        typeGrandIF = 2;
                        typeT = 4;
                        break;
                    case 56:
                        RaspredGradient = 0;
                        Poristostb = 0;
                        ModelPoristostb = 0;
                        typeGrandIF = 3;
                        typeT = 4;
                        break;
                    case 57:
                        RaspredGradient = 0.1;
                        Poristostb = 0.8;
                        ModelPoristostb = 1;
                        typeGrandIF = 3;
                        typeT = 4;


                        break;
                    case 58:
                        RaspredGradient = 0.1;
                        Poristostb = 0.8;
                        ModelPoristostb = 2;
                        typeGrandIF = 3;
                        typeT = 4;
                        break;
                    case 59:
                        RaspredGradient = 0.1;
                        Poristostb = 0.8;
                        ModelPoristostb = 3;
                        typeGrandIF = 3;
                        typeT = 4;
                        break;

                    case 60:
                        RaspredGradient = 0;
                        Poristostb = 0;
                        ModelPoristostb = 0;
                        typeGrandIF = 3;
                        typeT = 4;
                        l = 0.3; 
                        break;
                    case 61:
                        RaspredGradient = 0.1;
                        Poristostb = 0.8;
                        ModelPoristostb = 1;
                        typeGrandIF = 3;
                        typeT = 4;
                        l = 0.3;


                        break;
                    case 62:
                        RaspredGradient = 0.1;
                        Poristostb = 0.8;
                        ModelPoristostb = 2;
                        typeGrandIF = 3;
                        typeT = 4;
                        l = 0.3;
                        break;
                    case 63:
                        RaspredGradient = 0.1;
                        Poristostb = 0.8;
                        ModelPoristostb = 3;
                        typeGrandIF = 3;
                        typeT = 4;
                        l = 0.3;
                        break;
                    case 64:
                        RaspredGradient = 0;
                        Poristostb = 0;
                        ModelPoristostb = 0;
                        typeGrandIF = 3;
                        typeT = 4;
                        l = 0.5;
                        break;
                    case 65:
                        RaspredGradient = 0.1;
                        Poristostb = 0.8;
                        ModelPoristostb = 1;
                        typeGrandIF = 3;
                        typeT = 4;
                        l = 0.5;


                        break;
                    case 66:
                        RaspredGradient = 0.1;
                        Poristostb = 0.8;
                        ModelPoristostb = 2;
                        typeGrandIF = 3;
                        typeT = 4;
                        l = 0.5;
                        break;
                    case 67:
                        RaspredGradient = 0.1;
                        Poristostb = 0.8;
                        ModelPoristostb = 3;
                        typeGrandIF = 3;
                        typeT = 4;
                        l = 0.5;
                        break;
                    case 68:
                        RaspredGradient = 0;
                        Poristostb = 0;
                        ModelPoristostb = 0;
                        typeGrandIF = 3;
                        typeT = 4;
                        l = 0.7;
                        break;
                    case 69:
                        RaspredGradient = 0.1;
                        Poristostb = 0.8;
                        ModelPoristostb = 1;
                        typeGrandIF = 3;
                        typeT = 4;
                        l = 0.7;
                        break;
                    case 70:
                        RaspredGradient = 0.1;
                        Poristostb = 0.8;
                        ModelPoristostb = 2;
                        typeGrandIF = 3;
                        typeT = 4;
                        l = 0.7;
                        break;
                    case 71:
                        RaspredGradient = 0.1;
                        Poristostb = 0.8;
                        ModelPoristostb = 3;
                        typeGrandIF = 3;
                        typeT = 4;
                        l = 0.7;
                        break;

                }
                // создаю объект класса БАЛКА
                //PlaneSheremeteva_Peleha plast = new PlaneSheremeteva_Peleha(N, M, P, n, m, p);
                PlaneKirgofa plast = new PlaneKirgofa(N, M, P, n, m, p);
                // типы исследуемых моделей
                switch (model)
                {
                    case 1:
                        {
                            //allum
                            nu = 0.314;
                            G0 = 26490;
                            sigmas = 98.0665;
                            //sigmas = 0.003702;
                            //es = 0.001234;
                            //l = 0.5;
                            plast.InicializationPlast(nu, G0);
                            plast.InicializationExpFizicalNonlinProblem(sigmas);
                            //plast.InicializationNano(l);
                            break;
                        }
                    case 2:
                        {
                            //steel T=0
                            // для экспоненциального!!!
                            G0 = 83250;
                            sigmas = 728.522;
                            //es = 0.002917;
                            nu = 0.3;
                            plast.InicializationPlast(nu, G0);
                            plast.InicializationExpFizicalNonlinProblem(sigmas);
                            break;
                        }
                    case 3:
                        {
                            //steel T=300;
                            // для экспоненциального!!!
                            G0 = 74996;//коэффициент модуля сдвига для сжатия
                            sigmas = 630.642;
                            //es = 0.002803;
                            nu = 0.3;// коэффициент Пуассона
                            QTemp = 300; //Температурная нагрузка
                            alf = 0.000017; //коэффициент линейного температурного расширения
                            plast.InicializationPlast(nu, G0);
                            plast.InicializationExpFizicalNonlinProblem(sigmas);
                            plast.InicializationTempature(QTemp, alf);
                            break;
                        }
                    case 4:
                        {
                            //steel T=500;
                            // для экспоненциального!!!
                            G0 = 62755;
                            sigmas = 520.555;
                            //es = 0.002765;
                            nu = 0.3;
                            QTemp = 500; //Температурная нагрузка
                            alf = 0.000018; //коэффициент линейного температурного расширения
                            plast.InicializationPlast(nu, G0);
                            plast.InicializationExpFizicalNonlinProblem(sigmas);
                            plast.InicializationTempature(QTemp, alf);
                            break;
                        }
                    case 5:
                        {
                            //steel T=0
                            // для экспоненциального!!!
                            G0 = 83250;
                            sigmas = 728.522;
                            //es = 0.002917;
                            nu = 0.3;
                            l = 0.5;
                            plast.InicializationPlast(nu, G0);
                            plast.InicializationExpFizicalNonlinProblem(sigmas);
                            plast.InicializationNano(l);

                            break;
                        }
                    case 6:
                        {
                            //steel T=300;
                            // для экспоненциального!!!
                            G0 = 74996;//коэффициент модуля сдвига для сжатия
                            sigmas = 630.642;
                            //es = 0.002803;
                            nu = 0.3;// коэффициент Пуассона
                            QTemp = 300; //Температурная нагрузка
                            alf = 0.000017; //коэффициент линейного температурного расширения
                            l = 0.5;
                            plast.InicializationPlast(nu, G0);
                            plast.InicializationExpFizicalNonlinProblem(sigmas);
                            plast.InicializationNano(l);
                            plast.InicializationTempature(QTemp, alf);
                            break;
                        }
                    case 7:
                        {
                            //steel T=500;
                            // для экспоненциального!!!
                            G0 = 62755;
                            sigmas = 520.555;
                            //es = 0.002765;
                            nu = 0.3;
                            QTemp = 500; //Температурная нагрузка
                            alf = 0.000018; //коэффициент линейного температурного расширения
                            l = 0.5;
                            plast.InicializationPlast(nu, G0);
                            plast.InicializationExpFizicalNonlinProblem(sigmas);
                            plast.InicializationNano(l);
                            plast.InicializationTempature(QTemp, alf);
                            break;
                        }
                    case 8:
                        {
                            //steel T=0
                            // для линейного упрочнения
                            G0 = 75000;
                            G1 = 2211.53;
                            es = 0.002917;
                            nu = 0.3;
                            plast.InicializationPlast(nu, G0);
                            plast.InicializationLinearyHard(G1, es);
                            break;
                        }
                    case 9:
                        {
                            //steel T=300;
                            // для линейного упрочнения!!
                            G0 = 69230.76;
                            G1 = 1442.3;
                            es = 0.002803;
                            nu = 0.3;// коэффициент Пуассона
                            QTemp = 300; //Температурная нагрузка
                            alf = 0.000017; //коэффициент линейного температурного расширения
                            plast.InicializationPlast(nu, G0);
                            plast.InicializationLinearyHard(G1, es);
                            plast.InicializationTempature(QTemp, alf);
                            break;
                        }
                    case 10:
                        {
                            //steel T=500;
                            // для линейного упрочнения!!
                            G0 = 60000;
                            G1 = 576.92;
                            es = 0.002765;
                            nu = 0.3;
                            QTemp = 500; //Температурная нагрузка
                            alf = 0.000018; //коэффициент линейного температурного расширения
                            plast.InicializationPlast(nu, G0);
                            plast.InicializationLinearyHard(G1, es);
                            plast.InicializationTempature(QTemp, alf);
                            break;
                        }
                    case 11:
                        {
                            //steel T=0
                            // для линейного упрочнения
                            G0 = 75000;
                            G1 = 2211.53;
                            es = 0.002917;
                            nu = 0.3;
                            l = 0.5;
                            plast.InicializationPlast(nu, G0);
                            plast.InicializationLinearyHard(G1, es);
                            plast.InicializationNano(l);
                            break;
                        }
                    case 12:
                        {
                            //steel T=300;
                            // для линейного упрочнения!!
                            G0 = 69230.76;
                            G1 = 1442.3;
                            es = 0.002803;
                            nu = 0.3;// коэффициент Пуассона
                            QTemp = 300; //Температурная нагрузка
                            alf = 0.000017; //коэффициент линейного температурного расширения
                            l = 0.5;
                            plast.InicializationPlast(nu, G0);
                            plast.InicializationLinearyHard(G1, es);
                            plast.InicializationNano(l);
                            plast.InicializationTempature(QTemp, alf);
                            break;
                        }
                    case 13:
                        { 
                            //steel T=500;
                            // для линейного упрочнения!!
                            G0 = 60000;
                            G1 = 576.92;
                            es = 0.002765;
                            nu = 0.3;
                            QTemp = 500; //Температурная нагрузка
                            alf = 0.000018; //коэффициент линейного температурного расширения
                            l = 0.5;
                            plast.InicializationPlast(nu, G0);
                            plast.InicializationLinearyHard(G1, es);
                            plast.InicializationNano(l);
                            plast.InicializationTempature(QTemp, alf);
                            break;
                        }
                    case 14:
                        {
                            //allum АМг61
                            nu = 0.314;// коэффициент Пуассона
                            sigmas = 193;
                            es = 0.002382;
                            G0 = 27000;//коэффициент модуля сдвига для сжатия
                                       //sigmas1 = 204;
                                       //es1 = 0.002518;
                                       //G1 = 27000;
                            sigmas1 = 0;
                            es1 = 0;
                            G1 = 0;
                            l = 0;
                            RaspredGradient = 0;
                            Poristostb = 0;
                            ModelPoristostb = 0;

                            break;
                        }
                    case 15:
                        {
                            //allum АМг61
                            nu = 0.314;// коэффициент Пуассона
                            sigmas = 193;
                            es = 0.002382;
                            G0 = 27000;//коэффициент модуля сдвига для сжатия
                            nu2 = 0.314;
                            sigmas1 = 204;
                            es1 = 0.008;//0.002518
                            G1 = 8500;
                            l = 0;
                            RaspredGradient = 0;
                            Poristostb = 0;
                            ModelPoristostb = 0;

                            break;
                        }
                    case 16:
                        {
                            //зернистый композит
                            //G1 и G0 пологается брать одинаковым так как оберазмеривание идёт по наименьшему модулю сдвига
                            nu = 0.398;// коэффициент Пуассона
                            sigmas = 100;
                            es = 0.0021235;
                            G0 = 9668;//коэффициент модуля сдвига для сжатия
                            nu2 = 0.276;// коэффициент Пуассона
                            sigmas1 = 100;
                            es1 = 0.0034478;
                            G1 = 9668;//коэффициент модуля сдвига для сжатия
                            l = 0.7;
                            RaspredGradient = 0;
                            Poristostb = 0.4;
                            ModelPoristostb = 2;
                            break;
                        }
                    case 17:
                        {
                            //зернистый композит
                            //G1 и G0 пологается брать одинаковым так как оберазмеривание идёт по наименьшему модулю сдвига
                            nu = 0.276;// коэффициент Пуассона
                            sigmas = 100;
                            es = 0.0034478;
                            G0 = 9668;//коэффициент модуля сдвига для сжатия
                            nu2 = 0.398;// коэффициент Пуассона
                            sigmas1 = 100;
                            es1 = 0.0021235;
                            G1 = 9668;//коэффициент модуля сдвига для сжатия
                                      //G1 = 15697;
                            l = 0.7;
                            RaspredGradient = 0;
                            Poristostb = 0.4;
                            ModelPoristostb = 2;
                            break;
                        }
                    case 18:
                        {
                            //гипотетический алюминий
                            nu = 0.314;// коэффициент Пуассона
                            sigmas = 98.0665;
                            es = 0.001234;
                            G0 = 26490; //коэффициент модуля сдвига для сжатия
                                        //гипотетческое сигма 15 % от существенного.
                            nu2 = 0.314;// коэффициент Пуассона
                            sigmas1 = 112.776475;
                            es1 = 0.001234;
                            G1 = 26490;
                            l = 0.7;
                            RaspredGradient = 0;
                            Poristostb = 0.4;
                            ModelPoristostb = 2;
                            break;
                        }
                    case 19:
                        {
                            //гипотетический алюмминий с усилением на сжатие
                            nu = 0.314;// коэффициент Пуассона
                            sigmas = 98.0665;
                            es = 0.001234;
                            G0 = 26490;//коэффициент модуля сдвига для сжатия
                                       //гипотетческое сигма 15 % от существенного.
                            nu2 = 0.314;// коэффициент Пуассона
                            sigmas1 = 112.776475;
                            es1 = 0.008;
                            G1 = 26490;//G1 = G2 для соответсвия размерности
                            l = 0.7;
                            RaspredGradient = 0;
                            Poristostb = 0.4;
                            ModelPoristostb = 3;
                            break;
                        }
                    case 20:
                        {
                            // для экспоненциального!!!
                            G0 = 103707.4647;//коэффициент модуля сдвига для сжатия
                            nu = 0.3;// коэффициент Пуассона
                            plast.InicializationPlast(nu, G0);
                            sigmas = 710.605;
                            //es = 0.00263539;
                            plast.InicializationExpFizicalNonlinProblem(sigmas);
                            //пористость
                            Poristostb = 0.4;
                            ModelPoristostb = 1;
                            plast.InicializationPorisoty(Poristostb, ModelPoristostb);
                            alf = 0.000018; //коэффициент линейного температурного расширения
                            plast.InicializationTempature("C:/Users/proto/Documents/GitHub/Plane1/Rtemp_2100.txt", alf);
                            break;
                        }
                    case 21:
                        {
                            // для экспоненциального!!!
                            G0 = 103707.4647;//коэффициент модуля сдвига для сжатия
                            nu = 0.3;// коэффициент Пуассона
                            plast.InicializationPlast(nu, G0);
                            sigmas = 710.605;
                            //es = 0.00263539;
                            plast.InicializationExpFizicalNonlinProblem(sigmas);
                            //пористость
                            //Poristostb = 0.4;
                            //ModelPoristostb = 1;
                            //plast.InicializationPorisoty(Poristostb, ModelPoristostb);
                            //коэффициент линейного температурного расширения
                            //l = 0.5;
                            //plast.InicializationNano(l);
                            break;
                        }
                    case 22:
                        {
                            // для экспоненциального!!!
                            G0 = 103707.4647;//коэффициент модуля сдвига для сжатия
                            nu = 0.3;// коэффициент Пуассона
                            plast.InicializationPlast(nu, G0);
                            sigmas = 710.605;
                            //es = 0.00263539;
                            plast.InicializationExpFizicalNonlinProblem(sigmas);
                            //пористость
                            //Poristostb = 0.4;
                            //ModelPoristostb = 1;
                            //plast.InicializationPorisoty(Poristostb, ModelPoristostb);
                            alf = 0.000018; //коэффициент линейного температурного расширения
                            plast.InicializationTempatureWithDefform("C:/Users/proto/Documents/GitHub/Plane1/Rtemp_2100.txt",
                                "C:/Users/proto/Documents/GitHub/Plane1/Распр s(e,t).txt", alf, 273.15, 873.15, 0, 0.012, bet);
                            //l = 0;
                            //plast.InicializationNano(l);
                            break;
                        }
                    case 23:
                        {
                            // для экспоненциального!!!
                            G0 = 103707.4647;//коэффициент модуля сдвига для сжатия
                            nu = 0.3;// коэффициент Пуассона
                            plast.InicializationPlast(nu, G0);
                            sigmas = 710.605;
                            //es = 0.00263539;
                            plast.InicializationExpFizicalNonlinProblem(sigmas);
                            //пористость
                            //Poristostb = 0.4;
                            //ModelPoristostb = 1;
                            //plast.InicializationPorisoty(Poristostb, ModelPoristostb);
                            alf = 0.000018; //коэффициент линейного температурного расширения
                            plast.InicializationTempatureWithDefform("C:/Users/proto/Documents/GitHub/Plane1/локальная_400_15х15.txt",
                                "C:/Users/proto/Documents/GitHub/Plane1/Распр s(e,t).txt", alf, 273.15, 873.15, 0, 0.012, bet);
                            //l = 0;
                            //plast.InicializationNano(l);
                            break;
                        }
                    case 24:
                        {
                            // для экспоненциального!!!
                            G0 = 103707.4647;//коэффициент модуля сдвига для сжатия
                            nu = 0.3;// коэффициент Пуассона
                            plast.InicializationPlast(nu, G0);
                            sigmas = 710.605;
                            //es = 0.00263539;
                            plast.InicializationExpFizicalNonlinProblem(sigmas);
                            //пористость
                            Poristostb = 0.4;
                            ModelPoristostb = 1;
                            plast.InicializationPorisoty(Poristostb, ModelPoristostb);
                            alf = 0.000018; //коэффициент линейного температурного расширения
                            plast.InicializationTempatureWithDefform("C:/Users/proto/Documents/GitHub/Plane1/temp_500.txt",
                            "C:/Users/proto/Documents/GitHub/Plane1/Распр s(e,t).txt", alf, 273.15, 873.15, 0, 0.012, bet);
                            //l = 0;
                            //plast.InicializationNano(l);
                            break;
                        }
                    case 25:
                        {
                            // для экспоненциального!!!
                            G0 = 103707.4647;//коэффициент модуля сдвига для сжатия
                            nu = 0.3;// коэффициент Пуассона
                            plast.InicializationPlast(nu, G0);
                            sigmas = 710.605;
                            //es = 0.00263539;
                            plast.InicializationExpFizicalNonlinProblem(sigmas);
                            //пористость
                            Poristostb = 0.4;
                            ModelPoristostb = 1;
                            plast.InicializationPorisoty(Poristostb, ModelPoristostb);
                            alf = 0.000018; //коэффициент линейного температурного расширения
                            plast.InicializationTempatureWithDefform("C:/Users/proto/Documents/GitHub/Plane1/temp_600.txt",
                                "C:/Users/proto/Documents/GitHub/Plane1/Распр s(e,t).txt", alf, 273.15, 873.15, 0, 0.012, bet); //30 30 7
                            //l = 0;
                            //plast.InicializationNano(l);
                            break;
                        }
                    case 26:
                        {
                            //!!!!!!СНачало добавлется температура потом градиент потом пористость, никак иначе
                            //sus 304
                            nu = 0.3262;
                            G0 = 75795.5;
                            //G1 = 2026;
                            alf = 0.00001233;
                            //sigmas = 98.0665;
                            //sigmas = 0.003702;
                            sigmas = 220;
                            //es = 0.00089822;
                            l = 0.5;
                            plast.InicializationPlast(nu, G0);
                            plast.InicializationExpFizicalNonlinProblem(sigmas);
                            //plast.InicializationMultimodulProblem(nu, G0);
                            //plast.InicializationMultimoduWithExpFizicalNonlinlProblem(sigmas + 100);

                            if (C > 0)
                            {
                                //plast.InicializationHumidity(C, 0, Cgamma, 0.0005);
                                plast.InicializationHumidity("C:/Users/proto/Documents/GitHub/Plane1/локальная_1_21x21_30x30x19_0-1_(20).txt", 0.0005,C);
                            }
                            switch (typeT)
                            {
                                case 1:
                                    plast.InicializationTempature(100, alf);
                                    break;
                                case 2:
                                    plast.InicializationTempature(200, alf);
                                    break;
                                case 3:
                                    //plast.InicializationTempature(0, 100, alf, 0.4); //z = -1 T = 300, z = 0 T = 800,
                                    plast.InicializationTempature("C:/Users/proto/Documents/GitHub/Plane1/локальная_400_21x21_30x30x19_300-400_(20).txt", alf);
                                    for (int x = 0; x < plast.T.GetLength(0); x++)
                                        for (int y = 0; y < plast.T.GetLength(1); y++)
                                            for (int z = 0; z < plast.T.GetLength(2); z++)
                                                plast.T[x, y, z] = plast.T[x, y, z] - 200;
                                    break;
                                case 4:
                                    //plast.InicializationTempature(100, 0, alf, 0.4); //z = -1 T = 800, z = 0 T = 300,
                                    plast.InicializationTempatureInversion("C:/Users/proto/Documents/GitHub/Plane1/локальная_400_21x21_30x30x19_300-400_(20).txt", alf);
                                    for (int x = 0; x < plast.T.GetLength(0); x++)
                                        for (int y = 0; y < plast.T.GetLength(1); y++)
                                            for (int z = 0; z < plast.T.GetLength(2); z++)
                                                plast.T[x, y, z] = plast.T[x, y, z] - 200;
                                    break;
                                case 5:
                                    plast.InicializationTempature("C:/Users/proto/Documents/GitHub/Plane1/локальная_400_30_30_19.txt", alf); //30 30 7
                                    break;
                                case 6:
                                    plast.InicializationTempature("C:/Users/proto/Documents/GitHub/Plane1/локальная_2_400_21x21_30x30x19_300-400_(20).txt", alf); 
                                    break;
                                default: break;
                            }
                            if (typeT != 0)
                                plast.InicializationTempatureWithENUG0(0, 0.3262, -0.00020002, 0.0000003797, 0, 0, 201040, 0.0003079, -0.0000006534, 0, 0, 0.000012330, 0.0008086, 0, 0);
                            
                            if (RaspredGradient > 0 && typeT != 0)
                                plast.InicializationVariationGradient(RaspredGradient, 0, 0.24, 0, 0, 0, 0, 348430, -0.0003070, 0.0000002160, 0.00000000008946, 0, 0.0000058723, 0.0006095, 0, 0);
                            else if(RaspredGradient > 0)
                                plast.InicializationVariationGradient(RaspredGradient, 0.24, 348430);
                            if (Poristostb > 0 && ModelPoristostb != 0)
                                plast.InicializationPorisoty(Poristostb, ModelPoristostb);
                            if (typeT != 0)
                                plast.BezrazmernostbT();
                            if(l!=0)
                                plast.InicializationNano(l);
                            break;
                        }
                    default: break;
                }
                if (plast.Error) return;


                //для нахождения предельной точки пластического шарнира
                bool FlagPlastSharn = true;
                bool FlagFirstNonlinear = true;
                if (flagNagruzki)
                {
                    U1 = 1; // начальное знавение нагрузки 
                    U2 = 1; // коненчное значение нагрузки 
                    FlagPlastSharn = false;
                    FlagFirstNonlinear = false;
                }
                switch (typeGrandIF)
                {
                    case 1:
                        TypeBorder1 = 0;
                        TypeBorder2 = 0;
                        TypeBorder3 = 0;
                        TypeBorder4 = 0;
                        break;
                    case 2:
                        TypeBorder1 = 0;
                        TypeBorder2 = 1;
                        TypeBorder3 = 0;
                        TypeBorder4 = 1;
                        break;
                    case 3:
                        TypeBorder1 = 1;
                        TypeBorder2 = 1;
                        TypeBorder3 = 1;
                        TypeBorder4 = 1;
                        break;
                    default:

                        break;
                }


                Load F = new Load(N, M, typeQ, Q); // массив для нагрузки
                                                   //создаем объект
                Stopwatch stopwatch = new Stopwatch();
                //засекаем время начала операции
                stopwatch.Start();
                //выполняем какую-либо операцию
                //for (int a = 5; a <= 3; a++)
                //{
                //WriteMassivInFile(plast.T, 5, "Temp " + 5);
                //WriteMassivInFile(plast.T, 3, "Temp " + 3);

                int flag1 = 0;
                string Form = String.Format("TQ{0}_TIF{1}_bg{2}_Met{3}_n{4}m{5}p{6}_Q{7}_N{8}M{9}P{10}_t0{11}t1{12}dis{13}_L{14}_model{15}_RG{16}_Modpor{17}_MFP{18}_C{19}_Cgam{20}_TypeT{21}",
                            typeQ, typeGrandIF, a, TypeMethod, n, m, p, Q, N, M, 2 * P - 1, t0, t1, dissipation, l, model, RaspredGradient, ModelPoristostb, ModelFunctionСurvilinearPlane,C,Cgamma,typeT);
                //}
                if (true)
                {

                    for (int i = U1; i <= U2; i++)
                    {
                        F.LoadValue = Q * i;
                        F.FullLoad();
                        ///0-Жёсткая заделка
                        ///1-Шарнир 

                        switch (TypeMethod)
                        {
                            case 1:
                                plast.StaticDecisionMethodVariationIteration(F, 10000, 0, 10000, TypeBorder1, TypeBorder2, TypeBorder3, TypeBorder4);
                                break;
                            case 2:
                                plast.StaticDecisionMethodVariationIteration(F, 1000, 1000, 1000, TypeBorder1, TypeBorder2, TypeBorder3, TypeBorder4);
                                break;
                            case 3:
                                plast.StaticDecisionMethodBubnovGalercin(F, a, 0, 1000, TypeBorder1, TypeBorder2, TypeBorder3, TypeBorder4);
                                break;
                            case 4:
                                plast.StaticDecisionMethodBubnovGalercin(F, a, 1000, 1000, TypeBorder1, TypeBorder2, TypeBorder3, TypeBorder4);
                                break;
                            case 5:
                                plast.DinamiDecisionMethodRugeKyta4(F, t0, t1, 0.0004, 1, TypeBorder1, TypeBorder2, TypeBorder3, TypeBorder4, 200);
                                break;
                            case 6:
                                plast.DinamiDecisionMethodVariationIteration(F, t0, t1, 0.0004, 1, TypeBorder1, TypeBorder2, TypeBorder3, TypeBorder4, 1);
                                break;
                            case 7:
                                plast.StaticDecisionMethodVariationIterartionTwoLevel(F, 1000, 0, 1000, TypeBorder1, TypeBorder2, TypeBorder3, TypeBorder4);
                                break;
                            case 8:
                                plast.StaticDecisionMethodVariationIterartionTwoLevel(F, 1000, 1000, 1000, TypeBorder1, TypeBorder2, TypeBorder3, TypeBorder4);
                                break;
                            case 9:
                                plast.StaticDecisionMethodKantorovichVlasovTwoLevel(F, 0, 0, 1000, TypeBorder1, TypeBorder2, TypeBorder3, TypeBorder4);
                                break;
                            case 10:
                                plast.StaticDecisionMethodKantorovichVlasovTwoLevel(F, 0, 1000, 1000, TypeBorder1, TypeBorder2, TypeBorder3, TypeBorder4);
                                break;
                            case 11:
                                plast.StaticDecisionMethodKantorovichVlasovOneLevel(F, 0, 0, 1000, TypeBorder1, TypeBorder2, TypeBorder3, TypeBorder4);
                                break;
                            case 12:
                                plast.StaticDecisionMethodKantorovichVlasovOneLevel(F, 0, 1000, 1000, TypeBorder1, TypeBorder2, TypeBorder3, TypeBorder4);
                                break;
                            case 13:
                                plast.StaticDecisionMethodVariationIteration(F, 1000, 0, 0, TypeBorder1, TypeBorder2, TypeBorder3, TypeBorder4);
                                break;
                            default: break;
                        }

                        ///Вывод информации в файл)))))!!!!!!!!!!!!!!!!!!!!!
                        ///!!!!!!!!!!!!!!!!!!!!!!!!!!
                        ///!!!!!!!!!!!!!!!!!!!!!!!!!!!


                        // эпюр прогиба
                        for (int x = 1; x < N + 1; x++)
                            PrintWX[x - 1] = plast.W[x, M / 2];

                        //// эпюр моментов
                        //for (int x = 1; x < N + 1; x++)
                        //    PrintW2X[x - 1] = (plast.W[x - 1, M / 2] - 2 * plast.W[x, M / 2] + plast.W[x + 1, M / 2]) / (dx * dx);


                        //for (int y = 1; y < M + 1; y++)
                        //    PrintW2Y[y - 1] = (plast.W[N / 2, y - 1] - 2 * plast.W[N / 2, y] + plast.W[N / 2, y + 1]) / (dy * dy);

                        //Двойной интеграл
                        //double[] IntMasxX = new double[M];
                        //double Integral = 0;
                        //double dx = n / (N - 1);
                        //double dy = m / (M - 1);
                        //for (int y = 0; y < M; y++)
                        //{
                        //    for (int x = 0; x < N; x += 2)
                        //    {
                        //        if (N - x == 2)
                        //        {
                        //            IntMasxX[y] += (plast.W[N + 1, y + 1] + plast.W[N, y + 1]) * dx / 2;
                        //            break;
                        //        }
                        //        //формула Симсона
                        //        IntMasxX[y] += (plast.W[x + 1, y + 1] + 4 * plast.W[x + 2, y + 1] + plast.W[x + 3, y + 1]) * 2 * dx / 6;

                        //    }
                        //}
                        //for (int y = 0; y < M; y+=2)
                        //{
                        //    if (M - y == 2)
                        //    {
                        //        Integral += (IntMasxX[M - 1] + IntMasxX[M - 2]) * dy / 2;
                        //        break;
                        //    }
                        //    Integral+= (IntMasxX[y] + 4 * IntMasxX[y+1] + IntMasxX[y+2]) * 2 * dy / 6;
                        //}
                        
                        //WriteNumberInFile(Integral, "Integral" + Form);
                        // 

                        //    PrintW2Y[y - 1] = (plast.W[N / 2, y - 1] - 2 * plast.W[N / 2, y] + plast.W[N / 2, y + 1]) / (dy * dy);



                        WriteNumberInFile(plast.MaximumCountIterationOfMethodAllDecision, "IM_" + Form);

                        WriteNumberInFile(plast.CountIterationPhisicalNoneleneary, "FI_" + Form);

                        WriteNumberInFile(plast.mW, "mW_" + Form);

                        WriteMassivInFile(PrintWX, "PW" + Form);

                        //WriteNumberInFile(plast.Eii[14,14,0], "0.5_0.5_0" + Form);
                        //WriteNumberInFile(plast.Eii[6, 14, 0], "0.2_0.5_0" + Form);
                        //WriteNumberInFile(plast.Eii[6, 6, 0], "0.2_0.2_0" + Form);
                        //WriteNumberInFile(plast.Eii[14, 14, 1], "0.5_0.5_1" + Form);
                        //WriteNumberInFile(plast.Eii[6, 14, 1], "0.2_0.5_1" + Form);
                        //WriteNumberInFile(plast.Eii[6, 6, 1], "0.2_0.2_1" + Form);

                        // Вывод в консоле

                        Console.WriteLine(String.Format("{0} {1:0.00} {2} {3} {4} {5} {6}", i, Q * i, plast.mW, plast.CountIterationPhisicalNoneleneary,
                            plast.CountIterationABS, plast.MaxCountIterationSpendOfMetods, plast.MaximumCountIterationOfMethodAllDecision));


                        //Вывод связанный со временем


                        //WriteNumberInFile(plast.SignalMaxW[plast.SignalMaxW.Length-1], String.Format("TQ{0}_{1}{2}{3}{4}_B(a){5}M{6}n{7}Q{8}WX_{9}_{10}_{11}t0{12}t1{11}dis{12}", 
                        //    typeQ, TypeBorder1, TypeBorder2, TypeBorder3, TypeBorder4, a, TypeMethod, n, Q, N, M, 2 * P - 1, t0, t1, dissipation));


                        //WriteMassivInFile(plast.SignalMaxW, String.Format("TQ{0}X2{1}Y2{2}X1{3}Y1{4}B(a){5}M{6}n{7}Q{8}WX_{9}_{10}_{11}t0{12}t1{11}dis{12}", 
                        //    typeQ, TypeBorder1, TypeBorder2, TypeBorder3, TypeBorder4, a, TypeMethod, n, Q, N, M, 2 * P - 1, t0, t1, dissipation));

                        //Console.WriteLine(String.Format("{0} {1:0.0}", i, Q * i) + " " + plast.SignalMaxW[plast.SignalMaxW.Length-1]);

                        //Проверка на возникноверние пластического шарнира
                        if (FlagPlastSharn)
                        {
                            for (int i1 = 0; i1 < plast.Eii.GetLength(0); i1++)
                                for (int j1 = 0; j1 < plast.Eii.GetLength(1); j1++)
                                {
                                    flag1 = 0;
                                    for (int k1 = (plast.Eii.GetLength(2) + 1) / 2; k1 < plast.Eii.GetLength(2); k1++)
                                        if (plast.Eii[i1, j1, k1] >= plast.eS)
                                            flag1++;
                                    if (flag1 == (plast.Eii.GetLength(2) - 1) / 2 && FlagPlastSharn == true)
                                    {
                                        FlagPlastSharn = false;
                                        Console.WriteLine(Q * i);
                                        WriteNumberInFile(Q * i, plast.mW, "QPSh_" + Form);
                                        WriteMassivInFile(plast.Eii, 0, String.Format("EiiPSh_{0}_", 0) + Form);
                                        WriteMassivInFile(plast.Eii, sloi, String.Format("EiiPSh_{0}_", sloi) + Form);
                                        WriteMassivInFile(plast.Eii, P - 3, String.Format("EiiPSh_{0}_", P - 3) + Form);
                                        double[,] granb = new double[P * 2 - 1, N];
                                        for (int i2 = 0; i2 < P * 2 - 1; i2++)
                                            for (int k2 = 0; k2 < N; k2++)
                                                granb[i2, k2] = plast.Eii[k2, 0, i2];
                                        WriteMassivInFile(granb, "granb 0", true);
                                        MassivePointsFisicalNonlineary[2, countb] = Q * i;
                                        MassivePointsFisicalNonlineary[3, countb] = plast.mW;
                                    }
                                }
                        }
                        //Печать данных для определённой нагрузки
                        if (flagNagruzki)
                        {
                            WriteNumberInFile(Q * i, "Q_" + Form);
                            WriteNumberInFile(plast.mW, "Q_W" + Form);
                            WriteMassivInFile(plast.Eii, 0, String.Format("Eii_{0}_", 0) + Form);
                            WriteMassivInFile(plast.Eii, sloi, String.Format("Eii_{0}_", sloi) + Form);
                            WriteMassivInFile(plast.Eii, P - 2, String.Format("Eii_{0}_", P - 2) + Form);
                            WriteMassivInFile(plast.Eii, 2*P-2, String.Format("Eii_{0}_", 2 * P - 2) + Form);
                            WriteMassivInFile(plast.Eii, 2 * P - 2-sloi, String.Format("Eii_{0}_", 2 * P - 2 - sloi) + Form);
                            WriteMassivInFile(plast.Eii, P, String.Format("Eii_{0}_", P) + Form);
                            Console.WriteLine(Q * i);
                            double[,] granb = new double[P * 2 - 1, N];
                            for (int i2 = 0; i2 < P * 2 - 1; i2++)
                                for (int k2 = 0; k2 < N; k2++)
                                    granb[i2, k2] = plast.Eii[k2, 0, i2];
                            WriteMassivInFile(granb, "granb 0", true);
                            WriteMassivInFile(plast.PrintW, "Progib", true);
                        }
                        //Первое проявление нелинейности
                        if (FlagFirstNonlinear)
                        {
                            for (int i1 = 0; i1 < plast.Eii.GetLength(0); i1++)
                            {
                                for (int j1 = 0; j1 < plast.Eii.GetLength(1); j1++)
                                {
                                    for (int k1 = (plast.Eii.GetLength(2) + 1) / 2; k1 < plast.Eii.GetLength(2); k1++)
                                    {
                                        if (plast.Eii[i1, j1, k1] >= plast.eS && FlagFirstNonlinear == true)
                                        {
                                            FlagFirstNonlinear = false;
                                            Console.WriteLine(Q * i);
                                            WriteNumberInFile(Q * i, plast.mW, "QFN_" + Form);
                                            WriteMassivInFile(plast.Eii, 0, String.Format("EiiFN_{0}_", 0) + Form);
                                            WriteMassivInFile(plast.Eii, sloi, String.Format("EiiFN_{0}_", sloi) + Form);
                                            WriteMassivInFile(plast.Eii, P - 3, String.Format("EiiFN_{0}_", P - 3) + Form);
                                            MassivePointsFisicalNonlineary[0, countb] = Q * i;
                                            MassivePointsFisicalNonlineary[1, countb] = plast.mW;
                                        }
                                    }
                                }
                            }
                            
                        }
                        MasiiveOfList[countb].Add(plast.mW);
                        //условие выхода из цикла расчёта
                        if (plast.CountIterationPhisicalNoneleneary >= 1000 || plast.mW > 0.3 || (i != 0 && plast.mW == 0)) break;
                        
                        //plast.ShowMassiv2(plast.W);
                        //plast.ShowMassiv2(plast.IM);
                        //plast.ShowMassiv3CentrR3(plast.Eii, 0);
                        //plast.ShowMassiv3CentrR3(plast.Eii, 4);
                        //plast.RELoad();
                    }
                }
                stopwatch.Stop();
                Console.WriteLine("Затрачиваемое время " + stopwatch.ElapsedMilliseconds + " мс");
                Console.WriteLine(Form);

                Console.WriteLine(bet);
                Thread.Sleep(1000);
                countb++;
                
            }
            WriteMassivInFile(MassivePointsFisicalNonlineary, "FuullFisicleNonlineary",false);
            WriteMasssiveListInFile(MasiiveOfList);
            //RasspredDiform(plast.Eii, 1, "1");
            Console.WriteLine("КОНЕЦ ПРОГРАММЫ");
                Console.ReadKey();
            
        }
    }
}