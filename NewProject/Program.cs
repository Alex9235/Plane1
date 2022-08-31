using Microsoft.VisualBasic;
using System;
using System.Globalization;
using System.IO;
using System.Linq;
using System.Linq.Expressions;
using System.Diagnostics;

namespace Decoder
{
    class Program
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
        public class Plast
        {
            //количество точек
            //по оси x
            public int N { get; }
            //по оси y
            public int M { get; }
            //по оси z
            public int P { get; }

            public double n { get; }
            public double m { get; }
            public double p { get; }

            //шаг разбиение 
            //по оси x
            protected double dx { get; }
            //по оси y
            protected double dy { get; }
            //по оси z
            protected double dz { get; }
            public double[,] DZ { get; }

            //шаг по времени
            protected double Deltat;

            public double SigmS;  //предел текучести
            public double eS; //деформация текучести

            protected double l; //Коэффициент связанный с нано

            private bool FlagMultiModularIntensiveDeformations;
            public double SigmS1;  //предел текучести учёт разномодульности
            public double eS1; //деформация текучести учёт разномодульности

            private bool FlagLinearyAndExponentModelSigmaEpsilon1;
            private bool FlagLinearyAndExponentModelSigmaEpsilon2;

            private bool FlagFisicalNonLineary;

            //добавление криволиненой поверхности
            private bool FlagСurvilinearPlane;
            private int ModelFunctionСurvilinearPlane;
            //учёт температуры
            private bool FlagTemperatureProblem;
            public double[,,] T { get; } // температуное поле
            private double QTemp; //внутренний источник тепла
            private double[,,] Temp; //массив внутренннего источника тепла
            private double alf; //коэффициент линейного температурного расширения
            //временная задача
            private bool FlagTimeProblem;
            private double Dissipation; //коэффициент диссипации
            //учёт градиента
            private bool FlagGradientProblem;
            private double ConstRaspredGradient;
            private double GradientE;
            private double GradientNU;
            //учёт пористости
            private bool FlagPorysotyProblem;
            private double ConstPoriststb;
            private double ModelPoriststb;

            public bool Error { get; private set; } // глобальная ошибка

            public double lam1 { get; } // отношение длины к ширине
            public double lam2 { get; } // отношение длины к толщине
            public double lam { get; } // отношение длины к толщине
            public double[,,] NU { get; } //Коэффициент Пуассона
            public double K1 { get; private set; }// коэффициент объемного растяжения 
            public double K2 { get; private set; }// коэффициент объемного сжатия 
            public double Gk1 { get; }// коэффициент линейности без учёта физ нелинейности расятжений
            public double Gk2 { get; private set; }// коэффициент линейности без учёта физ нелинейности сжатия
            
            //линенйное упрочнение
            private bool FlagLinearyHard1;
            private bool FlagLinearyHard2;
            public double G01 { get; private set; }
            public double G11 { get; private set; }

            public double G0 { get; private set; } // модуль сдвига 1   
            public double G1 { get; private set; } // модуль сдвига 2

            public double[,,] E { get; }// Модуль Юнга

            //double[,] F;//Усилие
            public double[,] W { get; }// Прогиб
            public double[,] PrintW { get; }// Для вывода

            // массив моментов    
            public double[,] DM { get; }
            public double[,] RM { get; }
            public double[,] SM { get; }
            public double[,] IM { get; }// Массив интеграллов с темпретурой 
            //.......................
            public double[,,] Eii { get; } // массив деформаций

            public double mW { get; set; } //максимальное значение прогиба
            private double nu;  //коэффициент пуассона

            //ИНИЦИАЛИЗАЦИИ Учёта усложнения модели
            public void InicializationPlast(double nu, double G0)
            {
                if ( G0 <= 0 || nu < 0)
                {
                    ErrorEnterConstructorParametrs();
                    return;
                }
                this.G0 = G0;
                this.nu = nu;
                this.K1 = 2 * (1 + nu) / (3 * (1 - 2 * nu));
            }
            public void InicializationExpFizicalNonlinProblem(double SigmaS0, double eS0)
            {
                if (SigmaS0 <= 0 || eS0 <= 0)
                {
                    ErrorEnterConstructorParametrs();
                    return;
                }

                this.SigmS = SigmaS0 * lam1 * lam1 / G0;
                this.eS = eS0 * lam1 * lam1;
                FlagLinearyAndExponentModelSigmaEpsilon1 = true;
                Console.WriteLine("OneModuleFizcalExpProblem");
            }
            public void InicializationMultimoduWithExpFizicalNonlinlProblem(double nu2, double G1,double SigmaS1, double eS1 )
            {
                if (nu2<=0 || G1<=0 || SigmaS1 <= 0 || eS1 <= 0)
                {
                    ErrorEnterConstructorParametrs();
                    return;
                }
                this.FlagLinearyAndExponentModelSigmaEpsilon2 = true;
                this.FlagMultiModularIntensiveDeformations = true;
                this.G1 = G1;
                this.SigmS1 = SigmaS1 * lam1 * lam1 / G0;
                this.eS1 = eS1 * lam1 * lam1;
                this.K2 = Gk2 * 2 * (1 + nu2) / (3 * (1 - 2 * nu2));
                Console.WriteLine("MultiModuleFizcalExpProblem");
            }
            public void InicializationMultimodulProblem(double nu2, double G1)
            {
                if (nu2 <= 0 || G1 <= 0 )
                {
                    ErrorEnterConstructorParametrs();
                    return;
                }
                this.FlagMultiModularIntensiveDeformations = true;
                this.G1 = G1;
                this.Gk2 = G1 / G0;
                this.K2 = Gk2 * 2 * (1 + nu2) / (3 * (1 - 2 * nu2));

            }
            public void InicializationNano(double l)
            {
                if (l <= 0 || l > 1)
                {
                    ErrorEnterConstructorParametrs();
                    return ;
                }
                this.l = l;
            }
            public void InicializationTempature(double QTemp, double alf,double left, double right, double top, double down, double front, double back)
            {
                if (FlagСurvilinearPlane)
                {
                    ErrorEnterConstructorParametrs();
                    Console.WriteLine("Temperature can't job whis CurvilinearPlane");
                    return;
                }
                if (alf<=0 || alf >1)
                {
                    ErrorEnterConstructorParametrs();
                    return;
                }
                this.FlagTemperatureProblem=true;
                this.alf = alf;
                double r = alf * lam1 * lam1;
                this.QTemp = QTemp * alf * lam1 * lam1;
                LoadTFULL(this.QTemp);
                BorderT(left * r, right * r, top * r,  down * r,  front * r, back * r);
                LoadT();
            }
            public void InicializationTempature(double QTemp, double alf)
            {
                if (FlagСurvilinearPlane)
                {
                    ErrorEnterConstructorParametrs();
                    Console.WriteLine("Temperature can't job whis CurvilinearPlane");
                    return;
                }
                if (alf <= 0 || alf > 1)
                {
                    ErrorEnterConstructorParametrs();
                    return;
                }
                this.FlagTemperatureProblem = true;
                this.alf = alf;
                double r = alf * lam1 * lam1;
                this.QTemp = QTemp * r;
                LoadTFULL(this.QTemp);
                BorderT(this.QTemp, this.QTemp, this.QTemp, this.QTemp, this.QTemp, this.QTemp);
            }  
            public void InicializationTempature(String puth, double alf)
            {
                // 
                if (FlagСurvilinearPlane)
                {
                    ErrorEnterConstructorParametrs();
                    Console.WriteLine("Temperature can't job whis CurvilinearPlane");
                    return;
                }
                if (alf <= 0 || alf > 1)
                {
                    ErrorEnterConstructorParametrs();
                    return;
                }

                this.FlagTemperatureProblem = true;
                this.alf = alf;
                double r = alf * lam1 * lam1;

                System.Threading.Thread.CurrentThread.CurrentCulture = new System.Globalization.CultureInfo("en-US");// смена точки на запятую
                FileStream file1 = new FileStream(puth, FileMode.Open); //создаем файловый поток
                StreamReader reader = new StreamReader(file1);//создаем «потоковый писатель» и связываем его с файловым потоком 
                string s= "";
                while (s != "% Data")
                    s = reader.ReadLine();
                int x;
                for (int z = 0; z < P; z++)
                
                    for (int y = 0; y < M; y++)
                    {
                        s = reader.ReadLine();
                        x = 0;
                        foreach (var number in s.Split())
                        {
                            if (number != "")
                            {
                                T[x, y, z] = Convert.ToDouble(number) *r;
                                x++;
                            }
                        }
                    }
                reader.Close(); //закрываем поток. Не закрыв поток, в файл ничего не запишется так как первичная запись идёт в ОЗУ
            }
            public void InicializationTempatureLevelTwo(double[,] Temp)
            {
                //задаётся температурное распределение через массив
                //  что делать с коэффициентом альфа
                //задаётся связанность зависимости интенсивности напряжения от деформации

                if (FlagСurvilinearPlane)
                {
                    ErrorEnterConstructorParametrs();
                    Console.WriteLine("Temperature can't job whis CurvilinearPlane");
                    return;
                }
                if (alf <= 0 || alf > 1)
                {
                    ErrorEnterConstructorParametrs();
                    return;
                }
                this.FlagTemperatureProblem = true;
                this.alf = alf;
                this.QTemp = QTemp * alf * lam1 * lam1;
                LoadTFULL(this.QTemp);
                LoadT();
            }

            public void InicializationTimeProblem(double Dissipation)
            {
                if (Dissipation < 0)
                {
                    ErrorEnterConstructorParametrs();
                    return;
                }
                FlagTimeProblem = true;
            }
            public void InicializationVariationGradient(double ConstRaspredGradient, double NU1, double E1) 
            {
                if (ConstRaspredGradient <= 0)
                {
                    ErrorEnterConstructorParametrs();
                    return ;
                }
                FlagGradientProblem = true; 
                this.GradientE = E1 / G0;
                this.GradientNU = NU1;
            }
            public void InicializationPorisoty(double CountPorisoty, int ModelPorisoty) 
            {
                if (CountPorisoty <= 0 || CountPorisoty > 1 || ModelPorisoty < 1 || ModelPorisoty > 4)
                {
                    ErrorEnterConstructorParametrs();
                    return;
                }
                FlagPorysotyProblem = true; 
                this.ConstPoriststb = CountPorisoty;
                this.ModelPoriststb = ModelPorisoty;
            }
            public void InicializationCurvelinearyPlane(int ModelCurvilinearPlane)
            {

                if(ModelCurvilinearPlane <= 0 || ModelCurvilinearPlane > 3)
                {
                    ErrorEnterConstructorParametrs();
                    return;
                }
                this.ModelFunctionСurvilinearPlane = ModelCurvilinearPlane;
                this.FlagСurvilinearPlane = true;
                LoadDzForСurvilinearPlane();


            }
            public void InicializationLinearyHard(double G1, double es0) 
            {
                if (G1 <= 0 || es0 <= 0)
                {
                    ErrorEnterConstructorParametrs();
                    return;
                }
                FlagLinearyHard1 = true;
                this.G01 = G1;
                this.eS = es0 * lam1 * lam1; ;
            }


            //конструктор Пластинки
            public Plast(int N, int M, int P, double n, double m, double p)//
            {

                //количество точек разбиения
                this.N = N;
                this.M = M;
                this.P = 2 * P - 1;
                //длины
                this.n = n;
                this.m = m;
                this.p = p;

                if (N < 2 || M < 2 || P < 2 || n <= 0 || m <= 0 || p <= 0)
                {
                    ErrorEnterConstructorParametrs();
                    return;
                }
                //массивы физических параметров
                // -1,0....n,n+1 ( учёт счёта производной)
                E = new double[N, M, this.P];
                NU = new double[N, M, this.P];
                DM = new double[N + 2, M + 2];
                RM = new double[N + 2, M + 2];
                SM = new double[N + 2, M + 2];
                IM = new double[N + 2, M + 2];
                W = new double[N + 2, M + 2];
                PrintW = new double[N, M];
                Eii = new double[N, M, this.P];
                T = new double[N, M, this.P];
                Temp = new double[N, M, this.P];
                DZ = new double[N, M];

                //геометрические параметры
                lam = n / m;
                lam1 = n / p;
                lam2 = m / p;
                dx = 1.0f / (N - 1);
                dy = 1.0f / (M - 1);
                dz = 2.0f / (this.P - 1);
                LoadDZ(dz);
                this.FlagFisicalNonLineary = true;
                this.FlagLinearyAndExponentModelSigmaEpsilon1 = false;
                this.FlagLinearyAndExponentModelSigmaEpsilon2 = false;
                this.FlagMultiModularIntensiveDeformations = false;

                // физические константы
                this.G0 = 0;
                this.nu = 0;
                this.SigmS = 0;
                this.eS = 0;
                this.Gk1 = 1; // потому что при учёте безразмерности 
                this.K1 = 0;
                this.Gk2 = 0;
                this.K2 = 0;
                // учёт нано
                this.l = 0;

                //cвязанность с температурой
                this.FlagTemperatureProblem = false;
                this.alf = 0;
                this.QTemp = 0;
                //учёты распределения градиента
                this.FlagGradientProblem = false;   
                this.ConstRaspredGradient = 0;
                this.GradientE = 0;
                this.GradientNU = 0;

                mW = 0; //максимальный прогиб
                

                //Временной параметр
                this.FlagTimeProblem = false;
                //Пористость
                this.FlagPorysotyProblem = false;
                this.ConstPoriststb = 0;
                this.ModelPoriststb = 0;
                //Криволинейная поверхность пластины в расположении top и down
                this.FlagСurvilinearPlane = false;
                this.ModelFunctionСurvilinearPlane = 0;
                //линейное упрочнение
                this.FlagLinearyHard1 = false;
                this.FlagLinearyHard2 = false; //!!!! пока не написано
                this.G01 = 0;
                this.G11 = 0;//!!!! пока не написано
            }

            private void ErrorEnterConstructorParametrs()
            {
                Error = true;
                Console.WriteLine("ErorEnterParameters!!!!!");
            }
            ///Методы связанные с физической нелинейностью пластинки пластинки
            //Вывод на консоль
            public void ShowMassiv(double[] A)
            {
                Console.WriteLine();
                for (int i = 0; i < A.Count(); i++)
                    Console.WriteLine(A[i]);
                Console.WriteLine();
            }
            public void ShowMassiv2(double[,] A)
            {
                Console.WriteLine();
                for (int j = A.GetLength(1) - 1; j >= 0; j--)
                {
                    for (int i = 0; i < A.GetLength(0); i++)
                        Console.Write("{0:0.00000000} ", A[i, j]);
                    Console.WriteLine();
                }
                Console.WriteLine();
            }

            public void ShowMassiv3CentrR3(double[,,] A, int z)
            {
                Console.WriteLine();
                for (int j = A.GetLength(0) - 1; j >= 0; j--)
                {
                    for (int i = 0; i < A.GetLength(1); i++)
                        Console.Write(A[i, j, z] + " ");
                    Console.WriteLine();
                }
                Console.WriteLine();
            }
            //всё что связанно с криволинейностью внешней поверхности
            //функция криволиненой поверхности
            private double FunctionСurvilinearPlane(double x, double y)
            {
                //максимум 0.125
                switch(ModelFunctionСurvilinearPlane)
                {
                    case 1:
                        return -(Math.Pow(Math.Cos(x * Math.PI) * Math.Cos(y * Math.PI), 2) - 1) / 8;
                    case 2:
                        return Math.Pow(Math.Cos(x * Math.PI) * Math.Cos(y * Math.PI), 2) / 8;
                    case 3:
                        return -(Math.Pow(Math.Sin((x-y) * Math.PI) * Math.Sin((x+y) * Math.PI), 2) - 1) / 8;
                    default:
                        return 0;

                }
            }
            //расчёт всех dz зависимых от x и y
            private void LoadDZ(double Dz)
            {
                for (int x = 0; x < N; x++)
                    for (int y = 0; y < M; y++)
                        DZ[x, y] = Dz;
            }
            private void LoadDzForСurvilinearPlane()
            {
                for (int x = 0; x < N; x++)
                    for (int y = 0; y < M; y++)
                        DZ[x, y] = dz * (FunctionСurvilinearPlane(dx * x, dy * y) + 1);
            }
            // интенсивность деформации
            private double ei(int x, int y, int z, double d2wx, double d2wy, double d2wxy)
            {
                double lam11 = lam * lam;   
                if (FlagMultiModularIntensiveDeformations)
                    return (double)(2.0 * DZ[x,y] * (z - (P - 1) / 2) *
                        Math.Sqrt((NU[x, y, z] / Math.Pow(1 - NU[x, y, z], 2.0) + 1) *
                        Math.Pow(d2wx + lam11 * d2wy, 2.0) + 3 * lam11 *
                        (Math.Pow(d2wxy, 2.0) - d2wx * d2wy)) / 3);
                else
                    return (double)(2.0 * Math.Abs(DZ[x,y] * (z - (P - 1) / 2)) *
                        Math.Sqrt((NU[x, y, z] / Math.Pow(1 - NU[x, y, z], 2.0) + 1) *
                        Math.Pow(d2wx + lam11 * d2wy, 2.0) + 3 * lam11 *
                        (Math.Pow(d2wxy, 2.0) - d2wx * d2wy)) / 3);
                

            }
            // зависимости напряжений от деформаций
            private double SigmaEii(double eii)
            {
                if (eii < 0)
                    if (FlagLinearyAndExponentModelSigmaEpsilon2)
                        return (double)(-SigmS1 * (1 - Math.Exp(eii / eS1)));
                    else
                        return (double)(3 * Gk2 * eii);
                else
                    if (FlagLinearyAndExponentModelSigmaEpsilon1)
                        return (double)(SigmS * (1 - Math.Exp(-eii / eS)));
                    else
                        if (FlagLinearyHard1)
                            if (eii <= eS)
                                return 3 * eii;
                            else
                                return (eS*(G0-G01)+G01*eii)*3/G0;

                        else
                            return (double)(3 * eii);
            }
            private double SigmaEiiT(double eii,double T)
            {
                if (eii < 0)
                    if (FlagLinearyAndExponentModelSigmaEpsilon2)
                        return (double)(-SigmS1 * (1 - Math.Exp(eii / eS1)));
                    else
                        return (double)(3 * Gk2 * eii);
                else
                    if (FlagLinearyAndExponentModelSigmaEpsilon1)
                    return (double)(SigmS * (1 - Math.Exp(-eii / eS)));
                else
                        if (FlagLinearyHard1)
                    if (eii <= eS)
                        return 3 * eii;
                    else
                        return (eS * (G0 - G01) + G01 * eii) * 3 / G0;

                else
                    return (double)(3 * eii);
            }


            // подсчёт G
            private double G(double eii)
            {
                if (eii == 0)
                    return 1;
                else
                    return SigmaEii(eii) / (3 * eii);
            }
            private void LoadEii()
            {
                for (int y = 0; y < M; y++)
                {
                    for (int x = 0; x < N; x++)
                    {
                        double d2wx = (W[x, y + 1] - 2 * W[x + 1, y + 1] + W[x + 2, y + 1]) / (dx * dx);
                        double d2wy = (W[x + 1, y] - 2 * W[x + 1, y + 1] + W[x + 1, y + 2]) / (dy * dy);
                        double d2wxy = (W[x + 2, y + 2] - W[x + 2, y] - W[x, y + 2] + W[x, y]) / (4 * dx * dy);
                        for (int z = 0; z < P; z++)
                            Eii[x, y, z] = ei(x, y, z, d2wx, d2wy, d2wxy);
                    }
                }
            }
            // функция для оперделения пористости
            private double PoristFunctionToHeigth(double x,double y,double z)
            {
                switch (ModelPoriststb)
                {
                    case 1:
                        return 1 - ConstPoriststb;
                    case 2:
                        if (FlagСurvilinearPlane) return 1 - (Math.Abs(z))* ConstPoriststb / (1 + FunctionСurvilinearPlane(x, y));
                        return 1 - ConstPoriststb * Math.Abs(z);
                    case 3:
                        if (FlagСurvilinearPlane) return 1 - ((1 + FunctionСurvilinearPlane(x, y))- Math.Abs(z))*ConstPoriststb / (1 + FunctionСurvilinearPlane(x, y));
                        else return 1 - ConstPoriststb * (1 - Math.Abs(z));
                    case 4:
                        return 1 - ConstPoriststb * Math.Cos(z * Math.PI / 2);
                    default:
                        return 0;
                }
            }

            // расчёт массива коэффициентов Пуссона
            private void LoadNU()
            {
                for (int x = 0; x < N; x++)
                {
                    for (int y = 0; y < M; y++)
                    {
                        for (int z = 0; z < P; z++)
                        {
                            double Gg = G(Eii[x, y, z]);
                            if (Eii[x, y, z] < 0)
                                NU[x, y, z] = (3 * K2 - 2 * Gg) / (6 * K2 + 2 * Gg);
                            else
                                NU[x, y, z] = (3 * K1 - 2 * Gg) / (6 * K1 + 2 * Gg);
                        }
                    }
                }
            }
            private void PoritostbLoadNU()
            {
                for (int x = 0; x < N; x++)
                    for (int y = 0; y < M; y++)
                        for (int z = 0; z < P; z++)
                        {
                            NU[x, y, z] = NU[x,y,z]*PoristFunctionToHeigth(dx * x, dy * y, DZ[x, y] * (z - (P - 1) / 2));

                            //if (z == 0 && GradientNU != 0) NU[x, y, z] = GradientNU * PoristFunctionToHeigth(-1);

                            //if (z == P - 1) NU[x, y, z] = NU[x, y, z] * PoristFunctionToHeigth(1);
                        }
            }
            private void GradientLoadNU()
            {
                for (int x = 0; x < N; x++)
                    for (int y = 0; y < M; y++)
                        for (int z = 0; z < P; z++)
                        {
                            NU[x, y, z] = (GradientNU + (NU[x, y, z] - GradientNU) *
                                Math.Pow(1.0 / 2 + (DZ[x, y] * (z - ((P - 1) / 2)) / (2 * (1 + FunctionСurvilinearPlane(dx * x, dy * y)))), ConstRaspredGradient));
                        }
            }
            // расчёт массива Модулей Юнга
            private void LoadE()
            {
                for (int x = 0; x < N; x++)
                {
                    for (int y = 0; y < M; y++)
                    {
                        for (int z = 0; z < P; z++)
                        {
                            double Gg = G(Eii[x, y, z]);
                            if (Eii[x, y, z] < 0)
                                E[x, y, z] = 9 * K2 * Gg / (3 * K2 + Gg);
                            else
                                E[x, y, z] = 9 * K1 * Gg / (3 * K1 + Gg);
                        }
                    }
                }
            }

            private void GradientLoadE()
            {
                for (int x = 0; x < N; x++)
                    for (int y = 0; y < M; y++)
                        for (int z = 0; z < P; z++)
                        {
                            E[x, y, z] = (GradientE + (E[x, y, z] - GradientE) *
                                Math.Pow(1.0 / 2 + (DZ[x, y] * (z - ((P - 1) / 2)) / (2 * (1 + FunctionСurvilinearPlane(dx * x, dy * y)))), ConstRaspredGradient));
                            // G=E/2(1+v) E/G=2(1+v) 7.927519
                            //if (z == 0 && GradientE != 0) E[x, y, z] = GradientE * PoristFunctionToHeigth(-1);
                            //if (z == P - 1) E[x, y, z] = E[x, y, z] * PoristFunctionToHeigth(1);
                        }
            }
            private void PoritostbLoadE()
            {
                for (int x = 0; x < N; x++)
                    for (int y = 0; y < M; y++)
                        for (int z = 0; z < P; z++)
                        {
                            E[x, y, z] = E[x, y, z]*PoristFunctionToHeigth(dx * x, dy * y, DZ[x, y] * (z - (P - 1) / 2));
                            // G=E/2(1+v) E/G=2(1+v) 7.927519
                            //if (z == 0 && GradientE != 0) E[x, y, z] = GradientE * PoristFunctionToHeigth(-1);
                            //if (z == P - 1) E[x, y, z] = E[x, y, z] * PoristFunctionToHeigth(1);
                        }
            }

            // Метод Гаусса-Жордана
            private void MethodGordanGauss(double[,] A)
            {
                for (int j = 0; j < A.GetLength(0); j++) //столбцы
                {
                    double cA = A[j, j];
                    for (int i = j; i < A.GetLength(0); i++)//строки
                    {
                        if (i == j)
                        {
                            for (int k = j; k < A.GetLength(1); k++)
                            {
                                A[i, k] = A[i, k] / cA;// вся строка делится на диагональный элемент
                            }
                        }
                        else
                        {
                            if (A[i, j] != 0)
                            {
                                double Ak = A[i, j];
                                for (int k = j; k < A.GetLength(1); k++)
                                {
                                    A[i, k] = A[i, k] - Ak * A[j, k];
                                }
                            }
                        }
                    }
                }
                for (int j = A.GetLength(0) - 1; j > 0; j--) //столбцы
                {
                    double cA = A[j, j];
                    for (int i = j - 1; i >= 0; i--)//строки
                    {

                        if (A[i, j] != 0)
                        {
                            double Ak = A[i, j];
                            A[i, A.GetLength(0)] = A[i, A.GetLength(0)] - Ak * A[j, A.GetLength(0)];
                            A[i, j] = 0;
                        }

                    }
                }
            }

            // расчёт интегралов, метод сведения 3-х мерной задачи к двумерной 
            private void LoadD()
            {
                double S;
                //матрица для метода гаусса с кубической интерполяцие
                double[,] A = new double[4, 5];
                for (int x = 0; x < N; x++)
                    for (int y = 0; y < M; y++)
                    {
                        S = 0;
                        //Формула симсона
                        for (int z = 1; z < P; z = z + 2)
                            S = S + (E[x, y, z - 1] * ((z - (P - 1) / 2 - 1) * (z - (P - 1) / 2 - 1) / (1 - NU[x, y, z - 1] * NU[x, y, z - 1]) + l * l / (1 + NU[x, y, z - 1]))
                                     + 4 * E[x, y, z] * ((z - (P - 1) / 2) * (z - (P - 1) / 2) / (1 - NU[x, y, z] * NU[x, y, z]) + l * l / (1 + NU[x, y, z]))
                                     + E[x, y, z + 1] * ((z - (P - 1) / 2 + 1) * (z - (P - 1) / 2 + 1) / (1 - NU[x, y, z + 1] * NU[x, y, z + 1]) + l * l / (1 + NU[x, y, z + 1]))) * DZ[x, y] * DZ[x, y] * DZ[x,y] / 3;
                        DM[x + 1, y + 1] = S;
                    }

                // нахождение за граничных значений методом кубической интерполяции
                // по оси х
                //ShowMassiv2(DM);
                for (int x = 1; x < N + 1; x++)
                {
                    for (int i = 0; i < 4; i++)
                    {
                        //Заполнение массива А
                        A[i, 0] = 1;
                        A[i, 1] = (double)Math.Pow(i * dy, 1.0);
                        A[i, 2] = (double)Math.Pow(i * dy, 2.0);
                        A[i, 3] = (double)Math.Pow(i * dy, 3.0);
                        A[i, 4] = DM[x, i + 1];
                    }
                    MethodGordanGauss(A);
                    //Находим значение спомощью кубического уравнения
                    DM[x, 0] = (double)(A[3, 4] * (double)Math.Pow(-dy, 3.0) + A[2, 4] * (double)Math.Pow(-dy, 2.0) + A[1, 4] * (-dy) + A[0, 4]);

                    for (int i = 0; i < 4; i++)
                    {
                        //Заполнение массива А
                        A[i, 0] = 1;
                        A[i, 1] = (double)Math.Pow(1 - (i * dy), 1.0);
                        A[i, 2] = (double)Math.Pow(1 - (i * dy), 2.0);
                        A[i, 3] = (double)Math.Pow(1 - (i * dy), 3.0);
                        A[i, 4] = DM[x, M - i];
                    }
                    MethodGordanGauss(A);
                    DM[x, M + 1] = (double)(A[3, 4] * (double)Math.Pow(1 + dy, 3.0) + A[2, 4] * (double)Math.Pow(1 + dy, 2.0) + A[1, 4] * (1 + dy) + A[0, 4]);
                }
                //ShowMassiv2(DM);
                // по оси y
                for (int y = 1; y < M + 1; y++)
                {
                    for (int i = 0; i < 4; i++)
                    {
                        //Заполнение массива А
                        A[i, 0] = 1;
                        A[i, 1] = (double)Math.Pow(i * dx, 1.0);
                        A[i, 2] = (double)Math.Pow(i * dx, 2.0);
                        A[i, 3] = (double)Math.Pow(i * dx, 3.0);
                        A[i, 4] = DM[i + 1, y];
                    }
                    MethodGordanGauss(A);
                    //Находим значение спомощью кубического уравнения
                    DM[0, y] = (double)(A[3, 4] * (double)Math.Pow(-dx, 3.0) + A[2, 4] * (double)Math.Pow(-dx, 2.0) + A[1, 4] * (-dx) + A[0, 4]);

                    for (int i = 0; i < 4; i++)
                    {
                        //Заполнение массива А
                        A[i, 0] = 1;
                        A[i, 1] = (double)Math.Pow(1 - (i * dx), 1.0);
                        A[i, 2] = (double)Math.Pow(1 - (i * dx), 2.0);
                        A[i, 3] = (double)Math.Pow(1 - (i * dx), 3.0);
                        A[i, 4] = DM[N - i, y];
                    }
                    MethodGordanGauss(A);
                    DM[N + 1, y] = (double)(A[3, 4] * (double)Math.Pow(dx * N, 3.0) + A[2, 4] * (double)Math.Pow(dx * N, 2.0) + A[1, 4] * dx * N + A[0, 4]);
                }
                //ShowMassiv2(DM);
            }
            private void LoadR()
            {
                double S;
                //матрица для метода гаусса с кубической интерполяцие
                double[,] A = new double[4, 5];
                for (int x = 0; x < N; x++)
                    for (int y = 0; y < M; y++)
                    {
                        S = 0;
                        //Формула симсона
                        for (int z = 1; z < P; z = z + 2)
                            S = S + (E[x, y, z - 1] * (NU[x, y, z - 1] * (z - (P - 1) / 2 - 1) * (z - (P - 1) / 2 - 1) / (1 - NU[x, y, z - 1] * NU[x, y, z - 1]) + l * l / (1 + NU[x, y, z - 1]))
                                     + 4 * E[x, y, z] * (NU[x, y, z] * (z - (P - 1) / 2) * (z - (P - 1) / 2) / (1 - NU[x, y, z] * NU[x, y, z]) + l * l / (1 + NU[x, y, z]))
                                     + E[x, y, z + 1] * (NU[x, y, z + 1] * (z - (P - 1) / 2 + 1) * (z - (P - 1) / 2 + 1) / (1 - NU[x, y, z + 1] * NU[x, y, z + 1]) + l * l / (1 + NU[x, y, z + 1])))
                                     * DZ[x, y] * DZ[x, y] * DZ[x, y] / 3;
                        RM[x + 1, y + 1] = S;
                    }
                // нахождение за граничных значений методом кубической интерполяции
                // по оси х
                for (int x = 1; x < N + 1; x++)
                {
                    for (int i = 0; i < 4; i++)
                    {
                        //Заполнение массива А
                        A[i, 0] = 1;
                        A[i, 1] = (double)Math.Pow(i * dy, 1.0);
                        A[i, 2] = (double)Math.Pow(i * dy, 2.0);
                        A[i, 3] = (double)Math.Pow(i * dy, 3.0);
                        A[i, 4] = RM[x, i + 1];
                    }
                    MethodGordanGauss(A);
                    //Находим значение спомощью кубического уравнения
                    RM[x, 0] = (double)(A[3, 4] * (double)Math.Pow(-dy, 3.0) + A[2, 4] * (double)Math.Pow(-dy, 2.0) + A[1, 4] * (-dy) + A[0, 4]);

                    for (int i = 0; i < 4; i++)
                    {
                        //Заполнение массива А
                        A[i, 0] = 1;
                        A[i, 1] = (double)Math.Pow(1 - (i * dy), 1.0);
                        A[i, 2] = (double)Math.Pow(1 - (i * dy), 2.0);
                        A[i, 3] = (double)Math.Pow(1 - (i * dy), 3.0);
                        A[i, 4] = RM[x, M - i];
                    }
                    MethodGordanGauss(A);
                    RM[x, M + 1] = (double)(A[3, 4] * (double)Math.Pow(1 + dy, 3.0) + A[2, 4] * (double)Math.Pow(1 + dy, 2.0) + A[1, 4] * (1 + dy) + A[0, 4]);
                }
                //ShowMassiv2(DM);
                // по оси y
                for (int y = 1; y < M + 1; y++)
                {
                    for (int i = 0; i < 4; i++)
                    {
                        //Заполнение массива А
                        A[i, 0] = 1;
                        A[i, 1] = (double)Math.Pow(i * dx, 1.0);
                        A[i, 2] = (double)Math.Pow(i * dx, 2.0);
                        A[i, 3] = (double)Math.Pow(i * dx, 3.0);
                        A[i, 4] = RM[i + 1, y];
                    }
                    MethodGordanGauss(A);
                    //Находим значение спомощью кубического уравнения
                    RM[0, y] = (double)(A[3, 4] * (double)Math.Pow(-dx, 3.0) + A[2, 4] * (double)Math.Pow(-dx, 2.0) + A[1, 4] * (-dx) + A[0, 4]);

                    for (int i = 0; i < 4; i++)
                    {
                        //Заполнение массива А
                        A[i, 0] = 1;
                        A[i, 1] = (double)Math.Pow(1 - (i * dx), 1.0);
                        A[i, 2] = (double)Math.Pow(1 - (i * dx), 2.0);
                        A[i, 3] = (double)Math.Pow(1 - (i * dx), 3.0);
                        A[i, 4] = RM[N - i, y];
                    }
                    MethodGordanGauss(A);
                    RM[N + 1, y] = (double)(A[3, 4] * (double)Math.Pow(dx * N, 3.0) + A[2, 4] * (double)Math.Pow(dx * N, 2.0) + A[1, 4] * dx * N + A[0, 4]);
                }


            }
            private void LoadS()
            {
                double S;
                //матрица для метода гаусса с кубической интерполяцие
                double[,] A = new double[4, 5];
                for (int x = 0; x < N; x++)
                    for (int y = 0; y < M; y++)
                    {
                        S = 0;
                        //Формула симсона
                        for (int z = 1; z < P; z = z + 2)
                            S = S + (E[x, y, z - 1] * ((z - (P - 1) / 2 - 1) * (z - (P - 1) / 2 - 1) + l * l) / (1 + NU[x, y, z - 1])
                                     + 4 * E[x, y, z] * ((z - (P - 1) / 2) * (z - (P - 1) / 2) + l * l) / (1 + NU[x, y, z])
                                     + E[x, y, z + 1] * ((z - (P - 1) / 2 + 1) * (z - (P - 1) / 2 + 1) + l * l) / (1 + NU[x, y, z + 1]))*DZ[x, y] * DZ[x, y] * DZ[x, y] / 3;
                        SM[x + 1, y + 1] = S;
                    }
                // нахождение за граничных значений методом кубической интерполяции
                // по оси х
                for (int x = 1; x < N + 1; x++)
                {
                    for (int i = 0; i < 4; i++)
                    {
                        //Заполнение массива А
                        A[i, 0] = 1;
                        A[i, 1] = (double)Math.Pow(i * dy, 1.0);
                        A[i, 2] = (double)Math.Pow(i * dy, 2.0);
                        A[i, 3] = (double)Math.Pow(i * dy, 3.0);
                        A[i, 4] = SM[x, i + 1];
                    }
                    MethodGordanGauss(A);
                    //Находим значение спомощью кубического уравнения
                    SM[x, 0] = (double)(A[3, 4] * (double)Math.Pow(-dy, 3.0) + A[2, 4] * (double)Math.Pow(-dy, 2.0) + A[1, 4] * (-dy) + A[0, 4]);

                    for (int i = 0; i < 4; i++)
                    {
                        //Заполнение массива А
                        A[i, 0] = 1;
                        A[i, 1] = (double)Math.Pow(1 - (i * dy), 1.0);
                        A[i, 2] = (double)Math.Pow(1 - (i * dy), 2.0);
                        A[i, 3] = (double)Math.Pow(1 - (i * dy), 3.0);
                        A[i, 4] = SM[x, M - i];
                    }
                    MethodGordanGauss(A);
                    SM[x, M + 1] = (double)(A[3, 4] * (double)Math.Pow(1 + dy, 3.0) + A[2, 4] * (double)Math.Pow(1 + dy, 2.0) + A[1, 4] * (1 + dy) + A[0, 4]);
                }
                //ShowMassiv2(DM);
                // по оси y
                for (int y = 1; y < M + 1; y++)
                {
                    for (int i = 0; i < 4; i++)
                    {
                        //Заполнение массива А
                        A[i, 0] = 1;
                        A[i, 1] = (double)Math.Pow(i * dx, 1.0);
                        A[i, 2] = (double)Math.Pow(i * dx, 2.0);
                        A[i, 3] = (double)Math.Pow(i * dx, 3.0);
                        A[i, 4] = SM[i + 1, y];
                    }
                    MethodGordanGauss(A);
                    //Находим значение спомощью кубического уравнения
                    SM[0, y] = (double)(A[3, 4] * (double)Math.Pow(-dx, 3.0) + A[2, 4] * (double)Math.Pow(-dx, 2.0) + A[1, 4] * (-dx) + A[0, 4]);

                    for (int i = 0; i < 4; i++)
                    {
                        //Заполнение массива А
                        A[i, 0] = 1;
                        A[i, 1] = (double)Math.Pow(1 - (i * dx), 1.0);
                        A[i, 2] = (double)Math.Pow(1 - (i * dx), 2.0);
                        A[i, 3] = (double)Math.Pow(1 - (i * dx), 3.0);
                        A[i, 4] = SM[N - i, y];
                    }
                    MethodGordanGauss(A);
                    SM[N + 1, y] = (double)(A[3, 4] * (double)Math.Pow(dx * N, 3.0) + A[2, 4] * (double)Math.Pow(dx * N, 2.0) + A[1, 4] * dx * N + A[0, 4]);
                }
                //ShowMassiv2(SM);
            }
            private void LoadI()
            {
                double S;
                //матрица для метода гаусса с кубической интерполяцие
                double[,] A = new double[4, 5];
                for (int x = 0; x < N; x++)
                    for (int y = 0; y < M; y++)
                    {
                        S = 0;
                        //Формула симсона
                        for (int z = 1; z < P; z = z + 2)
                            S = S + (E[x, y, z - 1] * T[x, y, z - 1] * (z - (P - 1) / 2 - 1) / (1 - NU[x, y, z - 1])
                                     + 4 * E[x, y, z] * T[x, y, z] * (z - (P - 1) / 2) / (1 - NU[x, y, z])
                                     + E[x, y, z + 1] * T[x, y, z + 1] * (z - (P - 1) / 2 + 1) / (1 - NU[x, y, z + 1])) * DZ[x, y] * DZ[x, y] / 3;
                        IM[x + 1, y + 1] = S;
                    }
                // нахождение за граничных значений методом кубической интерполяции
                // по оси х
                for (int x = 1; x < N + 1; x++)
                {
                    for (int i = 0; i < 4; i++)
                    {
                        //Заполнение массива А
                        A[i, 0] = 1;
                        A[i, 1] = (double)Math.Pow(i * dy, 1.0);
                        A[i, 2] = (double)Math.Pow(i * dy, 2.0);
                        A[i, 3] = (double)Math.Pow(i * dy, 3.0);
                        A[i, 4] = IM[x, i + 1];
                    }
                    MethodGordanGauss(A);
                    //Находим значение спомощью кубического уравнения
                    IM[x, 0] = (double)(A[3, 4] * (double)Math.Pow(-dy, 3.0) + A[2, 4] * (double)Math.Pow(-dy, 2.0) + A[1, 4] * (-dy) + A[0, 4]);

                    for (int i = 0; i < 4; i++)
                    {
                        //Заполнение массива А
                        A[i, 0] = 1;
                        A[i, 1] = (double)Math.Pow(1 - (i * dy), 1.0);
                        A[i, 2] = (double)Math.Pow(1 - (i * dy), 2.0);
                        A[i, 3] = (double)Math.Pow(1 - (i * dy), 3.0);
                        A[i, 4] = IM[x, M - i];
                    }
                    MethodGordanGauss(A);
                    IM[x, M + 1] = (double)(A[3, 4] * (double)Math.Pow(1 + dy, 3.0) + A[2, 4] * (double)Math.Pow(1 + dy, 2.0) + A[1, 4] * (1 + dy) + A[0, 4]);
                }
                //ShowMassiv2(DM);
                // по оси y
                for (int y = 1; y < M + 1; y++)
                {
                    for (int i = 0; i < 4; i++)
                    {
                        //Заполнение массива А
                        A[i, 0] = 1;
                        A[i, 1] = (double)Math.Pow(i * dx, 1.0);
                        A[i, 2] = (double)Math.Pow(i * dx, 2.0);
                        A[i, 3] = (double)Math.Pow(i * dx, 3.0);
                        A[i, 4] = IM[i + 1, y];
                    }
                    MethodGordanGauss(A);
                    //Находим значение спомощью кубического уравнения
                    IM[0, y] = (double)(A[3, 4] * (double)Math.Pow(-dx, 3.0) + A[2, 4] * (double)Math.Pow(-dx, 2.0) + A[1, 4] * (-dx) + A[0, 4]);

                    for (int i = 0; i < 4; i++)
                    {
                        //Заполнение массива А
                        A[i, 0] = 1;
                        A[i, 1] = (double)Math.Pow(1 - (i * dx), 1.0);
                        A[i, 2] = (double)Math.Pow(1 - (i * dx), 2.0);
                        A[i, 3] = (double)Math.Pow(1 - (i * dx), 3.0);
                        A[i, 4] = IM[N - i, y];
                    }
                    MethodGordanGauss(A);
                    IM[N + 1, y] = (double)(A[3, 4] * (double)Math.Pow(dx * N, 3.0) + A[2, 4] * (double)Math.Pow(dx * N, 2.0) + A[1, 4] * dx * N + A[0, 4]);
                }
                //ShowMassiv2(IM);
            }

            //граничные условия распределение температуры (учитываются в LoadT())
            private void BorderT(double left, double right, double top, double down, double front, double back)
            {
                //левая грань
                for (int y = 0; y < M; y++)
                    for (int z = 0; z < P; z++)
                        T[0, y, z] = left;
                //правая грань
                for (int y = 0; y < M; y++)
                    for (int z = 0; z < P; z++)
                        T[N - 1, y, z] = right;
                //верхняя грань
                for (int y = 0; y < M; y++)
                    for (int x = 0; x < N; x++)
                        T[x, y, P - 1] = top;
                //нижняя грань
                for (int y = 0; y < M; y++)
                    for (int x = 0; x < N; x++)
                        T[x, y, 0] = down;
                //фронт 
                for (int x = 0; x < N; x++)
                    for (int z = 0; z < P; z++)
                        T[x, 0, z] = front;
                //задняя грань
                for (int x = 0; x < N; x++)
                    for (int z = 0; z < P; z++)
                        T[x, M - 1, z] = back;
            }
            // уравнение для расчёта точек температуры 7 точечного шаблона
            private double YravnenieForT(double[,,] A, int x, int y, int z)
            {
                double del = 2 / (dx * dx) + 2 * lam * lam / (dy * dy) + 2 * lam1 * lam1 / (dz * dz);
                return ((A[x + 1, y, z] + A[x - 1, y, z]) / (dx * dx) +
                    lam * lam * (A[x, y + 1, z] + A[x, y - 1, z]) / (dy * dy) +
                    lam1 * lam1 * (A[x, y, z + 1] + A[x, y, z - 1]) / (dz * dz) + Temp[x, y, z]) / del;
            }

            // Расчёт распределения температуры при её учёте в зависимости
            //от граничных условий с точностью до 6 знака
            private void LoadT()
            {
                double[,,] Tf = new double[N, M, P];

                //for (int x = 0; x < N; x++)
                //    for (int y = 0; y < M; y++)
                //        for (int z = 0; z < P; z++)
                //        {
                //            T[x, y, z] = 0;
                //            //if((x == 10 || x==11) && (y == 10 || y == 11) && z == 3) Temp[x, y, z] = QTemp;
                //            //else Temp[x, y, z] = 0;
                //            Temp[x, y, z] = QTemp;
                //        }
                for (int x = 0; x < N; x++)
                    for (int y = 0; y < M; y++)
                        for (int z = 0; z < P; z++)
                            Tf[x, y, z] = T[x, y, z];
                while (true)
                {
                    double a = 0, b = 0, maxd = 0;
                    for (int x = 1; x < N - 1; x++)
                        for (int y = 1; y < M - 1; y++)
                            for (int z = 1; z < P - 1; z++)
                            {
                                a = T[x, y, z];
                                Tf[x, y, z] = YravnenieForT(T, x, y, z);
                                b = Tf[x, y, z];
                                if (Math.Abs(b - a) > maxd) maxd = Math.Abs(b - a);
                            }
                    for (int x = 0; x < N; x++)
                        for (int y = 0; y < M; y++)
                            for (int z = 0; z < P; z++)
                                T[x, y, z] = Tf[x, y, z];
                    if (maxd < 0.000001) break;
                    //ShowMassiv2(T);
                }
            }
            private void LoadTFULL(double ConstTemp)
            {
                for (int x = 0; x < N; x++)
                    for (int y = 0; y < M; y++)
                        for (int z = 0; z < P; z++)
                            T[x, y, z] = ConstTemp;
            }
            private void LoadTFULL(double[,,] ConstTemp)
            {
                for (int x = 0; x < N; x++)
                    for (int y = 0; y < M; y++)
                        for (int z = 0; z < P; z++)
                            T[x, y, z] = ConstTemp[x,y,z];
            }
            // Заполнение W Начальными значениями
            private void LoadW()
            {
                for (int x = 0; x < W.GetLength(0); x++)
                    for (int y = 0; y < W.GetLength(1); y++)
                        W[x, y] = 0;
            }
            // загрузка начальных значений
            public void RELoad()
            {

                //обнуление прогибов
                LoadW();
                // для начала забиваем коэффициент пуассона 
                for (int x = 0; x < N; x++)
                    for (int y = 0; y < M; y++)
                        for (int z = 0; z < P; z++)
                            NU[x, y, z] = nu;
                //вычисление деформции
                LoadEii();
                //вычисление модулей юнга 
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
                //Обнуление счётчика
                CountIterationPhisicalNoneleneary = 0;
                MaxCountMethodOfIteration = 0;
            }
            // загрузка физ параметров относительно установленного прогиба.
            private void Load()
            {
                //вычисление деформции
                LoadEii();
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
            }
            ///.............................................................

            ///Вызов метода с различным типом физической нелинейности...............
            //переменные производных для рассчёта методом.
            private double
                lamd2, lamd4,
                D, R, S,
                Ddx, Ddy, Dd2x, Dd2y,
                Rd2x, Rd2y, Rdx, Rdy,
                Sdx, Sdy, S2dxdy,
                Id2x = 0, Id2y = 0;
            private double
                Wd4x, Wd4y, Wd3x, Wd3y, Wd2x, Wd2y,
                Wd2xd2y, Wd2xdy, Wd2ydx, Wdydx;

            private int DegreeOfDifficFunctionMethodBubnovaGalercin = 0; //Сложность вычисления метода Бубного-Галёркина
            private int CountIterationInMethods = 0;// Счётчик метода итерационного метода
            public int CountIterationABS = 0;
            public int MaxCountIterationSpendOfMetods = 0; // Максимальное количество итераций затриченных методом
            public int MaximumCountIterationOfMethodAllDecision = 0;
            private double BorderExitIterationMethods = 0; //граница выхода метода вариационных итераций
            private double BorderExitMethodABS = 0; //граница выхода метода Аграновского Баглая Смирнова
            private double BorderExitPhisicalIteration = 0; //граница выхода подсчёта упрого-пластических дифформаций

            //для проверки выхода через максимум для итерационных методов
            private double MaxWForMethodIteration1;
            private double MaxWForMethodIteration2;
            // Делегаты для смены метода расчёта уравнения и исаользования в его физической нелинейности
            delegate void DelegatForMethods(Load F, int TypeBorder1, int TypeBorder2, int TypeBorder3, int TypeBorder4);

            /// ...................................................................................................... 

            ///Динамическая задача методом рунгекутта 4 порядка....................................................................................................
            //Пересчёт граничных условий для методов
            private void BorderZadelkSharnForMethods(double[,] A, int TypeBorder1, int TypeBorder2, int TypeBorder3, int TypeBorder4)
            {

                //Верхняя 
                switch (TypeBorder1)
                {
                    case 0:
                        for (int x = 0; x < A.GetLength(0); x++)
                        {
                            //Граничные условия жёсткой заделки
                            A[x, A.GetLength(1) - 2] = 0;
                            A[x, A.GetLength(1) - 1] = A[x, A.GetLength(1) - 3];
                        }
                        break;
                    case 1:
                        for (int x = 0; x < A.GetLength(0); x++)
                        {
                            //Граничные условия Шарнирного опирания
                            A[x, A.GetLength(1) - 2] = 0;
                            A[x, A.GetLength(1) - 1] = -A[x, A.GetLength(1) - 3];
                        }
                        break;
                    case 2:

                        for (int x = 0; x < A.GetLength(0); x++)
                        {
                            if (x < A.GetLength(0) / 2)
                            {
                                //Граничные условия жёсткой заделки
                                A[x, A.GetLength(1) - 2] = 0;
                                A[x, A.GetLength(1) - 1] = A[x, A.GetLength(1) - 3];
                            }
                            else
                            {
                                //Граничные условия Шарнирного опирания
                                A[x, A.GetLength(1) - 2] = 0;
                                A[x, A.GetLength(1) - 1] = -A[x, A.GetLength(1) - 3];
                            }
                        }
                        break;
                    case 3:
                        for (int x = 0; x < A.GetLength(0); x++)
                        {
                            if (x < A.GetLength(0) / 2)
                            {
                                //Граничные условия Шарнирного опирания
                                A[x, A.GetLength(1) - 2] = 0;
                                A[x, A.GetLength(1) - 1] = -A[x, A.GetLength(1) - 3];
                            }
                            else
                            {
                                //Граничные условия жёсткой заделки
                                A[x, A.GetLength(1) - 2] = 0;
                                A[x, A.GetLength(1) - 1] = A[x, A.GetLength(1) - 3];
                            }
                        }
                        break;
                    default:
                        for (int x = 0; x < A.GetLength(0); x++)
                        {
                            //Граничные условия жёсткой заделки
                            A[x, A.GetLength(1) - 2] = 0;
                            A[x, A.GetLength(1) - 1] = A[x, A.GetLength(1) - 3];
                        }
                        break;
                }
                //Нижняя
                switch (TypeBorder3)
                {
                    case 0:
                        for (int x = 0; x < A.GetLength(0); x++)
                        {
                            //Граничные условия жёсткой заделки
                            A[x, 1] = 0;
                            A[x, 0] = A[x, 2];
                        }
                        break;
                    case 1:
                        for (int x = 0; x < A.GetLength(0); x++)
                        {
                            //Граничные условия Шарнирного опирания
                            A[x, 1] = 0;
                            A[x, 0] = -A[x, 2];
                        }
                        break;
                    case 2:

                        for (int x = 0; x < A.GetLength(0); x++)
                        {
                            if (x < A.GetLength(0) / 2)
                            {
                                //Граничные условия жёсткой заделки
                                A[x, 1] = 0;
                                A[x, 0] = A[x, 2];
                            }
                            else
                            {
                                //Граничные условия Шарнирного опирания
                                A[x, 1] = 0;
                                A[x, 0] = -A[x, 2];
                            }
                        }
                        break;
                    case 3:
                        for (int x = 0; x < A.GetLength(0); x++)
                        {
                            if (x < A.GetLength(0) / 2)
                            {
                                //Граничные условия Шарнирного опирания
                                A[x, 1] = 0;
                                A[x, 0] = -A[x, 2];
                            }
                            else
                            {
                                //Граничные условия жёсткой заделки
                                A[x, 1] = 0;
                                A[x, 0] = A[x, 2];
                            }
                        }
                        break;
                    default:
                        for (int x = 0; x < A.GetLength(0); x++)
                        {
                            //Граничные условия жёсткой заделки
                            A[x, 1] = 0;
                            A[x, 0] = A[x, 2];
                        }
                        break;
                }
                //Правая
                switch (TypeBorder2)
                {
                    case 0:
                        for (int y = 0; y < A.GetLength(1); y++)
                        {
                            //Граничные условия жёсткой заделки
                            A[A.GetLength(0) - 2, y] = 0;
                            A[A.GetLength(0) - 1, y] = A[A.GetLength(0) - 3, y];
                        }
                        break;
                    case 1:
                        for (int y = 0; y < A.GetLength(1); y++)
                        {
                            //Граничные условия жёсткой заделки
                            A[A.GetLength(0) - 2, y] = 0;
                            A[A.GetLength(0) - 1, y] = -A[A.GetLength(0) - 3, y];
                        }
                        break;
                    case 2:
                        for (int y = 0; y < A.GetLength(1); y++)
                        {
                            if (y < A.GetLength(1) / 2)
                            {
                                //Граничные условия жёсткой заделки
                                A[A.GetLength(0) - 2, y] = 0;
                                A[A.GetLength(0) - 1, y] = A[A.GetLength(0) - 3, y];
                            }
                            else
                            {
                                //Граничные условия Шарнирного опирания
                                A[A.GetLength(0) - 2, y] = 0;
                                A[A.GetLength(0) - 1, y] = -A[A.GetLength(0) - 3, y];
                            }
                        }
                        break;
                    case 3:
                        for (int y = 0; y < A.GetLength(1); y++)
                        {
                            if (y < A.GetLength(1) / 2)
                            {
                                //Граничные условия Шарнирного опирания
                                A[A.GetLength(0) - 2, y] = 0;
                                A[A.GetLength(0) - 1, y] = -A[A.GetLength(0) - 3, y];
                            }
                            else
                            {
                                //Граничные условия жёсткой заделки
                                A[A.GetLength(0) - 2, y] = 0;
                                A[A.GetLength(0) - 1, y] = A[A.GetLength(0) - 3, y];
                            }
                        }
                        break;
                    default:
                        for (int y = 0; y < A.GetLength(1); y++)
                        {
                            //Граничные условия жёсткой заделки
                            A[A.GetLength(0) - 2, y] = 0;
                            A[A.GetLength(0) - 1, y] = A[A.GetLength(0) - 3, y];
                        }
                        break;
                }
                //Левая
                switch (TypeBorder4)
                {
                    case 0:
                        for (int y = 0; y < A.GetLength(1); y++)
                        {
                            //Граничные условия жёсткой заделки
                            A[1, y] = 0;
                            A[0, y] = A[2, y];
                        }
                        break;
                    case 1:
                        for (int y = 0; y < A.GetLength(1); y++)
                        {
                            //Граничные условия жёсткой заделки
                            A[1, y] = 0;
                            A[0, y] = -A[2, y];
                        }
                        break;
                    case 2:
                        for (int y = 0; y < A.GetLength(1); y++)
                        {
                            if (y < A.GetLength(1) / 2)
                            {
                                //Граничные условия жёсткой заделки
                                A[1, y] = 0;
                                A[0, y] = A[2, y];
                            }
                            else
                            {
                                //Граничные условия Шарнирного опирания
                                A[1, y] = 0;
                                A[0, y] = -A[2, y];
                            }
                        }
                        break;
                    case 3:
                        for (int y = 0; y < A.GetLength(1); y++)
                        {
                            if (y < A.GetLength(1) / 2)
                            {
                                //Граничные условия Шарнирного опирания
                                A[1, y] = 0;
                                A[0, y] = -A[2, y];
                            }
                            else
                            {
                                //Граничные условия жёсткой заделки
                                A[1, y] = 0;
                                A[0, y] = A[2, y];
                            }
                        }
                        break;
                    default:
                        for (int y = 0; y < A.GetLength(1); y++)
                        {
                            //Граничные условия жёсткой заделки
                            A[1, y] = 0;
                            A[0, y] = A[2, y];
                        }
                        break;
                }
            }

            //массив сигнала прогиба в центре
            public double[] SignalMaxW;
            //диапозон времени
            private double t0;
            private double t1;
            //для замены, равна первой производной по времени
            private double[,] r;
            // для метода рунгекутта
            private double[,] k1;
            private double[,] l1;
            private double[,] k2;
            private double[,] l2;
            private double[,] k3;
            private double[,] l3;
            private double[,] k4;
            private double[,] l4;
            private double[,] f1;
            private double[,] f2;

            public void DinamiDecisionMethodRugeKyta4(Load F, double t0, double t1, double Deltat, int BorderExitPhisicalIteration,
                int TypeBorder1, int TypeBorder2, int TypeBorder3, int TypeBorder4, int WriteShag)
            {
                // проверка на то что вызов конструктора был с учётом времени
                if (!FlagTimeProblem)
                {
                    Console.WriteLine("Error - FlagTimeProblem = false");
                    return;
                }

                RELoad();
                this.BorderExitPhisicalIteration = BorderExitPhisicalIteration;
                this.BorderExitMethodABS = 0;
                this.Deltat = Deltat;
                this.t0 = t0;
                this.t1 = t1;

                //размерность массива для заполнение сигнала
                SignalMaxW = new double[(int)((t1 - t0) / (Deltat * WriteShag)) + 1];
                //реализация массивов используемых в методе рунгекутта 4
                r = new double[N + 2, M + 2];
                k1 = new double[N + 2, M + 2];
                l1 = new double[N + 2, M + 2];
                k2 = new double[N + 2, M + 2];
                l2 = new double[N + 2, M + 2];
                k3 = new double[N + 2, M + 2];
                l3 = new double[N + 2, M + 2];
                k4 = new double[N + 2, M + 2];
                l4 = new double[N + 2, M + 2];
                f1 = new double[N + 2, M + 2];
                f2 = new double[N + 2, M + 2];

                for (int x = 0; x < N + 2; x++)
                    for (int y = 0; y < M + 2; y++)
                        r[x, y] = 0;

                int i = 0;
                int WriteSignal = 0;
                while (true)
                {
                    MethodRugeKyta4(F, TypeBorder1, TypeBorder2, TypeBorder3, TypeBorder4);
                    this.t0 += Deltat;

                    if (this.t0 > this.t1) break;
                    //Console.WriteLine("{0:0.0000} {1:0.00000000}", this.t0, W[W.GetLength(0) / 2 + 1, W.GetLength(1) / 2 + 1]);

                    if (WriteSignal == WriteShag)
                    {
                        Console.WriteLine("{0:0.0000} {1:0.00000000}", this.t0, W[W.GetLength(0) / 2 + 1, W.GetLength(1) / 2 + 1]);
                        SignalMaxW[i / WriteShag] = W[W.GetLength(0) / 2 + 1, W.GetLength(1) / 2 + 1]; // центр пластинки
                        WriteSignal = 0;
                    }
                    i++;
                    WriteSignal++;
                    //Пересчёт физических параметров
                    Load();
                }
                //Заполнения значиний прогибов для печати
                for (int x = 0; x < N; x++)
                    for (int y = 0; y < M; y++)
                        PrintW[x, y] = W[x + 1, y + 1];
            }

            // Замена и приведение к системе диф уравнений

            private double OperatorL(Load F, int x, int y, double[,] StateW, double t)
            {

                lamd2 = lam * lam;
                lamd4 = Math.Pow(lam, 4.0);
                D = DM[x, y]; R = RM[x, y]; S = SM[x, y];
                Ddx = (DM[x + 1, y] - DM[x - 1, y]) / (2 * dx);
                Ddy = (DM[x, y + 1] - DM[x, y - 1]) / (2 * dy);
                Dd2x = (DM[x + 1, y] - 2 * DM[x, y] + DM[x - 1, y]) / (dx * dx);
                Dd2y = (DM[x, y + 1] - 2 * DM[x, y] + DM[x, y - 1]) / (dy * dy);
                Rd2x = (RM[x + 1, y] - 2 * RM[x, y] + RM[x - 1, y]) / (dx * dx);
                Rd2y = (RM[x, y + 1] - 2 * RM[x, y] + RM[x, y - 1]) / (dy * dy);
                Rdx = (RM[x + 1, y] - RM[x - 1, y]) / (2 * dx);
                Rdy = (RM[x, y + 1] - RM[x, y - 1]) / (2 * dy);
                Sdx = (SM[x + 1, y] - SM[x - 1, y]) / (2 * dx);
                Sdy = (SM[x, y + 1] - SM[x, y - 1]) / (2 * dy);
                S2dxdy = (SM[x + 1, y + 1] - SM[x + 1, y - 1] - SM[x - 1, y + 1] + SM[x - 1, y - 1]) / (4 * dx * dy);

                Wd4x = (StateW[x + 2, y] - 4 * StateW[x + 1, y] + 6 * StateW[x, y] - 4 * StateW[x - 1, y] + StateW[x - 2, y]) / Math.Pow(dx, 4.0);
                Wd4y = (StateW[x, y + 2] - 4 * StateW[x, y + 1] + 6 * StateW[x, y] - 4 * StateW[x, y - 1] + StateW[x, y - 2]) / Math.Pow(dy, 4.0);
                Wd3x = (StateW[x + 2, y] - 2 * StateW[x + 1, y] + 2 * StateW[x - 1, y] - StateW[x - 2, y]) / (2 * Math.Pow(dx, 3.0));
                Wd3y = (StateW[x, y + 2] - 2 * StateW[x, y + 1] + 2 * StateW[x, y - 1] - StateW[x, y - 2]) / (2 * Math.Pow(dy, 3.0));
                Wd2x = (StateW[x + 1, y] - 2 * StateW[x, y] + StateW[x - 1, y]) / Math.Pow(dx, 2.0);
                Wd2y = (StateW[x, y + 1] - 2 * StateW[x, y] + StateW[x, y - 1]) / Math.Pow(dy, 2.0);
                Wd2xd2y = (StateW[x + 1, y + 1] + StateW[x + 1, y - 1] + StateW[x - 1, y + 1] + StateW[x - 1, y - 1]
                    - 2 * StateW[x, y + 1] - 2 * StateW[x, y - 1] - 2 * StateW[x - 1, y] - 2 * StateW[x + 1, y] + 4 * StateW[x, y]) / (dx * dx * dy * dy);
                Wd2ydx = (StateW[x + 1, y + 1] + StateW[x + 1, y - 1] - StateW[x - 1, y + 1] - StateW[x - 1, y - 1]
                    + 2 * StateW[x - 1, y] - 2 * StateW[x + 1, y]) / (2 * dx * dy * dy);
                Wd2xdy = (StateW[x + 1, y + 1] - StateW[x + 1, y - 1] + StateW[x - 1, y + 1] - StateW[x - 1, y - 1]
                    - 2 * StateW[x, y + 1] + 2 * StateW[x, y - 1]) / (2 * dx * dx * dy);
                Wdydx = (StateW[x + 1, y + 1] - StateW[x + 1, y - 1] - StateW[x - 1, y + 1] + StateW[x - 1, y - 1]) / (4 * dx * dy);

                x--;
                y--;
                if (F.FlagTimeProblem)
                    F.FullLoadSint(t);
                return F.F[x, y] - (Wd4x * D + 2 * Wd3x * Ddx + (Dd2x + lamd2 * Rd2y) * Wd2x +
                    lamd4 * (Wd4y * D + 2 * Wd3y * Ddy + (Dd2y + Rd2x / lamd2) * Wd2y) +
                    2 * lamd2 * (Wd2xd2y * (R + S) + Wd2ydx * (Rdx + Sdx) + Wd2xdy * (Rdy + Sdy) + Wdydx * S2dxdy));

                //return F.F[x, y] - D *(Wd4x + 2 * Wd2xd2y + Wd4y);

            }
            private double FunctionG1(Load F, int x, int y, double[,] StateW, double[,] State_r, double t)
            {
                return State_r[x, y];
            }
            private double FunctionG2(Load F, int x, int y, double[,] StateW, double[,] State_r, double t)
            {
                return OperatorL(F, x, y, StateW, t) - Dissipation * State_r[x, y];
            }
            private void MethodRugeKyta4(Load F, int TypeBorder1, int TypeBorder2, int TypeBorder3, int TypeBorder4)
            {
                //Уравнение 
                // d^2w/dt^2 + e*dw/dt = L(w,t)
                //Система дифференциальных уравнений
                // dw/dt=r - G1
                // dr/dt=L(w,t)-er -G2

                // процедура метода рунгекутта 4 порядка
                for (int x = 0; x < N + 2; x++)
                    for (int y = 0; y < M + 2; y++)
                    {
                        f1[x, y] = W[x, y];
                        f2[x, y] = r[x, y];
                    }
                //Console.WriteLine("f1");
                //ShowMassiv2(f1);
                //Console.WriteLine("f2");
                //ShowMassiv2(f2);
                for (int x = 2; x < N; x++)
                    for (int y = 2; y < M; y++)
                    {
                        k1[x, y] = FunctionG1(F, x, y, f1, f2, t0);
                        l1[x, y] = FunctionG2(F, x, y, f1, f2, t0);
                    }
                //Console.WriteLine("k1");
                //ShowMassiv2(k1);
                //Console.WriteLine("l1");
                //ShowMassiv2(l1);
                for (int x = 0; x < N + 2; x++)
                    for (int y = 0; y < M + 2; y++)
                    {
                        f1[x, y] = W[x, y] + Deltat / 2 * k1[x, y];
                        f2[x, y] = r[x, y] + Deltat / 2 * l1[x, y];
                    }
                //пересчёт граничных условий
                BorderZadelkSharnForMethods(f1, TypeBorder1, TypeBorder2, TypeBorder3, TypeBorder4);
                //Console.WriteLine("f1");
                //ShowMassiv2(f1);
                //Console.WriteLine("f2");
                //ShowMassiv2(f2);

                for (int x = 2; x < N; x++)
                    for (int y = 2; y < M; y++)
                    {
                        k2[x, y] = FunctionG1(F, x, y, f1, f2, t0 + Deltat / 2);
                        l2[x, y] = FunctionG2(F, x, y, f1, f2, t0 + Deltat / 2);
                    }
                //Console.WriteLine("k2");
                //ShowMassiv2(k2);
                //Console.WriteLine("l2");
                //ShowMassiv2(l2);
                for (int x = 0; x < N + 2; x++)
                    for (int y = 0; y < M + 2; y++)
                    {
                        f1[x, y] = W[x, y] + Deltat / 2 * k2[x, y];
                        f2[x, y] = r[x, y] + Deltat / 2 * l2[x, y];
                    }
                //пересчёт граничных условий
                BorderZadelkSharnForMethods(f1, TypeBorder1, TypeBorder2, TypeBorder3, TypeBorder4);
                //Console.WriteLine("f1");
                //ShowMassiv2(f1);
                //Console.WriteLine("f2");
                //ShowMassiv2(f2);
                for (int x = 2; x < N; x++)
                    for (int y = 2; y < M; y++)
                    {
                        k3[x, y] = FunctionG1(F, x, y, f1, f2, t0 + Deltat / 2);
                        l3[x, y] = FunctionG2(F, x, y, f1, f2, t0 + Deltat / 2);
                    }
                //Console.WriteLine("k3");
                //ShowMassiv2(k3);
                //Console.WriteLine("l3");
                //ShowMassiv2(l3);
                for (int x = 0; x < N + 2; x++)
                    for (int y = 0; y < M + 2; y++)
                    {
                        f1[x, y] = W[x, y] + Deltat * k3[x, y];
                        f2[x, y] = r[x, y] + Deltat * l3[x, y];
                    }
                //Console.WriteLine("f1");
                //ShowMassiv2(f1);
                //Console.WriteLine("f2");
                //ShowMassiv2(f2);
                //пересчёт граничных условий
                BorderZadelkSharnForMethods(f1, TypeBorder1, TypeBorder2, TypeBorder3, TypeBorder4);
                for (int x = 2; x < N; x++)
                    for (int y = 2; y < M; y++)
                    {
                        k4[x, y] = FunctionG1(F, x, y, f1, f2, t0 + Deltat);
                        l4[x, y] = FunctionG2(F, x, y, f1, f2, t0 + Deltat);
                    }
                //Console.WriteLine("k4");
                //ShowMassiv2(k4);
                //Console.WriteLine("l4");
                //ShowMassiv2(l4);
                for (int x = 0; x < N + 2; x++)
                    for (int y = 0; y < M + 2; y++)
                    {
                        W[x, y] += Deltat * (k1[x, y] + 2 * k2[x, y] + 2 * k3[x, y] + k4[x, y]) / 6;
                        r[x, y] += Deltat * (l1[x, y] + 2 * l2[x, y] + 2 * l3[x, y] + l4[x, y]) / 6;
                    }
                BorderZadelkSharnForMethods(W, TypeBorder1, TypeBorder2, TypeBorder3, TypeBorder4);
                //Console.WriteLine("W");
                //ShowMassiv2(W);
                //Console.WriteLine("r");
                //ShowMassiv2(r);
            }
            ///..........................................................


            ///Динамический метод с использованием метода вариационных итераций.......

            private double[,] W0;
            private double[,] W1;
            private double[,] W2;
            private Load QDT;
            public void DinamiDecisionMethodVariationIteration(Load F, double t0, double t1, double Deltat, int BorderExitPhisicalIteration,
                int TypeBorder1, int TypeBorder2, int TypeBorder3, int TypeBorder4, int WriteShag)
            {
                BorderExitIterationMethods = 100;
                // проверка на то что вызов конструктора был с учётом времени
                if (!FlagTimeProblem)
                {
                    Console.WriteLine("Error - FlagTimeProblem = false");
                    return;
                }

                RELoad();

                this.BorderExitPhisicalIteration = BorderExitPhisicalIteration;
                this.BorderExitMethodABS = 0;
                this.Deltat = Deltat;
                this.t0 = t0;
                this.t1 = t1;

                //размерность массива для заполнение сигнала
                SignalMaxW = new double[(int)((t1 - t0) / (Deltat * WriteShag)) + 1];
                //реализация массивов используемых в методе рунгекутта 4
                W0 = new double[N + 2, M + 2];
                W1 = new double[N + 2, M + 2];
                W2 = new double[N + 2, M + 2];
                QDT = new Load(N, M);
                double MaximumDeltaW = 0;


                for (int x = 0; x < N + 2; x++)
                    for (int y = 0; y < M + 2; y++)
                        W0[x, y] = W[x, y];

                for (int x = 1; x < N + 1; x++)
                    for (int y = 1; y < M + 1; y++)
                        W1[x, y] = F.F[x - 1, y - 1] / 2;




                int i = 0;
                int WriteSignal = 0;
                while (true)
                {

                    while (true)
                    {
                        BorderZadelkSharnForMethods(W1, TypeBorder1, TypeBorder2, TypeBorder3, TypeBorder4);

                        for (int x = 2; x < N; x++)
                            for (int y = 2; y < M; y++)
                            {
                                W2[x, y] = (Dissipation * W0[x, y] / (2 * Deltat) + (2 * W1[x, y] - W0[x, y]) / (Deltat * Deltat) + OperatorL(F, x, y, W1, t0)) / (1 / (Deltat * Deltat) + Dissipation / (2 * Deltat));
                            }

                        for (int x = 0; x < N; x++)
                            for (int y = 0; y < M; y++)
                            {
                                QDT.F[x, y] = F.F[x, y] - (W0[x + 1, y + 1] - 2 * W1[x + 1, y + 1] + W2[x + 1, y + 1]) / (Deltat * Deltat) - Dissipation * (W2[x + 1, y + 1] - W0[x + 1, y + 1]) / (2 * Deltat);
                            }

                        MethodVariationIteration(QDT, TypeBorder1, TypeBorder2, TypeBorder3, TypeBorder4);

                        for (int x = 0; x < N + 2; x++)
                            for (int y = 0; y < M + 2; y++)
                            {
                                if (MaximumDeltaW < Math.Abs(W1[x, y] - W[x, y])) MaximumDeltaW = Math.Abs(W1[x, y] - W[x, y]);
                                W1[x, y] = W[x, y];
                            }
                        if (MaximumDeltaW < 0.000001)
                        {
                            MaximumDeltaW = 0;
                            break;
                        }
                    }
                    for (int x = 0; x < N + 2; x++)
                        for (int y = 0; y < M + 2; y++)
                            W0[x, y] = W1[x, y];

                    for (int x = 1; x < N + 1; x++)
                        for (int y = 1; y < M + 1; y++)
                            W1[x, y] = W2[x, y];



                    this.t0 += Deltat;

                    if (this.t0 > this.t1) break;
                    //Console.WriteLine("{0:0.0000} {1:0.00000000}", this.t0, W[W.GetLength(0) / 2 + 1, W.GetLength(1) / 2 + 1]);

                    if (WriteSignal == WriteShag)
                    {
                        Console.WriteLine("{0:0.0000} {1:0.00000000}", this.t0, W[W.GetLength(0) / 2 + 1, W.GetLength(1) / 2 + 1]);
                        SignalMaxW[i / WriteShag] = W[W.GetLength(0) / 2 + 1, W.GetLength(1) / 2 + 1]; // центр пластинки
                        WriteSignal = 0;
                    }
                    i++;
                    WriteSignal++;
                    //Пересчёт физических параметров
                    Load();
                }
                //Заполнения значиний прогибов для печати
                for (int x = 0; x < N; x++)
                    for (int y = 0; y < M; y++)
                        PrintW[x, y] = W[x + 1, y + 1];
            }
            ///..................................................................


            ///Расчёт  пластинки методом вариационных итераций............................
            public void StaticDecisionMethodVariationIteration(Load F, int BorderExitIterationMethods,
                int BorderExitMethodABS, int BorderExitPhisicalIteration, int TypeBorder1, int TypeBorder2, int TypeBorder3, int TypeBorder4)
            {
                RELoad();
                this.BorderExitIterationMethods = BorderExitIterationMethods;
                this.BorderExitMethodABS = BorderExitMethodABS;
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
            public void StaticDecisionMethodKantorovichVlasovOneLevel(Load F, int BorderExitIterationMethods,
                int BorderExitMethodABS, int BorderExitPhisicalIteration, int TypeBorder1, int TypeBorder2, int TypeBorder3, int TypeBorder4)
            {
                RELoad();
                this.BorderExitIterationMethods = 0;
                this.BorderExitMethodABS = BorderExitMethodABS;
                this.BorderExitPhisicalIteration = BorderExitPhisicalIteration;
                MaxCountIterationSpendOfMetods = 0;
                MaximumCountIterationOfMethodAllDecision = 0;
                if (BorderExitMethodABS == 0)
                    PhisicalNolenearyPart(F, TypeBorder1, TypeBorder2, TypeBorder3, TypeBorder4, false, MethodKantotrovichVlasovOneLevel);
                else
                    PhisicalNolenearyPart(F, TypeBorder1, TypeBorder2, TypeBorder3, TypeBorder4, true, MethodKantotrovichVlasovOneLevel);
                //Заполнения значиний прогибов для печати
                for (int i = 0; i < N; i++)
                    for (int j = 0; j < M; j++)
                        PrintW[i, j] = W[i + 1, j + 1];

            }

            //Параметры использованнае в методе Вариационных итераций
            // Интегальные функции играющие роль в виде не постоянных коэффициентов связанных 
            //с решение линейного уравнения  4 порядка
            //p1(x)A''''(x)+p2(x)A'''(x)+p3(x)A''(x)+p4(x)A'(x)+p5(x)A(x)=p6(x)
            private double[] p1;
            private double[] p2;
            private double[] p3;
            private double[] p4;
            private double[] p5;
            private double[] p6;
            double[] A;// по х
            double[] B;// по y


            //Делегат функция принадлежадщфеащие методу 
            delegate double Functions(double[] AB, int i, int j, bool xy, double[,] F);


            // Интегральные функции при линейном дифференциальном уравнени 
            // Переменные для производных
            private double ABd2x, ABdx, ABd2y, ABdy, ABd4x, ABd3x, ABd4y, ABd3y;
            private double F1(double[] AB, int i, int j, bool xy, double[,] F)
            {
                if (xy) return (double)(Math.Pow(AB[i + 1], 2.0) * DM[i, j] * Math.Pow(lam, 4.0));
                else return (double)Math.Pow(AB[j + 1], 2.0) * DM[i, j];
            }
            private double F2(double[] AB, int i, int j, bool xy, double[,] F)
            {
                Ddx = (DM[i + 1, j] - DM[i - 1, j]) / (2 * dx);
                Ddy = (DM[i, j + 1] - DM[i, j - 1]) / (2 * dy);
                if (xy) return 2 * (double)(Math.Pow(AB[i + 1], 2.0) * Ddy * Math.Pow(lam, 4.0));
                else return 2 * (double)Math.Pow(AB[j + 1], 2.0) * Ddx;
            }
            private double F3(double[] AB, int i, int j, bool xy, double[,] F)
            {
                Dd2x = (DM[i + 1, j] - 2 * DM[i, j] + DM[i - 1, j]) / (dx * dx);
                Dd2y = (DM[i, j + 1] - 2 * DM[i, j] + DM[i, j - 1]) / (dy * dy);
                Rd2x = (RM[i + 1, j] - 2 * RM[i, j] + RM[i - 1, j]) / (dx * dx);
                Rd2y = (RM[i, j + 1] - 2 * RM[i, j] + RM[i, j - 1]) / (dy * dy);
                Rdx = (RM[i + 1, j] - RM[i - 1, j]) / (2 * dx);
                Rdy = (RM[i, j + 1] - RM[i, j - 1]) / (2 * dy);
                Sdx = (SM[i + 1, j] - SM[i - 1, j]) / (2 * dx);
                Sdy = (SM[i, j + 1] - SM[i, j - 1]) / (2 * dy);
                if (xy)
                {
                    ABd2x = (AB[i + 2] - 2 * AB[i + 1] + AB[i]) / (dx * dx);
                    ABdx = (AB[i + 2] - AB[i]) / (2 * dx);
                    return AB[i + 1] * ((double)Math.Pow(lam, 4.0) * AB[i + 1] * (Dd2y + Rd2x / (lam * lam)) +
                        2 * lam * lam * (ABd2x * (RM[i, j] + SM[i, j]) + ABdx * (Rdx + Sdx)));
                }
                else
                {
                    ABd2y = (AB[j + 2] - 2 * AB[j + 1] + AB[j]) / (dy * dy);
                    ABdy = (AB[j + 2] - AB[j]) / (2 * dy);
                    return AB[j + 1] * (AB[j + 1] * (Dd2x + lam * lam * Rd2y) +
                        2 * lam * lam * (ABd2y * (RM[i, j] + SM[i, j]) + ABdy * (Rdy + Sdy)));
                }
            }
            private double F4(double[] AB, int i, int j, bool xy, double[,] F)
            {
                Rdx = (RM[i + 1, j] - RM[i - 1, j]) / (2 * dx);
                Rdy = (RM[i, j + 1] - RM[i, j - 1]) / (2 * dy);
                Sdx = (SM[i + 1, j] - SM[i - 1, j]) / (2 * dx);
                Sdy = (SM[i, j + 1] - SM[i, j - 1]) / (2 * dy);
                S2dxdy = (SM[i + 1, j + 1] - SM[i + 1, j - 1] - SM[i - 1, j + 1] + SM[i - 1, j - 1]) / (4 * dx * dy);
                if (xy)
                {
                    ABd2x = (AB[i + 2] - 2 * AB[i + 1] + AB[i]) / (dx * dx);
                    ABdx = (AB[i + 2] - AB[i]) / (2 * dx);
                    return 2 * lam * lam * AB[i + 1] * (ABd2x * (Rdy + Sdy) + ABdx * S2dxdy);
                }
                else
                {
                    ABd2y = (AB[j + 2] - 2 * AB[j + 1] + AB[j]) / (dy * dy);
                    ABdy = (AB[j + 2] - AB[j]) / (2 * dy);
                    return 2 * lam * lam * AB[j + 1] * (ABd2y * (Rdx + Sdx) + ABdy * S2dxdy);
                }
            }
            private double F5(double[] AB, int i, int j, bool xy, double[,] F)
            {
                Dd2x = (DM[i + 1, j] - 2 * DM[i, j] + DM[i - 1, j]) / (dx * dx);
                Dd2y = (DM[i, j + 1] - 2 * DM[i, j] + DM[i, j - 1]) / (dy * dy);
                Ddx = (DM[i + 1, j] - DM[i - 1, j]) / (2 * dx);
                Ddy = (DM[i, j + 1] - DM[i, j - 1]) / (2 * dy);
                Rd2x = (RM[i + 1, j] - 2 * RM[i, j] + RM[i - 1, j]) / (dx * dx);
                Rd2y = (RM[i, j + 1] - 2 * RM[i, j] + RM[i, j - 1]) / (dy * dy);

                if (xy)
                {
                    ABd4x = (AB[i + 3] - 4 * AB[i + 2] + 6 * AB[i + 1] - 4 * AB[i] + AB[i - 1]) / (dx * dx * dx * dx);
                    ABd3x = (AB[i + 3] - 2 * AB[i + 2] + 2 * AB[i] - AB[i - 1]) / (2 * dx * dx * dx);
                    ABd2x = (AB[i + 2] - 2 * AB[i + 1] + AB[i]) / (dx * dx);
                    return AB[i + 1] * (ABd4x * DM[i, j] + 2 * ABd3x * Ddx + ABd2x * (Dd2x + Rd2y * lam * lam));
                }
                else
                {
                    ABd4y = (AB[j + 3] - 4 * AB[j + 2] + 6 * AB[j + 1] - 4 * AB[j] + AB[j - 1]) / (dy * dy * dy * dy);
                    ABd3y = (AB[j + 3] - 2 * AB[j + 2] + 2 * AB[j] - AB[j - 1]) / (2 * dy * dy * dy);
                    ABd2y = (AB[j + 2] - 2 * AB[j + 1] + AB[j]) / (dy * dy);
                    return (double)Math.Pow(lam, 4.0) * AB[j + 1] * (ABd4y * DM[i, j] + 2 * ABd3y * Ddy +
                        ABd2y * (Dd2y + Rd2x / (lam * lam)));
                }
            }
            private double F6(double[] AB, int i, int j, bool xy, double[,] F)
            {
                if (alf != 0)
                {
                    Id2x = (IM[i + 1, j] - 2 * IM[i, j] + IM[i - 1, j]) / (dx * dx);
                    Id2y = (IM[i, j + 1] - 2 * IM[i, j] + IM[i, j - 1]) / (dy * dy);
                }
                if (xy)
                {
                    return -AB[i + 1] * (Id2x + 2 * lam * lam * Id2y - F[i - 1, j - 1]);
                }
                else
                {
                    return -AB[j + 1] * (Id2x + 2 * lam * lam * Id2y - F[i - 1, j - 1]);
                }
            }
            private void Pn(double[] p, double[] AB, bool xy, Functions Fn, double[,] F)
            {
                // расчёт функций стоящих при роизводных ЛДУ, 
                //в момент нахождения функциии по одной переменной

                // Размерность Массивов N и M
                int iLength = N;
                int jLength = M;
                double f1, f2;
                // Зависит от того по какой переменно идёт находжения функций Pn
                if (xy) { iLength = M; jLength = N; }

                for (int i = 0; i < iLength - 2; i++)
                {
                    // обнуление массива
                    p[i] = 0;
                    for (int j = 1; j < jLength; j++)
                    {
                        //Интегралы считаются на основе формулы трапеций
                        if (xy)
                        {
                            // AB-это ондна из функций A(x) или B(y)
                            // Вариация относительно B(y)
                            // интегрование по х
                            f1 = Fn(AB, j, i + 2, xy, F);
                            f2 = Fn(AB, j + 1, i + 2, xy, F);
                            p[i] = p[i] + (f1 + f2) * dx / 2;
                        }
                        else
                        {
                            // Вариация относительно A(x) 
                            // интегрование по у
                            f1 = Fn(AB, i + 2, j, xy, F);
                            f2 = Fn(AB, i + 2, j + 1, xy, F);
                            p[i] = p[i] + (f1 + f2) * dy / 2;
                        }
                    }
                }
            }
            private double[] SolutionProblemMetodGaussaForVariationIteration(double d, double[] p1, double[] p2, double[] p3, double[] p4, double[] p5, double[] f, int GrUsl1, int GrUsl2)
            {
                //Грачиные условия типа известных значений функции с двух концов 
                //p1(x)y''''(x)+p2(x)y'''(x)+p3(x)y''(x)+p4(x)y'(x)+p5(x)A(x)=f(x)
                //GrUsl=0 - (по умолчанию )жёсткая заделка на граница первая производная равна нулю
                //GrUsl=1 - шарнирное операние на границе вторая производная равна нулю
                //d - шаг разбиения
                //Размерность решения зависит от подаваемых функций. если размерно P функций N, то выходной будет N+4, 
                //так как 4 производная и функции P на границах не учитываются

                //Определяем размерность
                int Ra = p1.Length;
                //задаём массив решение 
                //0 1(г) 2 .... Ra+1 Ra+2(г) Ra+3 Cчёт от (2) до (Ra+1) 
                double[] A = new double[Ra + 6];
                // матрица для метода Гаусса
                double[,] Gauss = new double[Ra, Ra + 1];

                //Заполнение матрицы Гаусса
                //Создаётся 5 диагональная матрица
                for (int i = 0; i < Ra; i++)
                {
                    if (i > 1 && i < Ra - 2)
                    {
                        Gauss[i, i + 2] = p1[i] / (d * d * d * d) + p2[i] / (2 * d * d * d);
                        Gauss[i, i + 1] = p3[i] / (d * d) + p4[i] / (2 * d) - 4 * p1[i] / (d * d * d * d)
                            - p2[i] / (d * d * d);
                        Gauss[i, i] = 6 * p1[i] / (d * d * d * d) - 2 * p3[i] / (d * d) + p5[i];
                        Gauss[i, i - 1] = p3[i] / (d * d) - p4[i] / (2 * d) - 4 * p1[i] / (d * d * d * d)
                            + p2[i] / (d * d * d);
                        Gauss[i, i - 2] = p1[i] / (d * d * d * d) - p2[i] / (2 * d * d * d);
                    }
                    Gauss[i, Ra] = f[i];
                    if (i == 0)
                    {
                        if (GrUsl1 == 0)
                        {
                            Gauss[i, i + 2] = p1[i] / (d * d * d * d) + p2[i] / (2 * d * d * d);
                            Gauss[i, i + 1] = p3[i] / (d * d) + p4[i] / (2 * d) - 4 * p1[i] / (d * d * d * d)
                                - p2[i] / (d * d * d);
                            Gauss[i, i] = 6 * p1[i] / (d * d * d * d) - 2 * p3[i] / (d * d) + p5[i]
                                + p1[i] / (d * d * d * d) - p2[i] / (2 * d * d * d);
                        }
                        if (GrUsl1 != 0)
                        {
                            Gauss[i, i + 2] = p1[i] / (d * d * d * d) + p2[i] / (2 * d * d * d);
                            Gauss[i, i + 1] = p3[i] / (d * d) + p4[i] / (2 * d) - 4 * p1[i] / (d * d * d * d)
                                - p2[i] / (d * d * d);
                            Gauss[i, i] = 6 * p1[i] / (d * d * d * d) - 2 * p3[i] / (d * d) + p5[i]
                                - p1[i] / (d * d * d * d) + p2[i] / (2 * d * d * d);
                        }
                    }
                    if (i == 1)
                    {
                        Gauss[i, i + 2] = p1[i] / (d * d * d * d) + p2[i] / (2 * d * d * d);
                        Gauss[i, i + 1] = p3[i] / (d * d) + p4[i] / (2 * d) - 4 * p1[i] / (d * d * d * d)
                            - p2[i] / (d * d * d);
                        Gauss[i, i] = 6 * p1[i] / (d * d * d * d) - 2 * p3[i] / (d * d) + p5[i];
                        Gauss[i, i - 1] = p3[i] / (d * d) - p4[i] / (2 * d) - 4 * p1[i] / (d * d * d * d)
                            + p2[i] / (d * d * d);
                    }
                    if (i == Ra - 2)
                    {
                        Gauss[i, i + 1] = p3[i] / (d * d) + p4[i] / (2 * d) - 4 * p1[i] / (d * d * d * d)
                            - p2[i] / (d * d * d);
                        Gauss[i, i] = 6 * p1[i] / (d * d * d * d) - 2 * p3[i] / (d * d) + p5[i];
                        Gauss[i, i - 1] = p3[i] / (d * d) - p4[i] / (2 * d) - 4 * p1[i] / (d * d * d * d)
                            + p2[i] / (d * d * d);
                        Gauss[i, i - 2] = p1[i] / (d * d * d * d) - p2[i] / (2 * d * d * d);
                    }
                    if (i == Ra - 1)
                    {
                        if (GrUsl2 == 0)
                        {
                            Gauss[i, i] = 6 * p1[i] / (d * d * d * d) - 2 * p3[i] / (d * d) + p5[i] + p1[i] / (d * d * d * d) + p2[i] / (2 * d * d * d);
                            Gauss[i, i - 1] = p3[i] / (d * d) - p4[i] / (2 * d) - 4 * p1[i] / (d * d * d * d)
                                + p2[i] / (d * d * d);
                            Gauss[i, i - 2] = p1[i] / (d * d * d * d) - p2[i] / (2 * d * d * d);
                        }
                        if (GrUsl2 != 0)
                        {
                            Gauss[i, i] = 6 * p1[i] / (d * d * d * d) - 2 * p3[i] / (d * d) + p5[i] - p1[i] / (d * d * d * d) - p2[i] / (2 * d * d * d);
                            Gauss[i, i - 1] = p3[i] / (d * d) - p4[i] / (2 * d) - 4 * p1[i] / (d * d * d * d)
                                + p2[i] / (d * d * d);
                            Gauss[i, i - 2] = p1[i] / (d * d * d * d) - p2[i] / (2 * d * d * d);
                        }
                    }
                }
                //Нахождение прогибов методом гаусса-жордана
                MethodGordanGauss(Gauss);

                //копирование прогибов из главной диагонали матрицы

                for (int i = 0; i < A.Length - 6; i++)
                    A[i + 3] = (double)(Gauss[i, Ra]);

                return A;
            }

            //Метод изменения граничных условий при счёте МВИ
            private void BorderZadelSharnkMVIW(double[] AB, double d, double l, int GrUsl1, int GrUsl2)
            {
                //Граничные условия жёсткой заделки
                //Слева 
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
                // Далее забиваются крайние значения на основе кубической интерполяции для нахождения

                // производной 4 порядка на границе
                double[,] A = new double[4, 5];
                for (int i = 0; i < 4; i++)
                    for (int j = 0; j < 5; j++)
                        A[i, j] = 0;
                // первый слой 
                A[0, 0] = 1;
                // второй слой
                A[1, 0] = 1;
                A[1, 1] = (double)Math.Pow(-d, 1.0);
                A[1, 2] = (double)Math.Pow(-d, 2.0);
                A[1, 3] = (double)Math.Pow(-d, 3.0);
                A[1, 4] = AB[1];

                for (int i = 2; i < 4; i++)
                {
                    //Заполнение массива А
                    A[i, 0] = 1;
                    A[i, 1] = (double)Math.Pow((i - 1) * d, 1.0);
                    A[i, 2] = (double)Math.Pow((i - 1) * d, 2.0);
                    A[i, 3] = (double)Math.Pow((i - 1) * d, 3.0);
                    A[i, 4] = AB[i + 1];
                }
                MethodGordanGauss(A);

                //double[,] A = new double[3, 4];


                //for (int i = 0; i < 3; i++)
                //{
                //    //Заполнение массива А
                //    A[i, 0] = (double)Math.Pow(i * d, 1.0);
                //    A[i, 1] = (double)Math.Pow(i * d, 2.0);
                //    A[i, 2] = (double)Math.Pow(i * d, 3.0);
                //    A[i, 3] = AB[i + 2];
                //    if (i == 0)
                //    {
                //        A[i, 0] = (double)Math.Pow(-d, 1.0);
                //        A[i, 1] = (double)Math.Pow(-d, 2.0);
                //        A[i, 2] = (double)Math.Pow(-d, 3.0);
                //        A[i, 3] = AB[i + 1];
                //    }
                //}
                //MethodGordanGauss(A);

                //Находим значение спомощью кубического уравнения
                //AB[0] = (double)(A[2, 3] * (double)Math.Pow(-2 * d, 3.0) + A[1, 3] * (double)Math.Pow(-2 * d, 2.0) + A[0, 3] * (-2 * d));
                AB[0] = (double)(A[3, 4] * (double)Math.Pow(-2 * d, 3.0) + A[2, 4] * (double)Math.Pow(-2 * d, 2.0) + A[1, 4] * (-2 * d) + A[0, 4]);


                for (int i = 0; i < 4; i++)
                {
                    //Заполнение массива А
                    A[i, 0] = 1;
                    A[i, 1] = (double)Math.Pow(l + (i - 2) * d, 1.0);
                    A[i, 2] = (double)Math.Pow(l + (i - 2) * d, 2.0);
                    A[i, 3] = (double)Math.Pow(l + (i - 2) * d, 3.0);
                    A[i, 4] = AB[AB.Length - 5 + i];
                }
                MethodGordanGauss(A);
                AB[AB.Length - 1] = (double)(A[3, 4] * (double)Math.Pow(l + 2 * d, 3.0) + A[2, 4]
                    * (double)Math.Pow(l + 2 * d, 2.0) + A[1, 4] * (l + 2 * d) + A[0, 4]);
                //for (int i = 0; i < 3; i++)
                //{
                //    //Заполнение массива А
                //    A[i, 0] = (double)Math.Pow(1 + (i - 2) * d, 1.0);
                //    A[i, 1] = (double)Math.Pow(1 + (i - 2) * d, 2.0);
                //    A[i, 2] = (double)Math.Pow(1 + (i - 2) * d, 3.0);
                //    A[i, 3] = AB[AB.Length - 5 + i];
                //    if (i == 2)
                //    {
                //        A[i, 0] = (double)Math.Pow(1 + d, 1.0);
                //        A[i, 1] = (double)Math.Pow(1 + d, 2.0);
                //        A[i, 2] = (double)Math.Pow(1 + d, 3.0);
                //        A[i, 3] = AB[AB.Length - 2];
                //    }
                //}
                //MethodGordanGauss(A);

                //AB[AB.Length - 1] = (double)(A[2, 3] * (double)Math.Pow(2 * d + 1, 3.0) + A[1, 3]
                //    * (double)Math.Pow(2 * d + 1, 2.0) + A[0, 3] * (2 * d + 1));
            }
            private void BorderSharnMVIW(double[] AB, double d)
            {
                //Граничные условия Шарнирного опирания
                //Слева 
                AB[2] = 0;
                AB[1] = -AB[3];
                //Справа
                AB[AB.Length - 3] = 0;
                AB[AB.Length - 2] = -AB[AB.Length - 4];

                double[,] A = new double[4, 5];
                for (int i = 0; i < 4; i++)
                {
                    //Заполнение массива А
                    A[i, 0] = 1;
                    A[i, 1] = (double)Math.Pow((i - 1) * d, 1.0);
                    A[i, 2] = (double)Math.Pow((i - 1) * d, 2.0);
                    A[i, 3] = (double)Math.Pow((i - 1) * d, 3.0);
                    A[i, 4] = AB[i + 1];
                }
                MethodGordanGauss(A);
                //Находим значение спомощью кубического уравнения
                AB[0] = (double)(A[3, 4] * (double)Math.Pow(-2 * d, 3.0) + A[2, 4] * (double)Math.Pow(-2 * d, 2.0) + A[1, 4] * (-2 * d) + A[0, 4]);

                for (int i = 0; i < 4; i++)
                {
                    //Заполнение массива А
                    A[i, 0] = 1;
                    A[i, 1] = (double)Math.Pow(1 - (i - 1) * d, 1.0);
                    A[i, 2] = (double)Math.Pow(1 - (i - 1) * d, 2.0);
                    A[i, 3] = (double)Math.Pow(1 - (i - 1) * d, 3.0);
                    A[i, 4] = AB[AB.Length - 2 - i];
                }
                MethodGordanGauss(A);
                AB[AB.Length - 1] = (double)(A[3, 4] * (double)Math.Pow(d * (AB.Length - 1), 3.0) + A[2, 4]
                    * (double)Math.Pow(d * (AB.Length - 1), 2.0) + A[1, 4] * d * (AB.Length - 1) + A[0, 4]);
            }
            private void MethodKantotrovichVlasovOneLevel(Load F, int TypeBorder1, int TypeBorder2, int TypeBorder3, int TypeBorder4)
            {
                mW = 0;
                A = new double[N + 4];// по x
                B = new double[M + 4];// по y

                // переменная для счёта
                double Y = 0;
                //Заполнение функций в зависимости от граничных условий
                if (TypeBorder1 == 1 && TypeBorder3 == 1)
                {
                    for (int y = 0; y < M + 4; y++)
                    {
                        Y = dy * (y - 2);
                        B[y] = Math.Sin(Y * Math.PI);
                    }
                }
                if (TypeBorder1 == 1 && TypeBorder3 == 0)
                {
                    for (int y = 0; y < M + 4; y++)
                    {
                        Y = dy * (y - 2);
                        B[y] = Math.Sin(Y * Math.PI / 2) * Math.Sin(Y * Math.PI);
                    }
                }
                if (TypeBorder1 == 0 && TypeBorder3 == 1)
                {
                    for (int y = 0; y < M + 4; y++)
                    {
                        Y = dy * (y - 2);
                        B[y] = Math.Cos(Y * Math.PI / 2) * Math.Sin(Y * Math.PI);
                    }
                }
                if (TypeBorder1 == 0 && TypeBorder3 == 0)
                {
                    for (int y = 0; y < M + 4; y++)
                    {
                        Y = dy * (y - 2);
                        B[y] = Math.Sin(Y * Math.PI) * Math.Sin(Y * Math.PI);
                    }
                }
                else
                {
                    for (int y = 0; y < M + 4; y++)
                    {
                        Y = dy * (y - 2);
                        B[y] = Math.Sin(Y * Math.PI);
                    }
                }

                // задаём размерность массивам по переменной х
                p1 = new double[N - 2];
                p2 = new double[N - 2];
                p3 = new double[N - 2];
                p4 = new double[N - 2];
                p5 = new double[N - 2];
                p6 = new double[N - 2];

                // находим переменные функции стоящие при производных в ДУ.
                Pn(p1, B, false, F1, F.F);
                //ShowMassiv(p1);
                Pn(p2, B, false, F2, F.F);
                //ShowMassiv(p2);
                Pn(p3, B, false, F3, F.F);
                //ShowMassiv(p3);
                Pn(p4, B, false, F4, F.F);
                //ShowMassiv(p4);
                Pn(p5, B, false, F5, F.F);
                //ShowMassiv(p5);
                Pn(p6, B, false, F6, F.F);
                //ShowMassiv(p6);
                // находим решение уравнения методом гаусса жордана
                A = SolutionProblemMetodGaussaForVariationIteration(dx, p1, p2, p3, p4, p5, p6, TypeBorder4, TypeBorder2);
                BorderZadelSharnkMVIW(A, dx, n, TypeBorder4, TypeBorder2);
                //Заполнение массива определяющий прогиб и Нахожденеие максимального прогиба
                for (int i = 0; i < N + 2; i++)
                    for (int j = 0; j < M + 2; j++)
                    {
                        W[i, j] = A[i + 1] * B[j + 1];
                        if (Math.Abs(W[i, j]) > mW) mW = Math.Abs(W[i, j]); // Нахождение максимального прогиба
                    }

            }

            // метод вариационных итераций 
            private void MethodVariationIteration(Load F, int TypeBorder1, int TypeBorder2, int TypeBorder3, int TypeBorder4)
            {
                //Основные 2 Функции по каздой переменной
                //кол.т    0   1 .....  N-2  N-1(г)
                //A = 0 1 2(г) 3 ....    N  N+1(г) N+2 N+3
                //p = - - -(г) 0 1 2..  N-3  -(г)   -   -
                //DM = - 0 1(г) 2 ....   N-1  N(г)  N+1  - 
                A = new double[N + 4];// по x
                B = new double[M + 4];// по y
                // Задаём начальную функцию по y 
                for (int y = 0; y < M + 4; y++)
                    B[y] = (double)(Math.Sin(dy * (y - 2)));

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

                    // задаём размерность массивам по переменной х
                    p1 = new double[N - 2];
                    p2 = new double[N - 2];
                    p3 = new double[N - 2];
                    p4 = new double[N - 2];
                    p5 = new double[N - 2];
                    p6 = new double[N - 2];

                    // находим переменные функции стоящие при производных в ДУ.
                    Pn(p1, B, false, F1, F.F);
                    //ShowMassiv(p1);
                    Pn(p2, B, false, F2, F.F);
                    //ShowMassiv(p2);
                    Pn(p3, B, false, F3, F.F);
                    //ShowMassiv(p3);
                    Pn(p4, B, false, F4, F.F);
                    //ShowMassiv(p4);
                    Pn(p5, B, false, F5, F.F);
                    //ShowMassiv(p5);
                    Pn(p6, B, false, F6, F.F);
                    //ShowMassiv(p6);
                    // находим решение уравнения методом гаусса жордана
                    A = SolutionProblemMetodGaussaForVariationIteration(dx, p1, p2, p3, p4, p5, p6, TypeBorder4, TypeBorder2);
                    // относительно граничного условия и кубической интерполяцией задаём заграничные значения.
                    //граничные условия 
                    BorderZadelSharnkMVIW(A, dx, n, TypeBorder4, TypeBorder2);


                    //Переопределение массивов
                    p1 = new double[M - 2];
                    p2 = new double[M - 2];
                    p3 = new double[M - 2];
                    p4 = new double[M - 2];
                    p5 = new double[M - 2];
                    p6 = new double[M - 2];


                    Pn(p1, A, true, F1, F.F);
                    //ShowMassiv(p1);
                    Pn(p2, A, true, F2, F.F);
                    //ShowMassiv(p2);
                    Pn(p3, A, true, F3, F.F);
                    //ShowMassiv(p3);
                    Pn(p4, A, true, F4, F.F);
                    //ShowMassiv(p4);
                    Pn(p5, A, true, F5, F.F);
                    //ShowMassiv(p5);
                    Pn(p6, A, true, F6, F.F);
                    //ShowMassiv(p6);

                    B = SolutionProblemMetodGaussaForVariationIteration(dy, p1, p2, p3, p4, p5, p6, TypeBorder3, TypeBorder1);
                    //граничные условия
                    BorderZadelSharnkMVIW(B, dy, m, TypeBorder3, TypeBorder1);


                    CountIterationInMethods++;
                    for (int i = 0; i < N + 2; i++)
                        for (int j = 0; j < M + 2; j++)
                        {
                            if (Math.Abs(W[i, j] - A[i + 1] * B[j + 1]) > MaxWForMethodIteration1) MaxWForMethodIteration1 = Math.Abs(W[i, j] - A[i + 1] * B[j + 1]);// нахождение максмального отклонения найденного решения.
                            W[i, j] = A[i + 1] * B[j + 1];
                            if (Math.Abs(W[i, j]) > MaxWForMethodIteration2) MaxWForMethodIteration2 = Math.Abs(W[i, j]); // Нахождение максимального прогиба
                        }
                    if (Math.Abs(MaxWForMethodIteration1) < 0.0000001) break;
                }
                if (MaximumCountIterationOfMethodAllDecision < CountIterationInMethods) MaximumCountIterationOfMethodAllDecision = CountIterationInMethods;
                MaxCountIterationSpendOfMetods += CountIterationInMethods;
                mW = MaxWForMethodIteration2;
            }
            /// ................................................................


            ///Метод Канторовича - Власова....................................................................
            public void StaticDecisionMethodVariationIterartionTwoLevel(Load F, int BorderExitIterationMethods,
                int BorderExitMethodABS, int BorderExitPhisicalIteration, int TypeBorder1, int TypeBorder2, int TypeBorder3, int TypeBorder4)
            {
                RELoad();
                this.BorderExitIterationMethods = BorderExitIterationMethods;
                this.BorderExitMethodABS = BorderExitMethodABS;
                this.BorderExitPhisicalIteration = BorderExitPhisicalIteration;
                this.MaxCountIterationSpendOfMetods = 0;
                if (BorderExitMethodABS == 0)
                    PhisicalNolenearyPart(F, TypeBorder1, TypeBorder2, TypeBorder3, TypeBorder4, false, MethodVariationIterartionTwoLevel);
                else
                    PhisicalNolenearyPart(F, TypeBorder1, TypeBorder2, TypeBorder3, TypeBorder4, true, MethodVariationIterartionTwoLevel);
                //Заполнения значиний прогибов для печати
                for (int i = 0; i < N; i++)
                    for (int j = 0; j < M; j++)
                        PrintW[i, j] = W[i + 1, j + 1];

            }
            public void StaticDecisionMethodKantorovichVlasovTwoLevel(Load F, int BorderExitIterationMethods,
                int BorderExitMethodABS, int BorderExitPhisicalIteration, int TypeBorder1, int TypeBorder2, int TypeBorder3, int TypeBorder4)
            {
                RELoad();
                this.BorderExitIterationMethods = 0;
                this.BorderExitMethodABS = BorderExitMethodABS;
                this.BorderExitPhisicalIteration = BorderExitPhisicalIteration;
                this.MaxCountIterationSpendOfMetods = 0;
                if (BorderExitMethodABS == 0)
                    PhisicalNolenearyPart(F, TypeBorder1, TypeBorder2, TypeBorder3, TypeBorder4, false, MethodKantorovichVlasovTwoLevel);
                else
                    PhisicalNolenearyPart(F, TypeBorder1, TypeBorder2, TypeBorder3, TypeBorder4, true, MethodKantorovichVlasovTwoLevel);
                //Заполнения значиний прогибов для печати
                for (int i = 0; i < N; i++)
                    for (int j = 0; j < M; j++)
                        PrintW[i, j] = W[i + 1, j + 1];

            }

            private double[] C1;
            private double[] C2;
            private double[] C3;
            private double[] C4;
            private double[] C5;
            private double[] C6;
            private double[] C7;
            private double[] C8;
            private double[] C9;
            private double[] C10;
            private double[] C11;

            private double[] H1;
            private double[] H2;
            private double[] H3;
            private double[] H4;
            private double[] H5;
            private double[] H6;
            private double[] H7;
            private double[] H8;
            private double[] H9;
            private double[] H10;
            private double[] H11;

            double[] A1;// по х
            double[] B1;// по y
            double[] A2;// по х
            double[] B2;// по y
            // для сравнения с предыдущим шагом итерации для решения системы диф уравнений
            double[] AB1;// по х
            double[] AB2;// по y


            private double AnBnd2x, AnBndx, AnBnd2y, AnBndy, AnBnd4x, AnBnd3x, AnBnd4y, AnBnd3y;

            //Делегат функция принадлежадщфеащие методу 
            delegate double FunctionsA1A2B1B2(double[] A1B1, double[] A2B2, int i, int j, bool AB, double[,] F);
            private double L1(double[] A1B1, double[] A2B2, int i, int j, bool AB, double[,] F)
            {
                if (AB) return (double)(A1B1[i + 1] * A2B2[i + 1] * DM[i, j] * Math.Pow(lam, 4.0));
                else return (double)(A1B1[j + 1] * A2B2[j + 1] * DM[i, j]);
            }
            private double L2(double[] A1B1, double[] A2B2, int i, int j, bool AB, double[,] F)
            {
                Ddx = (DM[i + 1, j] - DM[i - 1, j]) / (2 * dx);
                Ddy = (DM[i, j + 1] - DM[i, j - 1]) / (2 * dy);
                if (AB) return 2 * (double)(A1B1[i + 1] * A2B2[i + 1] * Ddy * Math.Pow(lam, 4.0));
                else return 2 * (double)(A1B1[j + 1] * A2B2[j + 1] * Ddx);
            }
            private double L3(double[] A1B1, double[] A2B2, int i, int j, bool AB, double[,] F)
            {
                Dd2x = (DM[i + 1, j] - 2 * DM[i, j] + DM[i - 1, j]) / (dx * dx);
                Dd2y = (DM[i, j + 1] - 2 * DM[i, j] + DM[i, j - 1]) / (dy * dy);
                Rd2x = (RM[i + 1, j] - 2 * RM[i, j] + RM[i - 1, j]) / (dx * dx);
                Rd2y = (RM[i, j + 1] - 2 * RM[i, j] + RM[i, j - 1]) / (dy * dy);
                Rdx = (RM[i + 1, j] - RM[i - 1, j]) / (2 * dx);
                Rdy = (RM[i, j + 1] - RM[i, j - 1]) / (2 * dy);
                Sdx = (SM[i + 1, j] - SM[i - 1, j]) / (2 * dx);
                Sdy = (SM[i, j + 1] - SM[i, j - 1]) / (2 * dy);

                if (AB)
                {
                    AnBnd2x = (A2B2[i + 2] - 2 * A2B2[i + 1] + A2B2[i]) / (dx * dx);
                    AnBndx = (A2B2[i + 2] - A2B2[i]) / (2 * dx);
                    return A1B1[i + 1] * ((double)Math.Pow(lam, 4.0) * A2B2[i + 1] * (Dd2y + Rd2x / (lam * lam)) +
                        2 * lam * lam * (AnBnd2x * (RM[i, j] + SM[i, j]) + AnBndx * (Rdx + Sdx)));
                }
                else
                {
                    AnBnd2y = (A2B2[j + 2] - 2 * A2B2[j + 1] + A2B2[j]) / (dy * dy);
                    AnBndy = (A2B2[j + 2] - A2B2[j]) / (2 * dy);
                    return A1B1[j + 1] * (A2B2[j + 1] * (Dd2x + lam * lam * Rd2y) +
                        2 * lam * lam * (AnBnd2y * (RM[i, j] + SM[i, j]) + AnBndy * (Rdy + Sdy)));
                }
            }
            private double L4(double[] A1B1, double[] A2B2, int i, int j, bool AB, double[,] F)
            {
                Rdx = (RM[i + 1, j] - RM[i - 1, j]) / (2 * dx);
                Rdy = (RM[i, j + 1] - RM[i, j - 1]) / (2 * dy);
                Sdx = (SM[i + 1, j] - SM[i - 1, j]) / (2 * dx);
                Sdy = (SM[i, j + 1] - SM[i, j - 1]) / (2 * dy);
                S2dxdy = (SM[i + 1, j + 1] - SM[i + 1, j - 1] - SM[i - 1, j + 1] + SM[i - 1, j - 1]) / (4 * dx * dy);
                if (AB)
                {
                    AnBnd2x = (A2B2[i + 2] - 2 * A2B2[i + 1] + A2B2[i]) / (dx * dx);
                    AnBndx = (A2B2[i + 2] - A2B2[i]) / (2 * dx);
                    return 2 * lam * lam * A1B1[i + 1] * (AnBnd2x * (Rdy + Sdy) + AnBndx * S2dxdy);
                }
                else
                {
                    AnBnd2y = (A2B2[j + 2] - 2 * A2B2[j + 1] + A2B2[j]) / (dy * dy);
                    AnBndy = (A2B2[j + 2] - A2B2[j]) / (2 * dy);
                    return 2 * lam * lam * A1B1[j + 1] * (AnBnd2y * (Rdx + Sdx) + AnBndy * S2dxdy);
                }
            }
            private double L5(double[] A1B1, double[] A2B2, int i, int j, bool AB, double[,] F)
            {
                Dd2x = (DM[i + 1, j] - 2 * DM[i, j] + DM[i - 1, j]) / (dx * dx);
                Dd2y = (DM[i, j + 1] - 2 * DM[i, j] + DM[i, j - 1]) / (dy * dy);
                Ddx = (DM[i + 1, j] - DM[i - 1, j]) / (2 * dx);
                Ddy = (DM[i, j + 1] - DM[i, j - 1]) / (2 * dy);
                Rd2x = (RM[i + 1, j] - 2 * RM[i, j] + RM[i - 1, j]) / (dx * dx);
                Rd2y = (RM[i, j + 1] - 2 * RM[i, j] + RM[i, j - 1]) / (dy * dy);

                if (AB)
                {
                    AnBnd4x = (A2B2[i + 3] - 4 * A2B2[i + 2] + 6 * A2B2[i + 1] - 4 * A2B2[i] + A2B2[i - 1]) / (dx * dx * dx * dx);
                    AnBnd3x = (A2B2[i + 3] - 2 * A2B2[i + 2] + 2 * A2B2[i] - A2B2[i - 1]) / (2 * dx * dx * dx);
                    AnBnd2x = (A2B2[i + 2] - 2 * A2B2[i + 1] + A2B2[i]) / (dx * dx);
                    return A1B1[i + 1] * (AnBnd4x * DM[i, j] + 2 * AnBnd3x * Ddx + AnBnd2x * (Dd2x + Rd2y * lam * lam));
                }
                else
                {
                    AnBnd4y = (A2B2[j + 3] - 4 * A2B2[j + 2] + 6 * A2B2[j + 1] - 4 * A2B2[j] + A2B2[j - 1]) / (dy * dy * dy * dy);
                    AnBnd3y = (A2B2[j + 3] - 2 * A2B2[j + 2] + 2 * A2B2[j] - A2B2[j - 1]) / (2 * dy * dy * dy);
                    AnBnd2y = (A2B2[j + 2] - 2 * A2B2[j + 1] + A2B2[j]) / (dy * dy);
                    return (double)Math.Pow(lam, 4.0) * A1B1[j + 1] * (AnBnd4y * DM[i, j] + 2 * AnBnd3y * Ddy +
                        AnBnd2y * (Dd2y + Rd2x / (lam * lam)));
                }
            }
            private double L6(double[] A1B1, double[] A2B2, int i, int j, bool AB, double[,] F)
            {
                if (alf != 0)
                {
                    Id2x = (IM[i + 1, j] - 2 * IM[i, j] + IM[i - 1, j]) / (dx * dx);
                    Id2y = (IM[i, j + 1] - 2 * IM[i, j] + IM[i, j - 1]) / (dy * dy);
                }
                if (AB)
                {
                    return A1B1[i + 1] * (Id2x + 2 * lam * lam * Id2y - F[i - 1, j - 1]);
                }
                else
                {
                    return A1B1[j + 1] * (Id2x + 2 * lam * lam * Id2y - F[i - 1, j - 1]);
                }
            }

            private void IntegrationL(double[] Cn, double[] A1B1, double[] A2B2, bool AB, FunctionsA1A2B1B2 FunctionL, double[,] F)
            {
                // расчёт функций стоящих при роизводных ЛДУ, 
                //в момент нахождения функциии по одной переменной

                // Размерность Массивов N и M
                int iLength = N;
                int jLength = M;
                double f1, f2;
                // Зависит от того по какой переменно идёт находжения функций Pn
                if (AB) { iLength = M; jLength = N; }

                for (int i = 0; i < iLength - 2; i++)
                {
                    // обнуление массива
                    Cn[i] = 0;
                    for (int j = 1; j < jLength; j++)
                    {
                        //Интегралы считаются на основе формулы трапеций
                        if (AB)
                        {
                            // AB-это ондна из функций A(x) или B(y)
                            // Вариация относительно B(y)
                            // интегрование по х
                            f1 = FunctionL(A1B1, A2B2, j, i + 2, AB, F);
                            f2 = FunctionL(A1B1, A2B2, j + 1, i + 2, AB, F);
                            Cn[i] = Cn[i] + (f1 + f2) * dx / 2;
                        }
                        else
                        {
                            // Вариация относительно A(x) 
                            // интегрование по у
                            f1 = FunctionL(A1B1, A2B2, i + 2, j, AB, F);
                            f2 = FunctionL(A1B1, A2B2, i + 2, j + 1, AB, F);
                            Cn[i] = Cn[i] + (f1 + f2) * dy / 2;
                        }
                    }
                }
            }
            private double[] SolutionProblemMetodGaussaForKantorovichaVlasova(double d, double[] AnBn, double[] C1, double[] C2, double[] C3, double[] C4, double[] C5,
                double[] C6, double[] C7, double[] C8, double[] C9, double[] C10, double[] C11, int GrUsl1, int GrUsl2)
            {
                //Грачиные условия типа известных значений функции с двух концов 
                //p1(x)y''''(x)+p2(x)y'''(x)+p3(x)y''(x)+p4(x)y'(x)+p5(x)A(x)=f(x)
                //GrUsl=0 - (по умолчанию )жёсткая заделка на граница первая производная равна нулю
                //GrUsl=1 - шарнирное операние на границе вторая производная равна нулю
                //d - шаг разбиения
                //Размерность решения зависит от подаваемых функций. если размерно P функций N, то выходной будет N+4, 
                //так как 4 производная и функции P на границах не учитываются

                //Определяем размерность
                int Ra = C1.Length;
                //задаём массив решение 
                //0 1(г) 2 .... Ra+1 Ra+2(г) Ra+3 Cчёт от (2) до (Ra+1) 
                double[] A = new double[Ra + 6];
                // матрица для метода Гаусса
                double[,] Gauss = new double[Ra, Ra + 1];
                //Переменные для производных по переменной 
                double ABd1 = 0, ABd2 = 0, ABd3 = 0, ABd4 = 0;



                //Заполнение матрицы Гаусса
                //Создаётся 5 диагональная матрица
                for (int i = 0; i < Ra; i++)
                {
                    if (i > 1 && i < Ra - 2)
                    {
                        Gauss[i, i + 2] = C1[i] / (d * d * d * d) + C2[i] / (2 * d * d * d);
                        Gauss[i, i + 1] = C3[i] / (d * d) + C4[i] / (2 * d) - 4 * C1[i] / (d * d * d * d)
                            - C2[i] / (d * d * d);
                        Gauss[i, i] = 6 * C1[i] / (d * d * d * d) - 2 * C3[i] / (d * d) + C5[i];
                        Gauss[i, i - 1] = C3[i] / (d * d) - C4[i] / (2 * d) - 4 * C1[i] / (d * d * d * d)
                            + C2[i] / (d * d * d);
                        Gauss[i, i - 2] = C1[i] / (d * d * d * d) - C2[i] / (2 * d * d * d);
                    }
                    ABd1 = (AnBn[i + 1 + 3] - AnBn[i - 1 + 3]) / (2 * d);
                    ABd2 = (AnBn[i + 1 + 3] - 2 * AnBn[i + 3] + AnBn[i - 1 + 3]) / (d * d);
                    ABd3 = (AnBn[i + 2 + 3] - 2 * AnBn[i + 1 + 3] + 2 * AnBn[i - 1 + 3] - AnBn[i - 2 + 3]) / (2 * d * d * d);
                    ABd4 = (AnBn[i + 2 + 3] - 4 * AnBn[i + 1 + 3] + 6 * AnBn[i + 3] - 4 * AnBn[i - 1 + 3] + AnBn[i - 2 + 3]) / (d * d * d * d);
                    Gauss[i, Ra] = -C11[i] - C6[i] * ABd4 - C7[i] * ABd3 - C8[i] * ABd2 - C9[i] * ABd1 - C10[i] * AnBn[i + 3];
                    if (i == 0)
                    {
                        if (GrUsl1 == 0)
                        {
                            Gauss[i, i + 2] = C1[i] / (d * d * d * d) + C2[i] / (2 * d * d * d);
                            Gauss[i, i + 1] = C3[i] / (d * d) + C4[i] / (2 * d) - 4 * C1[i] / (d * d * d * d)
                                - C2[i] / (d * d * d);
                            Gauss[i, i] = 6 * C1[i] / (d * d * d * d) - 2 * C3[i] / (d * d) + C5[i]
                                + C1[i] / (d * d * d * d) - C2[i] / (2 * d * d * d);
                        }
                        if (GrUsl1 != 0)
                        {
                            Gauss[i, i + 2] = C1[i] / (d * d * d * d) + C2[i] / (2 * d * d * d);
                            Gauss[i, i + 1] = C3[i] / (d * d) + C4[i] / (2 * d) - 4 * C1[i] / (d * d * d * d)
                                - C2[i] / (d * d * d);
                            Gauss[i, i] = 6 * C1[i] / (d * d * d * d) - 2 * C3[i] / (d * d) + C5[i]
                                - C1[i] / (d * d * d * d) + C2[i] / (2 * d * d * d);
                        }
                    }
                    if (i == 1)
                    {
                        Gauss[i, i + 2] = C1[i] / (d * d * d * d) + C2[i] / (2 * d * d * d);
                        Gauss[i, i + 1] = C3[i] / (d * d) + C4[i] / (2 * d) - 4 * C1[i] / (d * d * d * d)
                            - C2[i] / (d * d * d);
                        Gauss[i, i] = 6 * C1[i] / (d * d * d * d) - 2 * C3[i] / (d * d) + C5[i];
                        Gauss[i, i - 1] = C3[i] / (d * d) - C4[i] / (2 * d) - 4 * C1[i] / (d * d * d * d)
                            + C2[i] / (d * d * d);
                    }
                    if (i == Ra - 2)
                    {
                        Gauss[i, i + 1] = C3[i] / (d * d) + C4[i] / (2 * d) - 4 * C1[i] / (d * d * d * d)
                            - C2[i] / (d * d * d);
                        Gauss[i, i] = 6 * C1[i] / (d * d * d * d) - 2 * C3[i] / (d * d) + C5[i];
                        Gauss[i, i - 1] = C3[i] / (d * d) - C4[i] / (2 * d) - 4 * C1[i] / (d * d * d * d)
                            + C2[i] / (d * d * d);
                        Gauss[i, i - 2] = C1[i] / (d * d * d * d) - C2[i] / (2 * d * d * d);
                    }
                    if (i == Ra - 1)
                    {
                        if (GrUsl2 == 0)
                        {
                            Gauss[i, i] = 6 * C1[i] / (d * d * d * d) - 2 * C3[i] / (d * d) + C5[i] + C1[i] / (d * d * d * d) + C2[i] / (2 * d * d * d);
                            Gauss[i, i - 1] = C3[i] / (d * d) - C4[i] / (2 * d) - 4 * C1[i] / (d * d * d * d)
                                + C2[i] / (d * d * d);
                            Gauss[i, i - 2] = C1[i] / (d * d * d * d) - C2[i] / (2 * d * d * d);
                        }
                        if (GrUsl2 != 0)
                        {
                            Gauss[i, i] = 6 * C1[i] / (d * d * d * d) - 2 * C3[i] / (d * d) + C5[i] - C1[i] / (d * d * d * d) - C2[i] / (2 * d * d * d);
                            Gauss[i, i - 1] = C3[i] / (d * d) - C4[i] / (2 * d) - 4 * C1[i] / (d * d * d * d)
                                + C2[i] / (d * d * d);
                            Gauss[i, i - 2] = C1[i] / (d * d * d * d) - C2[i] / (2 * d * d * d);
                        }
                    }
                }
                //Нахождение прогибов методом гаусса-жордана
                MethodGordanGauss(Gauss);

                //копирование прогибов из главной диагонали матрицы

                for (int i = 0; i < A.Length - 6; i++)
                    A[i + 3] = (double)(Gauss[i, Ra]);

                return A;
            }
            private double MaximumModuleBetweenTwoMassiveOneDimensional(double[] FirstMassive, double[] SecondMassive)
            {
                double maximum = 0;
                if (FirstMassive.Length != SecondMassive.Length) return 0;
                for (int i = 0; i < FirstMassive.Length; i++)
                    if (maximum > Math.Abs(FirstMassive[i] - SecondMassive[i]))
                        maximum = Math.Abs(FirstMassive[i] - SecondMassive[i]);
                return maximum;

            }
            private double[] ArrayAssignment(double[] Massive)
            {
                double[] newMassive = new double[Massive.Length];
                for (int i = 0; i < Massive.Length; i++)
                    newMassive[i] = Massive[i];
                return newMassive;
            }
            //Метод производный от метода вариационных итераций 
            private void MethodKantorovichVlasovTwoLevel(Load F, int TypeBorder1, int TypeBorder2, int TypeBorder3, int TypeBorder4)
            {
                //
                mW = 0;
                //Переменные для записи максимальной разницы по норме, настоящего решения с предыдущим. 
                double MaxModule1 = 0;
                double MaxModule2 = 0;

                A1 = new double[N + 4];// по x
                B1 = new double[M + 4];// по y
                A2 = new double[N + 4];// по x
                B2 = new double[M + 4];// по y

                // для сравнения с предыдущим шагом итерации для решения системы диф уравнений
                AB1 = new double[N + 4];
                AB2 = new double[N + 4];


                // задаём размерность массивам по переменной х
                C1 = new double[N - 2];
                C2 = new double[N - 2];
                C3 = new double[N - 2];
                C4 = new double[N - 2];
                C5 = new double[N - 2];
                C6 = new double[N - 2];
                C7 = new double[N - 2];
                C8 = new double[N - 2];
                C9 = new double[N - 2];
                C10 = new double[N - 2];
                C11 = new double[N - 2];
                // переменная для счёта
                double Y = 0;
                //Заполнение функций в зависимости от граничных условий
                if (TypeBorder1 == 1 && TypeBorder3 == 1)
                {
                    for (int y = 0; y < M + 4; y++)
                    {
                        Y = dy * (y - 2);
                        B1[y] = Math.Sin(Y * Math.PI);
                        B2[y] = Math.Sin(Y * 2 * Math.PI);
                    }
                }
                if (TypeBorder1 == 1 && TypeBorder3 == 0)
                {
                    for (int y = 0; y < M + 4; y++)
                    {
                        Y = dy * (y - 2);
                        B1[y] = Math.Sin(Y * Math.PI / 2) * Math.Sin(Y * Math.PI);
                        B2[y] = Math.Sin(Y * Math.PI * 3 / 2) * Math.Sin(3 * Y * Math.PI);
                    }
                }
                if (TypeBorder1 == 0 && TypeBorder3 == 1)
                {
                    for (int y = 0; y < M + 4; y++)
                    {
                        Y = dy * (y - 2);
                        B1[y] = Math.Cos(Y * Math.PI / 2) * Math.Sin(Y * Math.PI);
                        B2[y] = Math.Cos(Y * Math.PI * 3 / 2) * Math.Sin(3 * Y * Math.PI);
                    }
                }
                if (TypeBorder1 == 0 && TypeBorder3 == 0)
                {
                    for (int y = 0; y < M + 4; y++)
                    {
                        Y = dy * (y - 2);
                        B1[y] = Math.Sin(Y * Math.PI) * Math.Sin(Y * Math.PI);
                        B2[y] = Math.Sin(2 * Y * Math.PI) * Math.Sin(2 * Y * Math.PI);
                    }
                }
                else
                {
                    for (int y = 0; y < M + 4; y++)
                    {
                        Y = dy * (y - 2);
                        B1[y] = Math.Sin(Y * Math.PI);
                        B2[y] = Math.Sin(Y * 2 * Math.PI);
                    }
                }
                for (int x = 0; x < N + 4; x++)
                    A2[x] = 0;

                //Решение системы дифференциальных уравнений относительно x
                while (true)
                {
                    // находим переменные функции стоящие при производных в ДУ.
                    IntegrationL(C1, B1, B1, false, L1, F.F);
                    IntegrationL(C2, B1, B1, false, L2, F.F);
                    IntegrationL(C3, B1, B1, false, L3, F.F);
                    IntegrationL(C4, B1, B1, false, L4, F.F);
                    IntegrationL(C5, B1, B1, false, L5, F.F);
                    IntegrationL(C6, B1, B2, false, L1, F.F);
                    IntegrationL(C7, B1, B2, false, L2, F.F);
                    IntegrationL(C8, B1, B2, false, L3, F.F);
                    IntegrationL(C9, B1, B2, false, L4, F.F);
                    IntegrationL(C10, B1, B2, false, L5, F.F);
                    IntegrationL(C11, B1, B2, false, L6, F.F);

                    // находим решение уравнения методом гаусса жордана для функции А1
                    A1 = SolutionProblemMetodGaussaForKantorovichaVlasova(dx, A2, C1, C2, C3, C4, C5, C6, C7, C8, C9, C10, C11, TypeBorder4, TypeBorder2);
                    // относительно граничного условия и кубической интерполяцией задаём заграничные значения.
                    //граничные условия 
                    BorderZadelSharnkMVIW(A1, dx, n, TypeBorder4, TypeBorder2);
                    //нахождение разницы между предыдущей итерацией
                    MaxModule1 = MaximumModuleBetweenTwoMassiveOneDimensional(A1, AB1);

                    IntegrationL(C1, B2, B2, false, L1, F.F);
                    IntegrationL(C2, B2, B2, false, L2, F.F);
                    IntegrationL(C3, B2, B2, false, L3, F.F);
                    IntegrationL(C4, B2, B2, false, L4, F.F);
                    IntegrationL(C5, B2, B2, false, L5, F.F);
                    IntegrationL(C6, B2, B1, false, L1, F.F);
                    IntegrationL(C7, B2, B1, false, L2, F.F);
                    IntegrationL(C8, B2, B1, false, L3, F.F);
                    IntegrationL(C9, B2, B1, false, L4, F.F);
                    IntegrationL(C10, B2, B1, false, L5, F.F);
                    IntegrationL(C11, B2, B1, false, L6, F.F);

                    // находим решение уравнения методом гаусса жордана для функции А2
                    A2 = SolutionProblemMetodGaussaForKantorovichaVlasova(dx, A1, C1, C2, C3, C4, C5, C6, C7, C8, C9, C10, C11, TypeBorder4, TypeBorder2);
                    // относительно граничного условия и кубической интерполяцией задаём заграничные значения.
                    //граничные условия 
                    BorderZadelSharnkMVIW(A2, dx, n, TypeBorder4, TypeBorder2);
                    //нахождение разницы между предыдущей итерацией
                    MaxModule2 = MaximumModuleBetweenTwoMassiveOneDimensional(A2, AB2);

                    if (MaxModule1 >= MaxModule2)
                        if (MaxModule1 < 0.0000001) break;
                    if (MaxModule2 >= MaxModule1)
                        if (MaxModule2 < 0.0000001) break;
                    AB1 = ArrayAssignment(A1);
                    AB2 = ArrayAssignment(A2);

                }

                for (int i = 0; i < N + 2; i++)
                    for (int j = 0; j < M + 2; j++)
                    {
                        W[i, j] = A1[i + 1] * B1[j + 1] + A2[i + 1] * B2[j + 1];
                        if (Math.Abs(W[i, j]) > mW) mW = Math.Abs(W[i, j]); // Нахождение максимального прогиба
                    }
            }

            private void MethodVariationIterartionTwoLevel(Load F, int TypeBorder1, int TypeBorder2, int TypeBorder3, int TypeBorder4)
            {
                //Основные 2 Функции по каздой переменной
                //кол.т    0   1 .....  N-2  N-1(г)
                //A = 0 1 2(г) 3 ....    N  N+1(г) N+2 N+3
                //p = - - -(г) 0 1 2..  N-3  -(г)   -   -
                //DM = - 0 1(г) 2 ....   N-1  N(г)  N+1  - 

                //Переменные для записи максимальной разницы по норме, настоящего решения с предыдущим. 
                double MaxModule1 = 0;
                double MaxModule2 = 0;

                A1 = new double[N + 4];// по x
                B1 = new double[M + 4];// по y
                A2 = new double[N + 4];// по x
                B2 = new double[M + 4];// по y

                // для сравнения с предыдущим шагом итерации для решения системы диф уравнений
                AB1 = new double[N + 4];
                AB2 = new double[N + 4];

                // задаём размерность массивам по переменной х
                C1 = new double[N - 2];
                C2 = new double[N - 2];
                C3 = new double[N - 2];
                C4 = new double[N - 2];
                C5 = new double[N - 2];
                C6 = new double[N - 2];
                C7 = new double[N - 2];
                C8 = new double[N - 2];
                C9 = new double[N - 2];
                C10 = new double[N - 2];
                C11 = new double[N - 2];
                // задаём размерность массивам по переменной y
                H1 = new double[M - 2];
                H2 = new double[M - 2];
                H3 = new double[M - 2];
                H4 = new double[M - 2];
                H5 = new double[M - 2];
                H6 = new double[M - 2];
                H7 = new double[M - 2];
                H8 = new double[M - 2];
                H9 = new double[M - 2];
                H10 = new double[M - 2];
                H11 = new double[M - 2];

                // Задаём начальную функцию по y 
                for (int y = 0; y < M + 4; y++)
                {
                    B1[y] = (double)(Math.Sin(dy * (y - 2) * Math.PI));
                    B2[y] = (double)(Math.Sin(dy * (y - 2) * Math.PI));
                }
                for (int x = 0; x < N + 4; x++)
                    A2[x] = 0;

                //
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
                    //Решение системы дифференциальных уравнений относительно x
                    while (true)
                    {
                        // находим переменные функции стоящие при производных в ДУ.
                        IntegrationL(C1, B1, B1, false, L1, F.F);
                        IntegrationL(C2, B1, B1, false, L2, F.F);
                        IntegrationL(C3, B1, B1, false, L3, F.F);
                        IntegrationL(C4, B1, B1, false, L4, F.F);
                        IntegrationL(C5, B1, B1, false, L5, F.F);
                        IntegrationL(C6, B1, B2, false, L1, F.F);
                        IntegrationL(C7, B1, B2, false, L2, F.F);
                        IntegrationL(C8, B1, B2, false, L3, F.F);
                        IntegrationL(C9, B1, B2, false, L4, F.F);
                        IntegrationL(C10, B1, B2, false, L5, F.F);
                        IntegrationL(C11, B1, B2, false, L6, F.F);

                        // находим решение уравнения методом гаусса жордана
                        A1 = SolutionProblemMetodGaussaForKantorovichaVlasova(dx, A2, C1, C2, C3, C4, C5, C6, C7, C8, C9, C10, C11, TypeBorder4, TypeBorder2);
                        // относительно граничного условия и кубической интерполяцией задаём заграничные значения.
                        //граничные условия 
                        BorderZadelSharnkMVIW(A1, dx, n, TypeBorder4, TypeBorder2);
                        //нахождение разницы между предыдущей итерацией
                        MaxModule1 = MaximumModuleBetweenTwoMassiveOneDimensional(A1, AB1);

                        IntegrationL(C1, B2, B2, false, L1, F.F);
                        IntegrationL(C2, B2, B2, false, L2, F.F);
                        IntegrationL(C3, B2, B2, false, L3, F.F);
                        IntegrationL(C4, B2, B2, false, L4, F.F);
                        IntegrationL(C5, B2, B2, false, L5, F.F);
                        IntegrationL(C6, B2, B1, false, L1, F.F);
                        IntegrationL(C7, B2, B1, false, L2, F.F);
                        IntegrationL(C8, B2, B1, false, L3, F.F);
                        IntegrationL(C9, B2, B1, false, L4, F.F);
                        IntegrationL(C10, B2, B1, false, L5, F.F);
                        IntegrationL(C11, B2, B1, false, L6, F.F);

                        // находим решение уравнения методом гаусса жордана
                        A2 = SolutionProblemMetodGaussaForKantorovichaVlasova(dx, A1, C1, C2, C3, C4, C5, C6, C7, C8, C9, C10, C11, TypeBorder4, TypeBorder2);
                        // относительно граничного условия и кубической интерполяцией задаём заграничные значения.
                        //граничные условия 
                        BorderZadelSharnkMVIW(A2, dx, n, TypeBorder4, TypeBorder2);
                        //нахождение разницы между предыдущей итерацией
                        MaxModule2 = MaximumModuleBetweenTwoMassiveOneDimensional(A2, AB2);

                        if (MaxModule1 >= MaxModule2)
                            if (MaxModule1 < 0.0000001) break;
                        if (MaxModule2 >= MaxModule1)
                            if (MaxModule2 < 0.0000001) break;
                        AB1 = ArrayAssignment(A1);
                        AB2 = ArrayAssignment(A2);

                    }
                    //Решение системы дифференциальных уравнений относительно y

                    while (true)
                    {
                        // находим переменные функции стоящие при производных в ДУ.
                        IntegrationL(H1, A1, A1, true, L1, F.F);
                        IntegrationL(H2, A1, A1, true, L2, F.F);
                        IntegrationL(H3, A1, A1, true, L3, F.F);
                        IntegrationL(H4, A1, A1, true, L4, F.F);
                        IntegrationL(H5, A1, A1, true, L5, F.F);
                        IntegrationL(H6, A1, A2, true, L1, F.F);
                        IntegrationL(H7, A1, A2, true, L2, F.F);
                        IntegrationL(H8, A1, A2, true, L3, F.F);
                        IntegrationL(H9, A1, A2, true, L4, F.F);
                        IntegrationL(H10, A1, A2, true, L5, F.F);
                        IntegrationL(H11, A1, A2, true, L6, F.F);

                        // находим решение уравнения методом гаусса жордана
                        B1 = SolutionProblemMetodGaussaForKantorovichaVlasova(dy, B2, H1, H2, H3, H4, H5, H6, H7, H8, H9, H10, H11, TypeBorder3, TypeBorder1);
                        // относительно граничного условия и кубической интерполяцией задаём заграничные значения.
                        //граничные условия 
                        BorderZadelSharnkMVIW(B1, dy, n, TypeBorder3, TypeBorder1);
                        //нахождение разницы между предыдущей итерацией
                        MaxModule1 = MaximumModuleBetweenTwoMassiveOneDimensional(B1, AB1);

                        IntegrationL(H1, A2, A2, true, L1, F.F);
                        IntegrationL(H2, A2, A2, true, L2, F.F);
                        IntegrationL(H3, A2, A2, true, L3, F.F);
                        IntegrationL(H4, A2, A2, true, L4, F.F);
                        IntegrationL(H5, A2, A2, true, L5, F.F);
                        IntegrationL(H6, A2, A1, true, L1, F.F);
                        IntegrationL(H7, A2, A1, true, L2, F.F);
                        IntegrationL(H8, A2, A1, true, L3, F.F);
                        IntegrationL(H9, A2, A1, true, L4, F.F);
                        IntegrationL(H10, A2, A1, true, L5, F.F);
                        IntegrationL(H11, A2, A1, true, L6, F.F);

                        // находим решение уравнения методом гауHса жордана
                        B2 = SolutionProblemMetodGaussaForKantorovichaVlasova(dy, B1, H1, H2, H3, H4, H5, H6, H7, H8, H9, H10, H11, TypeBorder3, TypeBorder1);
                        // относительно граничного условия и кубической интерполяцией задаём заграничные значения.
                        //граничные условия 
                        BorderZadelSharnkMVIW(B2, dy, n, TypeBorder3, TypeBorder1);
                        //нахождение разницы между предыдущей итерацией
                        MaxModule2 = MaximumModuleBetweenTwoMassiveOneDimensional(B2, AB2);

                        if (MaxModule1 >= MaxModule2)
                            if (MaxModule1 < 0.0000001) break;
                        if (MaxModule2 >= MaxModule1)
                            if (MaxModule2 < 0.0000001) break;
                        AB1 = ArrayAssignment(A1);
                        AB2 = ArrayAssignment(A2);

                    }
                    CountIterationInMethods++;
                    for (int i = 0; i < N + 2; i++)
                        for (int j = 0; j < M + 2; j++)
                        {
                            if (Math.Abs(W[i, j] - (A1[i + 1] * B1[j + 1] + A2[i + 1] * B2[j + 1])) > MaxWForMethodIteration1) MaxWForMethodIteration1 = Math.Abs(W[i, j] - (A1[i + 1] * B1[j + 1] + A2[i + 1] * B2[j + 1]));// нахождение максмального отклонения найденного решения.
                            W[i, j] = A1[i + 1] * B1[j + 1] + A2[i + 1] * B2[j + 1];
                            if (Math.Abs(W[i, j]) > MaxWForMethodIteration2) MaxWForMethodIteration2 = Math.Abs(W[i, j]); // Нахождение максимального прогиба
                        }
                    if (Math.Abs(MaxWForMethodIteration1) < 0.000001) break;
                }
                MaxCountIterationSpendOfMetods += CountIterationInMethods;
                mW = MaxWForMethodIteration2;
            }
            ///..........................................................................................


            ///Метод Бубного-Галёркина
            public void StaticDecisionMethodBubnovGalercin(Load F, int DegreeOfDifficFunctionMethodBubnovaGalercin,
                int BorderExitMethodABS, int BorderExitPhisicalIteration, int TypeBorder1, int TypeBorder2, int TypeBorder3, int TypeBorder4)
            {
                RELoad();
                this.DegreeOfDifficFunctionMethodBubnovaGalercin = DegreeOfDifficFunctionMethodBubnovaGalercin;
                this.BorderExitMethodABS = BorderExitMethodABS;
                this.BorderExitPhisicalIteration = BorderExitPhisicalIteration;
                if (BorderExitMethodABS == 0)
                    PhisicalNolenearyPart(F, TypeBorder1, TypeBorder2, TypeBorder3, TypeBorder4, false, MethodBubnovGalerkin);
                else
                    PhisicalNolenearyPart(F, TypeBorder1, TypeBorder2, TypeBorder3, TypeBorder4, true, MethodBubnovGalerkin);
                //Заполнения значиний прогибов для печати
                for (int i = 0; i < N; i++)
                    for (int j = 0; j < M; j++)
                        PrintW[i, j] = W[i + 1, j + 1];

            }
            //Функция использующаяся в подсчёте
            private double[,] WBG;
            private double[,] WBGnm;
            private void FullingWBG(int TypeBorder3, int TypeBorder4, int i, int j)
            {
                WBG = new double[N + 4, M + 4];
                for (int x = 0; x < N + 4; x++)
                    for (int y = 0; y < M + 4; y++)
                    {
                        if (TypeBorder4 == 0)
                            WBG[x, y] = Math.Sin(Math.PI * i * (x - 2) * dx) * Math.Sin(Math.PI * i * (x - 2) * dx);
                        else
                            WBG[x, y] = Math.Sin(Math.PI * i * (x - 2) * dx);
                        if (TypeBorder3 == 0)
                            WBG[x, y] = WBG[x, y] * Math.Sin(Math.PI * j * (y - 2) * dy) * Math.Sin(Math.PI * j * (y - 2) * dy);
                        else
                            WBG[x, y] = WBG[x, y] * Math.Sin(Math.PI * j * (y - 2) * dy);
                    }
            }
            private void FullingWBGnm(int TypeBorder3, int TypeBorder4, int i, int j)
            {
                WBGnm = new double[N + 4, M + 4];
                for (int x = 0; x < N + 4; x++)
                    for (int y = 0; y < M + 4; y++)
                    {
                        if (TypeBorder4 == 0)
                            WBGnm[x, y] = Math.Sin(Math.PI * i * (x - 2) * dx) * Math.Sin(Math.PI * i * (x - 2) * dx);
                        else
                            WBGnm[x, y] = Math.Sin(Math.PI * i * (x - 2) * dx);
                        if (TypeBorder3 == 0)
                            WBGnm[x, y] = WBGnm[x, y] * Math.Sin(Math.PI * j * (y - 2) * dy) * Math.Sin(Math.PI * j * (y - 2) * dy);
                        else
                            WBGnm[x, y] = WBGnm[x, y] * Math.Sin(Math.PI * j * (y - 2) * dy);
                    }
            }

            //Оператор L из уравнения пластинки
            private double LBG(int x, int y)
            {

                lamd2 = lam * lam;
                lamd4 = (double)Math.Pow(lam, 4.0);
                D = DM[x, y]; R = RM[x, y]; S = SM[x, y];
                Ddx = (DM[x + 1, y] - DM[x - 1, y]) / (2 * dx);
                Ddy = (DM[x, y + 1] - DM[x, y - 1]) / (2 * dy);
                Dd2x = (DM[x + 1, y] - 2 * DM[x, y] + DM[x - 1, y]) / (dx * dx);
                Dd2y = (DM[x, y + 1] - 2 * DM[x, y] + DM[x, y - 1]) / (dy * dy);
                Rd2x = (RM[x + 1, y] - 2 * RM[x, y] + RM[x - 1, y]) / (dx * dx);
                Rd2y = (RM[x, y + 1] - 2 * RM[x, y] + RM[x, y - 1]) / (dy * dy);
                Rdx = (RM[x + 1, y] - RM[x - 1, y]) / (2 * dx);
                Rdy = (RM[x, y + 1] - RM[x, y - 1]) / (2 * dy);
                Sdx = (SM[x + 1, y] - SM[x - 1, y]) / (2 * dx);
                Sdy = (SM[x, y + 1] - SM[x, y - 1]) / (2 * dy);
                S2dxdy = (SM[x + 1, y + 1] - SM[x + 1, y - 1] - SM[x - 1, y + 1] + SM[x - 1, y - 1]) / (4 * dx * dy);

                x++;
                y++;
                Wd4x = (WBG[x + 2, y] - 4 * WBG[x + 1, y] + 6 * WBG[x, y] - 4 * WBG[x - 1, y] + WBG[x - 2, y]) / (double)Math.Pow(dx, 4.0);
                Wd4y = (WBG[x, y + 2] - 4 * WBG[x, y + 1] + 6 * WBG[x, y] - 4 * WBG[x, y - 1] + WBG[x, y - 2]) / (double)Math.Pow(dy, 4.0);
                Wd3x = (WBG[x + 2, y] - 2 * WBG[x + 1, y] + 2 * WBG[x - 1, y] - WBG[x - 2, y]) / (2 * (double)Math.Pow(dx, 3.0));
                Wd3y = (WBG[x, y + 2] - 2 * WBG[x, y + 1] + 2 * WBG[x, y - 1] - WBG[x, y - 2]) / (2 * (double)Math.Pow(dy, 3.0));
                Wd2x = (WBG[x + 1, y] - 2 * WBG[x, y] + WBG[x - 1, y]) / (double)Math.Pow(dx, 2.0);
                Wd2y = (WBG[x, y + 1] - 2 * WBG[x, y] + WBG[x, y - 1]) / (double)Math.Pow(dy, 2.0);
                Wd2xd2y = (WBG[x + 1, y + 1] + WBG[x + 1, y - 1] + WBG[x - 1, y + 1] + WBG[x - 1, y - 1]
                    - 2 * WBG[x, y + 1] - 2 * WBG[x, y - 1] - 2 * WBG[x - 1, y] - 2 * WBG[x + 1, y] + 4 * WBG[x, y]) / (dx * dx * dy * dy);
                Wd2ydx = (WBG[x + 1, y + 1] + WBG[x + 1, y - 1] - WBG[x - 1, y + 1] - WBG[x - 1, y - 1]
                    + 2 * WBG[x - 1, y] - 2 * WBG[x + 1, y]) / (2 * dx * dy * dy);
                Wd2xdy = (WBG[x + 1, y + 1] - WBG[x + 1, y - 1] + WBG[x - 1, y + 1] - WBG[x - 1, y - 1]
                    - 2 * WBG[x, y + 1] + 2 * WBG[x, y - 1]) / (2 * dx * dx * dy);
                Wdydx = (WBG[x + 1, y + 1] - WBG[x + 1, y - 1] - WBG[x - 1, y + 1] + WBG[x - 1, y - 1]) / (4 * dx * dy);

                return (Wd4x * D + 2 * Wd3x * Ddx + (Dd2x + lamd2 * Rd2y) * Wd2x +
                    lamd4 * (Wd4y * D + 2 * Wd3y * Ddy + (Dd2y + Rd2x / lamd2) * Wd2y) +
                    2 * lamd2 * (Wd2xd2y * (R + S) + Wd2ydx * (Rdx + Sdx) + Wd2xdy * (Rdy + Sdy) + Wdydx * S2dxdy)) * WBGnm[x, y];

            }
            //Оператор U
            private double UBG(double[,] F, int x, int y)
            {
                if (alf != 0)
                {
                    Id2x = (IM[x + 1, y] - 2 * IM[x, y] + IM[x - 1, y]) / (dx * dx);
                    Id2y = (IM[x, y + 1] - 2 * IM[x, y] + IM[x, y - 1]) / (dy * dy);
                }
                return 2 * lamd2 * Id2y + Id2x - F[x - 1, y - 1];
            }

            //Оператор U
            private void MethodBubnovGalerkin(Load F, int TypeBorder1, int TypeBorder2, int TypeBorder3, int TypeBorder4)
            {
                //C*C колличество уравнений с коэффициентами
                int C = DegreeOfDifficFunctionMethodBubnovaGalercin;
                //Матрица для расчёта коэффициентов методом жордана гаусса
                double[,] MatrixCoefficient = new double[C * C, C * C + 1];
                //Нулевые значиения расчитываемых интегралов.
                double U = 0, L = 0;
                int i = 0, j = 0;
                for (int n1 = 1; n1 <= C; n1++)
                    for (int m1 = 1; m1 <= C; m1++)
                    {
                        FullingWBGnm(TypeBorder3, TypeBorder4, n1, m1);
                        //ShowMassiv2(WBGnm);
                        for (int i1 = 1; i1 <= C; i1++)
                        {
                            for (int j1 = 1; j1 <= C; j1++)
                            {
                                FullingWBG(TypeBorder3, TypeBorder4, i1, j1);
                                //ShowMassiv2(WBG);C
                                for (int x = 1; x < N; x++)
                                    for (int y = 1; y < M; y++)
                                    {
                                        L = L + (LBG(x, y) + LBG(x + 1, y) + LBG(x, y + 1) + LBG(x + 1, y + 1)) * dx * dy / 4;
                                    }
                                MatrixCoefficient[i, j] = L;
                                L = 0;
                                j++;
                            }

                        }
                        for (int x = 1; x < N; x++)
                            for (int y = 1; y < M; y++)
                                U = U + (UBG(F.F, x, y) * WBGnm[x + 1, y + 1] + UBG(F.F, x + 1, y) * WBGnm[x + 2, y + 1] +
                                    UBG(F.F, x, y + 1) * WBGnm[x + 1, y + 2] + UBG(F.F, x + 1, y + 1) * WBGnm[x + 2, y + 2]) * dx * dy / 4;
                        MatrixCoefficient[i, j] = -U;
                        U = 0;
                        j = 0;
                        i++;
                    }

                MethodGordanGauss(MatrixCoefficient);
                for (int x = 0; x < N + 2; x++)
                    for (int y = 0; y < M + 2; y++)
                        W[x, y] = 0;
                int I = 0;
                for (i = 1; i <= C; i++)
                    for (j = 1; j <= C; j++)
                    {
                        FullingWBG(TypeBorder3, TypeBorder4, i, j);
                        for (int x = 0; x < N + 2; x++)
                            for (int y = 0; y < M + 2; y++)
                                W[x, y] = W[x, y] + (double)(MatrixCoefficient[I, C * C] * WBG[x + 1, y + 1]);
                        I++;
                    }
                //ShowMassiv2(MatrixCoefficient);
                mW = 0;
                for (int x = 0; x < N + 2; x++)
                    for (int y = 0; y < M + 2; y++)
                        if (mW < Math.Abs(W[x, y])) mW = Math.Abs(W[x, y]);

            }
            ///....................................

            /// Метод Аграновсвкого Баглая Смирнова.........................................

            // Счётчик метода  Аграновсвкого Баглая Смирнова
            // private int CountIterationABS;
            private Load Neviaska;
            private double MaxForABS1, MaxForABS2;
            // Уравнение пластинки  для Аграновского Баглая Смирнова
            public double EquationForMethodABS(double[,] F, int x, int y)
            {
                lamd2 = lam * lam;
                lamd4 = (double)Math.Pow(lam, 4.0);
                D = DM[x, y]; R = RM[x, y]; S = SM[x, y];
                Ddx = (DM[x + 1, y] - DM[x - 1, y]) / (2 * dx);
                Ddy = (DM[x, y + 1] - DM[x, y - 1]) / (2 * dy);
                Dd2x = (DM[x + 1, y] - 2 * DM[x, y] + DM[x - 1, y]) / (dx * dx);
                Dd2y = (DM[x, y + 1] - 2 * DM[x, y] + DM[x, y - 1]) / (dy * dy);
                Rd2x = (RM[x + 1, y] - 2 * RM[x, y] + RM[x - 1, y]) / (dx * dx);
                Rd2y = (RM[x, y + 1] - 2 * RM[x, y] + RM[x, y - 1]) / (dy * dy);
                Rdx = (RM[x + 1, y] - RM[x - 1, y]) / (2 * dx);
                Rdy = (RM[x, y + 1] - RM[x, y - 1]) / (2 * dy);
                Sdx = (SM[x + 1, y] - SM[x - 1, y]) / (2 * dx);
                Sdy = (SM[x, y + 1] - SM[x, y - 1]) / (2 * dy);
                S2dxdy = (SM[x + 1, y + 1] - SM[x + 1, y - 1] - SM[x - 1, y + 1] + SM[x - 1, y - 1]) / (4 * dx * dy);
                if (alf != 0)
                {
                    Id2x = (IM[x + 1, y] - 2 * IM[x, y] + IM[x - 1, y]) / (dx * dx);
                    Id2y = (IM[x, y + 1] - 2 * IM[x, y] + IM[x, y - 1]) / (dy * dy);
                }
                return (double)(((W[x + 2, y] - 4 * W[x + 1, y] + 6 * W[x, y] - 4 * W[x - 1, y] + W[x - 2, y]) * D / Math.Pow(dx, 4.0) +

                    (W[x + 2, y] - 2 * W[x + 1, y] + 2 * W[x - 1, y] - W[x - 2, y]) * Ddx / Math.Pow(dx, 3.0) +

                    (W[x + 1, y] - 2 * W[x, y] + W[x - 1, y]) * (Dd2x + lamd2 * Rd2y) / (dx * dx) +

                    (W[x, y + 2] - 4 * W[x, y + 1] + 6 * W[x, y] - 4 * W[x, y - 1] + W[x, y - 2]) * lamd4 * D / Math.Pow(dy, 4.0) +

                    (W[x, y + 2] - 2 * W[x, y + 1] + 2 * W[x, y - 1] - W[x, y - 2]) * lamd4 * Ddy / Math.Pow(dy, 3.0) +

                    (W[x, y + 1] - 2 * W[x, y] + W[x, y - 1]) * lamd4 * (Dd2y + Rd2x / lamd2) / (dy * dy) +

                    (W[x + 1, y + 1] + W[x + 1, y - 1] + W[x - 1, y + 1] + W[x - 1, y - 1]
                    - 2 * W[x, y + 1] - 2 * W[x, y - 1] - 2 * W[x - 1, y] - 2 * W[x + 1, y] + 4 * W[x, y]) * 2 * lamd2 * (R + S) / (dx * dx * dy * dy) +

                    (W[x + 1, y + 1] + W[x + 1, y - 1] - W[x - 1, y + 1] - W[x - 1, y - 1]
                    + 2 * W[x - 1, y] - 2 * W[x + 1, y]) * lamd2 * (Rdx + Sdx) / (dx * dy * dy) +

                    (W[x + 1, y + 1] - W[x + 1, y - 1] + W[x - 1, y + 1] - W[x - 1, y - 1]
                    - 2 * W[x, y + 1] + 2 * W[x, y - 1]) * lamd2 * (Rdy + Sdy) / (dx * dx * dy) +

                    (W[x + 1, y + 1] - W[x + 1, y - 1] - W[x - 1, y + 1] + W[x - 1, y - 1]) * lamd2 * S2dxdy / (2 * dx * dy) +

                    2 * lamd2 * Id2y + Id2x - F[x - 1, y - 1]));
            }
            public double EquationForMethodABS1(int x, int y)
            {
                lamd2 = lam * lam;
                lamd4 = (double)Math.Pow(lam, 4.0);
                D = DM[x, y]; R = RM[x, y]; S = SM[x, y];
                Ddx = (DM[x + 1, y] - DM[x - 1, y]) / (2 * dx);
                Ddy = (DM[x, y + 1] - DM[x, y - 1]) / (2 * dy);
                Dd2x = (DM[x + 1, y] - 2 * DM[x, y] + DM[x - 1, y]) / (dx * dx);
                Dd2y = (DM[x, y + 1] - 2 * DM[x, y] + DM[x, y - 1]) / (dy * dy);
                Rd2x = (RM[x + 1, y] - 2 * RM[x, y] + RM[x - 1, y]) / (dx * dx);
                Rd2y = (RM[x, y + 1] - 2 * RM[x, y] + RM[x, y - 1]) / (dy * dy);
                Rdx = (RM[x + 1, y] - RM[x - 1, y]) / (2 * dx);
                Rdy = (RM[x, y + 1] - RM[x, y - 1]) / (2 * dy);
                Sdx = (SM[x + 1, y] - SM[x - 1, y]) / (2 * dx);
                Sdy = (SM[x, y + 1] - SM[x, y - 1]) / (2 * dy);
                S2dxdy = (SM[x + 1, y + 1] - SM[x + 1, y - 1] - SM[x - 1, y + 1] + SM[x - 1, y - 1]) / (4 * dx * dy);
                if (alf != 0)
                {
                    Id2x = (IM[x + 1, y] - 2 * IM[x, y] + IM[x - 1, y]) / (dx * dx);
                    Id2y = (IM[x, y + 1] - 2 * IM[x, y] + IM[x, y - 1]) / (dy * dy);
                }
                return (double)(((W[x + 2, y] - 4 * W[x + 1, y] + 6 * W[x, y] - 4 * W[x - 1, y] + W[x - 2, y]) * D / Math.Pow(dx, 4.0) +

                    (W[x + 2, y] - 2 * W[x + 1, y] + 2 * W[x - 1, y] - W[x - 2, y]) * Ddx / Math.Pow(dx, 3.0) +

                    (W[x + 1, y] - 2 * W[x, y] + W[x - 1, y]) * (Dd2x + lamd2 * Rd2y) / (dx * dx) +

                    (W[x, y + 2] - 4 * W[x, y + 1] + 6 * W[x, y] - 4 * W[x, y - 1] + W[x, y - 2]) * lamd4 * D / Math.Pow(dy, 4.0) +

                    (W[x, y + 2] - 2 * W[x, y + 1] + 2 * W[x, y - 1] - W[x, y - 2]) * lamd4 * Ddy / Math.Pow(dy, 3.0) +

                    (W[x, y + 1] - 2 * W[x, y] + W[x, y - 1]) * lamd4 * (Dd2y + Rd2x / lamd2) / (dy * dy) +

                    (W[x + 1, y + 1] + W[x + 1, y - 1] + W[x - 1, y + 1] + W[x - 1, y - 1]
                    - 2 * W[x, y + 1] - 2 * W[x, y - 1] - 2 * W[x - 1, y] - 2 * W[x + 1, y] + 4 * W[x, y]) * 2 * lamd2 * (R + S) / (dx * dx * dy * dy) +

                    (W[x + 1, y + 1] + W[x + 1, y - 1] - W[x - 1, y + 1] - W[x - 1, y - 1]
                    + 2 * W[x - 1, y] - 2 * W[x + 1, y]) * lamd2 * (Rdx + Sdx) / (dx * dy * dy) +

                    (W[x + 1, y + 1] - W[x + 1, y - 1] + W[x - 1, y + 1] - W[x - 1, y - 1]
                    - 2 * W[x, y + 1] + 2 * W[x, y - 1]) * lamd2 * (Rdy + Sdy) / (dx * dx * dy) +

                    (W[x + 1, y + 1] - W[x + 1, y - 1] - W[x - 1, y + 1] + W[x - 1, y - 1]) * lamd2 * S2dxdy / (2 * dx * dy) +

                    2 * lamd2 * Id2y + Id2x));
            }
            private void MethodABS(Load F, int TypeBorder1, int TypeBorder2, int TypeBorder3, int TypeBorder4, DelegatForMethods Method)
            {
                double[,] NenugniMassiv = new double[N, M];
                Method(F, TypeBorder1, TypeBorder2, TypeBorder3, TypeBorder4);

                Neviaska = new Load(N, M);
                W1 = new double[N + 2, M + 2];

                MaxForABS1 = 0;
                MaxForABS2 = 0;
                //Console.WriteLine("T0");
                //ShowMassiv2(W);
                //Console.WriteLine("  " + mW);
                for (int i = 1; i <= BorderExitMethodABS; i++)
                {
                    MaxForABS1 = 0;
                    //Считаем невязку
                    //Console.WriteLine(i);
                    for (int x = 1; x < N - 1; x++)
                        for (int y = 1; y < M - 1; y++)
                        {
                            Neviaska.F[x, y] = -EquationForMethodABS(F.F, x + 1, y + 1);
                            //NenugniMassiv[x,y] = EquationForMethodABS1(x + 1, y + 1);
                        }


                    //Console.WriteLine("Neviazka");
                    //ShowMassiv2(Neviaska.F);
                    //Console.WriteLine("W");
                    //ShowMassiv2(NenugniMassiv);
                    //ShowMassiv2(W);

                    //Копируем прогиб
                    for (int x = 0; x < N + 2; x++)
                        for (int y = 0; y < M + 2; y++)
                            W1[x, y] = W[x, y];
                    //Считаем добавочную функцию тем же методом
                    Method(Neviaska, TypeBorder1, TypeBorder2, TypeBorder3, TypeBorder4);
                    //Console.WriteLine("T1");
                    //ShowMassiv2(W);

                    //Считается общий прогиб
                    for (int x = 0; x < N + 2; x++)
                        for (int y = 0; y < M + 2; y++)
                        {

                            W[x, y] = W[x, y] + W1[x, y];
                            if (Math.Abs(W[x, y] - W1[x, y]) > MaxForABS1) MaxForABS1 = Math.Abs(W[x, y] - W1[x, y]);// нахождение максмального отклонения найденного решения.
                            if (Math.Abs(W[x, y]) > MaxForABS2) MaxForABS2 = Math.Abs(W[x, y]); // Нахождение максимального прогиба
                        }
                    //Console.WriteLine("T0+T1");
                    //ShowMassiv2(W);

                    if (Math.Abs(MaxForABS1) < 0.000001) break;
                    CountIterationABS++;
                }
                //CountIterationInMethods += CountIterationABS;
                mW = MaxForABS2;
            }
            /// ........................................................................
            /// Метод пересчёта переменных параметров упрогости
            public double CountIterationPhisicalNoneleneary { get; set; } //Количесвто итераций затраченных на подсчёт переменных параметров упругости
            public double MaxCountMethodOfIteration { get; set; } // максимальное значение итераций затрачиваемых итерационным методом
            private void PhisicalNolenearyPart(Load F, int TypeBorder1, int TypeBorder2, int TypeBorder3, int TypeBorder4, bool FlagABS, DelegatForMethods Method)
            {
                double[,] Wlast = new double[N + 2, M + 2];
                double epsilon = 0;
                for (int i = 0; i < N + 2; i++)
                    for (int j = 0; j < M + 2; j++)
                        Wlast[i, j] = 0;
                CountIterationABS = 0;
                CountIterationPhisicalNoneleneary = 0;
                for (int fizi = 1; fizi <= BorderExitPhisicalIteration; fizi++)
                {

                    MaxWForMethodIteration1 = 0;
                    MaxWForMethodIteration2 = 0;
                    // Итерационный метод
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
                    //Максимальное колличество итераци затрачиваемых итерационным методом
                    //if (MaxCountMethodOfIteration <CountIterationInMethods) MaxCountMethodOfIteration = CountIterationInMethods;
                    //Счётчик подсчёта переменных параметров упругости
                    CountIterationPhisicalNoneleneary++;
                    //Console.WriteLine(CountIterationPhisicalNoneleneary);
                    for (int i = 0; i < N + 2; i++)
                        for (int j = 0; j < M + 2; j++)
                        {
                            if (Math.Abs(W[i, j] - Wlast[i, j]) > MaxWForMethodIteration1) MaxWForMethodIteration1 = Math.Abs(W[i, j] - Wlast[i, j]);// нахождение максмального отклонения найденного решения.
                            Wlast[i, j] = W[i, j];
                            //if (Math.Abs(W[i, j]) > MaxWForMethodIteration2) MaxWForMethodIteration2 = Math.Abs(W[i, j]); // Нахождение максимального прогиба
                        }

                    if (Math.Abs(MaxWForMethodIteration1) < epsilon) break;

                    //ShowMassiv2(W);
                    //Console.WriteLine(mW+" "+ CountIterationPhisicalNoneleneary);
                    // пересчёт всех значений с новыми прогибами
                    Load();
                    //Console.WriteLine(fizi);

                }
            }
            /// ...............................................................................

        }
        static void WriteMassivInFile(double[,] Mass, string j)
        {
            System.Threading.Thread.CurrentThread.CurrentCulture = new System.Globalization.CultureInfo("en-US");// смена точки на запятую
            FileStream file1 = new FileStream("File_" + j.ToString() + ".txt", FileMode.Create); //создаем файловый поток
            StreamWriter writer = new StreamWriter(file1);//создаем «потоковый писатель» и связываем его с файловым потоком 

            //Вывод массива данных в файл
            for (int i = 0; i < Mass.GetLength(0); i++)
            {
                for (int k = 0; k < Mass.GetLength(1); k++)
                {
                    writer.Write(String.Format("{0:0.00000000 }", Mass[i, k]));
                }
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
            WriteMassivInFile(Exit, s);
        }
       


        static void Main(string[] args)
        {
            // наилучшие 30 30 6
            // начальные данные 
            int N = 30;
            int M = 30;
            int P = 4;//(2*P-1)   
            double n = 1;
            double m = 1;
            double p = 0.05; //0.05
            double nu = 0;// коэффициент Пуассона
            double sigmas = 0; //предел текучести 
            double es = 0; //деформация текучести она получается за счёт предела текучести и модуля сдвига es=sigmas/(3*G0) 
            double G0 = 0;//коэффициент модуля сдвига для сжатия
            double QTemp = 0; //Температурная нагрузка
            double alf = 0; //коэффициент линейного расширения
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
            //Ессли задача связанна со временем
            double t0 = 0;
            double t1 = 0;
            double dissipation = 0;
            //градиент с материалом на ходящийся на нижней грани пластины
            double nu1 = 0.24; //коэффициент Пуассонта материала
            double E1 = 320000; //модуль юнга материала
            double RaspredGradient = 0; // коэффициент определяющий тип распределение градиента.
            double Poristostb = 0; // коэффициент пористости
            int ModelPoristostb = 0; // модель пористости (4)
            int ModelFunctionСurvilinearPlane = 0; // модель криволиненой поверхности пластинки
            // определяющие параметры используемы при счёте
            int model = 20;
            // тип граничных условий
            int typeGrandIF = 1;
            //Тип нагрузки
            int typeQ = 1;
            // Метод 
            int TypeMethod = 1;
            // Слой печатания
            int sloi = 1;
            //количество коэф бубнова галёркин
            int a = 9;
            // данные по керамике
            //Режим счёта с одной нагрузкой и множеством
            bool flagNagruzki = true;
            int U1 = 0; // начальное знавение нагрузки 
            int U2 = 400; // коненчное значение нагрузки 
            double Q = 130; // шаг нагрузки
            // https://www.lib.tpu.ru/fulltext/v/Bulletin_TPU/2006/v309/i2/06.pdf

            // создаю объект класса БАЛКА
            Plast plast = new Plast(N, M, P, n, m, p);
            // типы исследуемых моделей
            switch (model)
            {
                case 1:
                    //allum
                    nu = 0.314;
                    G0 = 26490;
                    sigmas = 98.0665;
                    es = 0.001234;
                    l = 0.5;
                    plast.InicializationPlast(nu, G0);
                    plast.InicializationExpFizicalNonlinProblem(sigmas, es);
                    plast.InicializationNano(l);
                    break;
                case 2:
                    //steel T=0
                    // для экспоненциального!!!
                    G0 = 83250;
                    sigmas = 728.522;
                    es = 0.002917;
                    nu = 0.3;
                    plast.InicializationPlast(nu, G0);
                    plast.InicializationExpFizicalNonlinProblem(sigmas, es);
                    break;
                case 3:
                    //steel T=300;
                    // для экспоненциального!!!
                    G0 = 74996;//коэффициент модуля сдвига для сжатия
                    sigmas = 630.642;
                    es = 0.002803;
                    nu = 0.3;// коэффициент Пуассона
                    QTemp = 300; //Температурная нагрузка
                    alf = 0.000017; //коэффициент линейного температурного расширения
                    plast.InicializationPlast(nu, G0);
                    plast.InicializationExpFizicalNonlinProblem(sigmas, es);
                    plast.InicializationTempature(QTemp, alf);
                    break;
                case 4:
                    //steel T=500;
                    // для экспоненциального!!!
                    G0 = 62755;
                    sigmas = 520.555;
                    es = 0.002765;
                    nu = 0.3;
                    QTemp = 500; //Температурная нагрузка
                    alf = 0.000018; //коэффициент линейного температурного расширения
                    plast.InicializationPlast(nu, G0);
                    plast.InicializationExpFizicalNonlinProblem(sigmas, es);
                    plast.InicializationTempature(QTemp, alf);
                    break;
                case 5:
                    //steel T=0
                    // для экспоненциального!!!
                    G0 = 83250;
                    sigmas = 728.522;
                    es = 0.002917; 
                    nu = 0.3;
                    l = 0.5;
                    plast.InicializationPlast(nu, G0);
                    plast.InicializationExpFizicalNonlinProblem(sigmas, es);
                    plast.InicializationNano(l);

                    break;
                case 6:
                    //steel T=300;
                    // для экспоненциального!!!
                    G0 = 74996;//коэффициент модуля сдвига для сжатия
                    sigmas = 630.642;
                    es = 0.002803;
                    nu = 0.3;// коэффициент Пуассона
                    QTemp = 300; //Температурная нагрузка
                    alf = 0.000017; //коэффициент линейного температурного расширения
                    l = 0.5;
                    plast.InicializationPlast(nu, G0);
                    plast.InicializationExpFizicalNonlinProblem(sigmas, es);
                    plast.InicializationNano(l);
                    plast.InicializationTempature(QTemp, alf);
                    break;
                case 7:
                    //steel T=500;
                    // для экспоненциального!!!
                    G0 = 62755;
                    sigmas = 520.555;
                    es = 0.002765;
                    nu = 0.3;
                    QTemp = 500; //Температурная нагрузка
                    alf = 0.000018; //коэффициент линейного температурного расширения
                    l = 0.5;
                    plast.InicializationPlast(nu, G0);
                    plast.InicializationExpFizicalNonlinProblem(sigmas, es);
                    plast.InicializationNano(l);
                    plast.InicializationTempature(QTemp, alf);
                    break;
                case 8:
                    //steel T=0
                    // для линейного упрочнения
                    G0 = 75000;
                    G1 = 2211.53;
                    es = 0.002917;
                    nu = 0.3;
                    plast.InicializationPlast(nu, G0);
                    plast.InicializationLinearyHard(G1, es);
                    break;
                case 9:
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
                case 10:
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
                case 11:
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
                case 12:
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
                case 13:
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
                case 14:
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
                case 15:
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
                case 16:
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
                case 17:
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
                case 18:
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
                case 19:
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
                case 20:
                    //steel T=300;
                    // для экспоненциального!!!
                    G0 = 74996;//коэффициент модуля сдвига для сжатия
                    sigmas = 630.642;
                    es = 0.002803;
                    nu = 0.3;// коэффициент Пуассона
                    QTemp = 300; //Температурная нагрузка
                    alf = 0.000017; //коэффициент линейного температурного расширения
                    plast.InicializationPlast(nu, G0);
                    plast.InicializationExpFizicalNonlinProblem(sigmas, es);
                    plast.InicializationTempature("тем_500.txt", alf);
                    break;
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
            switch(typeGrandIF)
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
            string Form = String.Format("TQ{0}_TIF{1}_bg{2}_Met{3}_n{4}m{5}p{6}_Q{7}_N{8}M{9}P{10}_t0{11}t1{12}dis{13}_L{14}_model{15}_RG{16}_Modpor{17}_MFP{18}",
                        typeQ, typeGrandIF, a, TypeMethod, n, m, p, Q, N, M, 2 * P - 1, t0, t1, dissipation, l, model, RaspredGradient, ModelPoristostb, ModelFunctionСurvilinearPlane);
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
                            plast.StaticDecisionMethodVariationIteration(F, 1000, 0, 1000, TypeBorder1, TypeBorder2, TypeBorder3, TypeBorder4);
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
                    //for (int x = 1; x < N + 1; x++)
                    //    PrintWX[x - 1] = plast.W[x, M / 2];

                    //// эпюр моментов
                    //for (int x = 1; x < N + 1; x++)
                    //    PrintW2X[x - 1] = (plast.W[x - 1, M / 2] - 2 * plast.W[x, M / 2] + plast.W[x + 1, M / 2]) / (dx * dx);


                    //for (int y = 1; y < M + 1; y++)
                    //    PrintW2Y[y - 1] = (plast.W[N / 2, y - 1] - 2 * plast.W[N / 2, y] + plast.W[N / 2, y + 1]) / (dy * dy);

                    

                    WriteNumberInFile(plast.MaximumCountIterationOfMethodAllDecision, "IM_" + Form);

                    WriteNumberInFile(plast.CountIterationPhisicalNoneleneary, "FI_" + Form);

                    WriteNumberInFile(plast.mW, "mW_" + Form);

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
                                    WriteNumberInFile(Q*i, "QPSh_"+Form);
                                    WriteNumberInFile(plast.mW, "QPSh_" + Form);
                                    WriteMassivInFile(plast.Eii, 0, String.Format("EiiPSh_{0}_", 0) + Form);
                                    WriteMassivInFile(plast.Eii, sloi, String.Format("EiiPSh_{0}_", sloi) + Form);
                                    WriteMassivInFile(plast.Eii, P - 3, String.Format("EiiPSh_{0}_", P - 3) + Form);
                                    double[,] granb = new double[P * 2 - 1, N];
                                    for (int i2 = 0; i2 < P * 2 - 1; i2++)
                                        for (int k2 = 0; k2 < N; k2++)
                                            granb[i2, k2] = plast.Eii[k2, 0, i2];
                                    WriteMassivInFile(granb, "granb 0");
                                }
                            }
                    }
                    //Печать данных для определённой нагрузки
                    if (flagNagruzki)
                    {
                        WriteNumberInFile(Q * i, "Q_" + Form);
                        WriteNumberInFile(plast.mW, "Q_" + Form);
                        WriteMassivInFile(plast.Eii, 0, String.Format("Eii_{0}_", 0) + Form);
                        WriteMassivInFile(plast.Eii, sloi, String.Format("Eii_{0}_", sloi) + Form);
                        WriteMassivInFile(plast.Eii, P - 2, String.Format("Eii_{0}_", P - 2) + Form);
                        Console.WriteLine(Q * i);
                        double[,] granb = new double[P * 2 - 1, N];
                        for (int i2 = 0; i2 < P * 2 - 1; i2++)
                            for (int k2 = 0; k2 < N; k2++)
                                granb[i2, k2] = plast.Eii[k2, 0, i2];
                        WriteMassivInFile(granb, "granb 0");
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
                                        WriteNumberInFile(Q * i, "QFN_" + Form);
                                        WriteNumberInFile(plast.mW, "QFN_" + Form);
                                        WriteMassivInFile(plast.Eii, 0, String.Format("EiiFN_{0}_", 0) + Form);
                                        WriteMassivInFile(plast.Eii, sloi, String.Format("EiiFN_{0}_", sloi) + Form);
                                        WriteMassivInFile(plast.Eii, P - 3, String.Format("EiiFN_{0}_", P - 3) + Form);

                                    }
                                }
                            }
                        }
                    }
                    //условие выхода из цикла расчёта
                    if (plast.CountIterationPhisicalNoneleneary >= 1000 || plast.mW > 0.3) break;
                    plast.RELoad();
                }
            }
            stopwatch.Stop();
            Console.WriteLine("Затрачиваемое время " + stopwatch.ElapsedMilliseconds + " мс");
            Console.WriteLine(Form);
            Console.WriteLine("КОНЕЦ ПРОГРАММЫ");
            Console.ReadKey();
        }
    }
}