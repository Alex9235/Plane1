using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace Decoder
{
    public class Plane
    {
        //количество точек
        //по оси x
        protected int N { get; set; }
        //по оси y
        protected int M { get; set; }
        //по оси z
        protected int P { get; set; }
        //размеры пластины в пространстве
        protected double n { get; set; }
        protected double m { get; set; }
        protected double p { get; set; }

        //шаг разбиение 
        //по оси x
        protected double dx { get; set; }
        //по оси y
        protected double dy { get; set; }
        //по оси z
        protected double dz { get; set; }
        //переменные используемые в рассчётах
        public double lam1 { get; protected set; } // отношение длины к ширине
        public double lam2 { get; protected set; } // отношение длины к толщине
        public double lam { get; protected set; } // отношение длины к толщине
        public double[,,] NU { get; protected set; } //Коэффициент Пуассона
        public double[,,] K1 { get; protected set; }// коэффициент объемного растяжения 
        public double Gk1 { get; protected set; }// коэффициент линейности без учёта физ нелинейности расятжений
        public double G0 { get; protected set; } // модуль сдвига  
        public double[,,] E { get; protected set; }// Модуль Юнга
        public double[,] W { get; protected set; }// Прогиб
        public double[,] PrintW { get; protected set; }// Для вывода
        public double[,,] Eii { get; protected set; } // массив деформаций
        public double mW { get; set; } //максимальное значение прогиба
        protected double nu;  //коэффициент пуассона

        //Экспоненциальная зависимость напряжения от деформаций
        protected bool FlagLinearyAndExponentModelSigmaEpsilon1;// Учёт экспоненциальной зависимости напряжения от деформаций
        public double SigmS;  //предел текучести
        public double eS; //деформация текучести
        //для случая с криволинейной поверхностью
        protected double[,] DZ { get; set; }

        //Связь с температурой
        protected bool FlagTemperatureProblem;
        protected double[,] DeformMapT; //Массив напряжений зависимых от диформаций и температуры
        protected bool FlagDeformMapT; //флаг учитывающий связь температуры
        protected double Temp0; //миниальное значение температуры в зависимости s(e,T)
        protected double Temp1; //максимальное значение температуры в зависимости s(e,T)
        protected double ei0; //миниальное значение деформаций в зависимости s(e,T)
        protected double ei1; //максимальное значение деформаций в зависимости s(e,T)
        public double[,,] T { get; protected set; } // температуное поле
        protected double[,,] Temp; //массив внутренннего источника тепла
        protected double[,,] alf; //коэффициент линейного температурного расширения
        protected double[,,] G000; //Если модуль сдига зависит от температуры
        protected double[,,] G001; //Если модуль сдига зависит от температуры
        protected double[,,] SigmaT1;
        protected double[,,] SigmaT2;
        protected bool FlagTemperatureWithENUG0;

        //Влажность
        protected bool FlagHumidity;
        public double[,,] C { get; protected set; } // Влажность
        protected double[,,] beta; //Коэффициент влажности

        //Временная задача
        protected bool FlagTimeProblem;//
        protected double Deltat;//шаг по времени
        protected double Dissipation; //коэффициент диссипации
        
        //Учёт размеро-зависимости
        protected double l; //Коэффициент связанный с нано
        
        //Разномодульная задача
        protected bool FlagMultiModularIntensiveDeformations;// Учёт экспоненциальной зависимости напряжения от деформаций
        public double SigmS1;  //предел текучести учёт разномодульности
        public double eS1; //деформация текучести учёт разномодульности
        protected bool FlagLinearyAndExponentModelSigmaEpsilon2;// Учёт экспоненциальной зависимости напряжения от деформаций 
        public double[,,] K2 { get; protected set; }// коэффициент объемного сжатия 
        public double Gk2 { get; protected set; }// коэффициент линейности без учёта физ нелинейности сжатия

        //Учёт физической нелинейности
        protected bool FlagFisicalNonLineary;

        //добавление криволиненой поверхности с симетрией относительно срединной плоскости.
        protected bool FlagСurvilinearPlane;
        protected int ModelFunctionСurvilinearPlane;

        //Учёт градиентного распределения материала в доль толщины
        //!!!! Не учитываютася деформации в водимом материале
        protected bool FlagGradientProblem;
        protected double ConstRaspredGradient;
        protected double[,,] GradientE;// Модуль юнга материала используемого в распределеним
        protected double[,,] GradientNU;// Коэффициент пуассона материала используемого в распределеним
        protected double[,,] GradientAlf;// Коэффициент линейнго температурного расширения используемого в распределеним
        protected double[,,] GradientBeta;// Коэффициент для влажности используемого в распределеним

        //учёт пористости
        protected bool FlagPorysotyProblem;
        protected double ConstPoriststb;
        protected double ModelPoriststb;

        public bool Error { get; protected set; } // глобальная ошибка

         //Задача линенйного упрочнение
        protected bool FlagLinearyHard1;
        protected bool FlagLinearyHard2;// Флаг с учёт задачи разномодульности
        public double G01 { get; protected set; }
        public double G11 { get; protected set; }//с учётом разномодульности
         public double G1 { get; protected set; } // модуль сдвига с учётом разномодульности
        //Выводимая ошибка параметров
        //УНИВЕРСАЛЬНЫЕ МЕТОДЫ
        protected void ErrorEnterConstructorParametrs()
        {
            Error = true;
            Console.WriteLine("ErorEnterParameters!!!!!");
        }
        protected void SymmetriInversionMassivCenterPlane(double[,,] Mass)
        {
            double Freedom = 0;
            for (int i = 0; i < Mass.GetLength(0); i++)
                for (int j = 0; j < Mass.GetLength(1); j++)
                    for (int k = 0; k < Mass.GetLength(2) / 2; k++)
                    {
                        Freedom= Mass[i,j, Mass.GetLength(2)-1-k];
                        Mass[i, j, Mass.GetLength(2) - 1 - k] = Mass[i, j, k];
                        Mass[i, j, k] = Freedom;
                    }
        }
        public void ShowMassiv(double[] A)
        {
            Console.WriteLine();
            for (int i = 0; i < A.Count(); i++)
                Console.WriteLine(A[i]);
            Console.WriteLine();
        }
        public void ShowMassiv(double[,] A, int i)
        {
            Console.WriteLine();
            for (int j = 0; j < A.GetLength(1); j++)
                Console.WriteLine("{0:0.00000000} ", A[i, j]);
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
        //проверка на симметрию
        public void Symetri(double[] A, double r)
        {
            bool sim = true;
            for (int i = 0; i < A.Count() / 2; i++)
                if (A[i] - A[A.Count() - 1 - i] > r) { sim = false; Console.WriteLine("d " + i + "  " + (A[i] - A[A.Count() - 1 - i])); }
            Console.WriteLine("sim = " + sim);
        }
        public void Symetri(double[,] A, double r)
        {
            bool sim = true;
            for (int i = 0; i < A.GetLength(0) - 1; i++)
                for (int j = i + 1; j < A.GetLength(0); j++)
                    if (A[i, j] - A[j, i] > r) { sim = false; Console.WriteLine("d " + i + "  " + j + "  " + (A[i, j] - A[j, i])); }
            if (sim)
            {
                for (int i = 0; i < A.GetLength(0) / 2; i++)
                    for (int j = 0; j < A.GetLength(0); j++)
                        if (A[i, j] - A[A.GetLength(0) - 1 - i, j] > r) { sim = false; Console.WriteLine(i + "  " + j + "  " + (A[i, j] - A[A.GetLength(0) - 1 - i, j])); }
            }
            Console.WriteLine("sim = " + sim);
        }
        public void Symetri(double[,,] A, int z, double r)
        {
            bool sim = true;
            for (int i = 0; i < A.GetLength(0) - 1; i++)
                for (int j = i + 1; j < A.GetLength(0); j++)
                    if (A[i, j, z] - A[j, i, z] > r) { sim = false; Console.WriteLine("d " + i + "  " + j + "  " + (A[i, j, z] - A[j, i, z])); }
            if (sim)
            {
                for (int i = 0; i < A.GetLength(0) / 2; i++)
                    for (int j = 0; j < A.GetLength(0); j++)
                        if (A[i, j, z] - A[A.GetLength(0) - 1 - i, j, z] > r) { sim = false; Console.WriteLine(i + "  " + j + "  " + (A[i, j, z] - A[A.GetLength(0) - 1 - i, j, z])); }
            }
            Console.WriteLine("sim = " + sim);
        }
        //Для заполнения массивоф нулями
        protected void FULZero(double[] A)
        {
            for (int i = 0; i < A.Length; i++)
                A[i] = 0;
        }
        protected void FULZero(double[,] A)
        {
            for (int i = 0; i < A.GetLength(0); i++)
                for (int j = 0; j < A.GetLength(1); j++)
                    A[i, j] = 0;
        }
        protected void FULZero(double[,,] A)
        {
            for (int i = 0; i < A.GetLength(0); i++)
                for (int j = 0; j < A.GetLength(1); j++)
                    for (int k = 0; k < A.GetLength(2); k++)
                        A[i, j, k] = 0;
        }
        // Метод Гаусса-Жордана
        protected void MethodGordanGauss(double[,] A)
        {
            //проверка на нулевой массив
            bool Proverka = false;
            for (int i = 0; i < A.GetLength(0); i++) //столбцы
                for (int j = 0; j < A.GetLength(1); j++)
                    if (A[i, j] > 0) Proverka = true;
            if (Proverka)
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
        }

        //ИНИЦИАЛИЗАЦИИ Учёта усложнения модели

        //Инициализация пластины
        public void InicializationPlast(double nu, double G0)
        {
            if (G0 <= 0 || nu < 0)
            {
                ErrorEnterConstructorParametrs();
                return;
            }
            this.G0 = G0;
            this.nu = nu;
            for(int i=0;i<N;i++)
                for(int j=0;j<M;j++)
                    for(int k=0; k<P;k++)
                        this.K1[i,j,k] = 2 * (1 + nu) / (3 * (1 - 2 * nu));

        }
        //Инициализация експоненциальной зависимости напряжения от деформаций
        public void InicializationExpFizicalNonlinProblem(double SigmaS0)
        {
            if (SigmaS0 <= 0)
            {
                ErrorEnterConstructorParametrs();
                return;
            }
            this.SigmS = SigmaS0 * lam1 * lam1 / G0;
            this.eS = SigmS / 3.0;
            FlagLinearyAndExponentModelSigmaEpsilon1 = true;
            Console.WriteLine("OneModuleFizcalExpProblem");
        }
        //Инициализация связанности с температурой(для различных случаев)
        public void InicializationTempature(double QTemp, double alf, double left, double right, double top, double down, double front, double back)
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
            LoadAlfFULL(alf);
            LoadTFULL(QTemp);
            BorderT(left, right, top, down, front, back);
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
            LoadAlfFULL(alf);
            LoadTFULL(QTemp);
        }
        public void InicializationTempature(double QTemp1,double QTemp2, double alf,double gamma)
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
            LoadAlfFULL(alf);
            LoadTFULL(QTemp1,QTemp2,gamma);
        }
        public void InicializationTempature(String PathTemp, double alf)
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

            LoadAlfFULL(alf);
            LoadTFULL(PathTemp, 1);
        } 
        public void InicializationTempatureInversion(String PathTemp, double alf)
        {
            InicializationTempature(PathTemp, alf);
            SymmetriInversionMassivCenterPlane(T);
        }
       
        private double FunFizicalConstWithTemp(double p_1,double p0, double p1, double p2, double p3, double Tem)
        {
            if (Tem == 0)
                return p0;
            else
                return p0 * (p_1 / Tem + 1 + p1 * Tem + p2 * Tem * Tem + p3 * Tem * Tem * Tem);   
        }
        public void InicializationTempatureWithENUG0(double NUp_1, double NUp0, double NUp1, double NUp2, double NUp3, 
            double Ep_1,double Ep0, double Ep1, double Ep2,double Ep3,
            double Alfp_1, double Alfp0, double Alfp1, double Alfp2,double Alfp3)
        {
            // 
            if (FlagСurvilinearPlane)
            {
                ErrorEnterConstructorParametrs();
                Console.WriteLine("Temperature can't job whis CurvilinearPlane");
                return;
            }
            if (FlagTemperatureProblem)
            {

                for (int x = 0; x < N; x++)
                {
                    for (int y = 0; y < M; y++)
                    {
                        for (int z = 0; z < P; z++)
                        {
                            alf[x, y, z] = FunFizicalConstWithTemp(Alfp_1, Alfp0, Alfp1, Alfp2, Alfp3, T[x, y, z]);
                            NU[x, y, z] = FunFizicalConstWithTemp(NUp_1, NUp0, NUp1, NUp2, NUp3, T[x, y, z]);
                            //Безразмерные G0 зависимые от температуры
                            G000[x, y, z] = FunFizicalConstWithTemp(Ep_1, Ep0, Ep1, Ep2, Ep3, T[x, y, z]) / (2 * (1 + NU[x, y, z]) * G0);
                            //безрамерное объёмное рассширение зависимое от температуры
                            K1[x, y, z] = 2 * G000[x, y, z] * (1 + NU[x, y, z]) / (3 * (1 - 2 * NU[x, y, z]));
                            SigmaT1[x, y, z] = 3 * (SigmS - SigmS * (1 - G000[x, y, z])) / (2 * (1 + NU[x, y, z]));
                            if (FlagMultiModularIntensiveDeformations)
                            {
                                G001[x, y, z] = FunFizicalConstWithTemp(Ep_1, Ep0, Ep1, Ep2, Ep3, T[x, y, z]) / (2 * (1 + NU[x, y, z]) * G0);
                                //безрамерное объёмное рассширение зависимое от температуры
                                K2[x, y, z] = 2 * G001[x, y, z] * (1 + NU[x, y, z]) / (3 * (1 - 2 * NU[x, y, z]));
                                SigmaT2[x, y, z] = 3 * (SigmS1 - SigmS1 * (1 - G001[x, y, z])) / (2 * (1 + NU[x, y, z]));
                            }
                        }
                    }
                }
                this.FlagTemperatureWithENUG0 = true; 
            }
            else
            {
                Console.WriteLine("Не временная задача");
                return;
            }
        }
        
        public void InicializationTempatureWithDefform(String PathTemp, String PathDeformMap, double alf, double Tem0, double Tem1, double Def0, double Def1, double a)
        {
            // 
            if (FlagСurvilinearPlane)
            {
                ErrorEnterConstructorParametrs();
                Console.WriteLine("Temperature can't job whis CurvilinearPlane so");
                return;
            }
            if (alf <= 0 || alf > 1)
            {
                ErrorEnterConstructorParametrs();
                return;
            }
            this.FlagDeformMapT = true;
            this.FlagTemperatureProblem = true;
            LoadAlfFULL(alf);
            double po = alf * lam1 * lam1;
            //значения вводятся в Кельвинах
            this.Temp0 = Tem0 * po;
            this.Temp1 = Tem1 * po;
            //значения вводятся не безразмерными
            this.ei0 = Def0 * lam1 * lam1;
            this.ei1 = Def1 * lam1 * lam1;
            LoadTFULL(PathTemp, a);
            LoadDeformMapT(PathDeformMap);


        }
        public void InicializationHumidity(double QWater1, double QWater2, double gamma, double bet)
        {
            if (FlagСurvilinearPlane)
            {
                ErrorEnterConstructorParametrs();
                Console.WriteLine("Temperature can't job whis CurvilinearPlane");
                return;
            }
            if (bet <= 0 || bet > 1 )
            {
                ErrorEnterConstructorParametrs();
                return;
            }

            this.FlagHumidity = true;
            if (!FlagPorysotyProblem && !FlagGradientProblem)
                LoadBetaFULL(bet);
            
            LoadCFULL(QWater1, QWater2,gamma);
            BezrazmernostbC();
        }
        public void InicializationHumidity(String PathC, double bet)
        {
            if (FlagСurvilinearPlane)
            {
                ErrorEnterConstructorParametrs();
                Console.WriteLine("Temperature can't job whis CurvilinearPlane");
                return;
            }
            if (bet <= 0 || bet > 1)
            {
                ErrorEnterConstructorParametrs();
                return;
            }

            this.FlagHumidity = true;
            if (!FlagPorysotyProblem && !FlagGradientProblem)
                LoadBetaFULL(bet);

            LoadCFULL(PathC,0.1);
            //SymmetriInversionMassivCenterPlane(C);
            BezrazmernostbC();
        }
        //Инициализация градиентного добавочного материала
        //влажность локальная_400_21x21_30x30x19_300-400_(20)
        public void InicializationVariationGradient(double ConstRaspredGradient, double NU1, double E1)
        {
            if (ConstRaspredGradient <= 0 || FlagTemperatureProblem)
            {
                ErrorEnterConstructorParametrs();
                return;
            }
            FlagGradientProblem = true;
            this.ConstRaspredGradient = ConstRaspredGradient;
            for (int x = 0; x < N; x++)
            {
                for (int y = 0; y < M; y++)
                {
                    for (int z = 0; z < P; z++)
                    {
                        GradientE[x, y, z] = E1 / G0;
                        GradientNU[x, y, z] = NU1;                        
                    }
                }
            }
        }
        public void InicializationVariationGradient(double ConstRaspredGradient, double NU1, double E1, double alf1)
        {
            if (!FlagTemperatureProblem)
            {
                ErrorEnterConstructorParametrs();
                return;
            }
            if (ConstRaspredGradient <= 0 || alf1 <= 0 || alf1 > 1)
            {
                ErrorEnterConstructorParametrs();
                return;
            }
            this.ConstRaspredGradient = ConstRaspredGradient;
            FlagGradientProblem = true;
            for (int x = 0; x < N; x++)
            {
                for (int y = 0; y < M; y++)
                {
                    for (int z = 0; z < P; z++)
                    {
                        GradientE[x, y, z] = E1 / G0;
                        GradientNU[x, y, z] = NU1;
                        GradientAlf[x, y, z] = alf1;
                    }
                }
            }
            Gradient(alf, GradientAlf);
            //Необходимо дописать температуру

        }
        public void InicializationVariationGradient(double ConstRaspredGradient,
            double p_1NU1, double p0NU1, double p1NU1, double p2NU1, double p3NU1,
            double p_1E1,double p0E1, double p1E1, double p2E1, double p3E1,
            double p_1alf1,double p0alf1, double p1alf1, double p2alf1, double p3alf1)
        {
            if (!FlagTemperatureProblem)
            {
                ErrorEnterConstructorParametrs();
                return;
            }
            if (ConstRaspredGradient <= 0)
            {
                ErrorEnterConstructorParametrs();
                return;
            }
            FlagGradientProblem = true;
            this.ConstRaspredGradient = ConstRaspredGradient;
            for (int x = 0; x < N; x++)
            {
                for (int y = 0; y < M; y++)
                {
                    for (int z = 0; z < P; z++)
                    {
                        GradientAlf[x, y, z] = FunFizicalConstWithTemp(p_1alf1, p0alf1, p1alf1, p2alf1, p3alf1, T[x, y, z]);
                        GradientNU[x, y, z] = FunFizicalConstWithTemp(p_1NU1, p0NU1, p1NU1, p2NU1, p3NU1, T[x, y, z]);
                        GradientE[x, y, z] = FunFizicalConstWithTemp(p_1E1, p0E1, p1E1, p2E1, p3E1, T[x, y, z])/ G0;
                    }
                }
            }
            Gradient(alf, GradientAlf);

        }

        //Инициализация пористости
        public void InicializationPorisoty(double CountPorisoty, int ModelPorisoty)
        {
            if (CountPorisoty <= 0 || CountPorisoty > 1 || ModelPorisoty < 1 || ModelPorisoty > 4 )
            {
                ErrorEnterConstructorParametrs();
                return;
            }
            FlagPorysotyProblem = true;
            this.ConstPoriststb = CountPorisoty;
            this.ModelPoriststb = ModelPorisoty;
            if (FlagTemperatureProblem) Porisoty(alf);
            if (FlagHumidity) Porisoty(beta);
        }
        //Инициализация криволинейнойной поверхности на гранях симметричных относительно срединной плоскости.
        public void InicializationCurvelinearyPlane(int ModelCurvilinearPlane)
        {

            if (ModelCurvilinearPlane <= 0 || ModelCurvilinearPlane > 3)
            {
                ErrorEnterConstructorParametrs();
                return;
            }
            this.ModelFunctionСurvilinearPlane = ModelCurvilinearPlane;
            this.FlagСurvilinearPlane = true;
            LoadDzForСurvilinearPlane();


        }
        //Инициализация линейного упрочнения
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
        //Инициализация мульти модульной задачи
        public void InicializationMultimodulProblem(double nu2, double G1)
        {
            if (nu2 <= 0 || G1 <= 0)
            {
                ErrorEnterConstructorParametrs();
                return;
            }
            this.FlagMultiModularIntensiveDeformations = true;
            this.G1 = G1;
            this.Gk2 = G1 / G0;
            for (int i = 0; i < N; i++)
                for (int j = 0; j < M; j++)
                    for (int k = 0; k < P; k++)
                        this.K2[i, j, k] = Gk2 * 2 * (1 + nu2) / (3 * (1 - 2 * nu2));


        }
        //инициализация мульти модульной задачи при экспоненциальной зависимости напряжения от деформации
        public void InicializationMultimoduWithExpFizicalNonlinlProblem(double SigmaS1)
        {
            if ( SigmaS1 <= 0 || !FlagMultiModularIntensiveDeformations)
            {
                ErrorEnterConstructorParametrs();
                return;
            }
            this.FlagLinearyAndExponentModelSigmaEpsilon2 = true;
            this.FlagMultiModularIntensiveDeformations = true;
            this.SigmS1 = SigmaS1 * lam1 * lam1 / G0;
            this.eS1 = (SigmS1 / 3.0)*G0/G1;
            Console.WriteLine("MultiModuleFizcalExpProblem");
        }
        //инициализация нано размерности
        public void InicializationNano(double l)
        {
            if (l <= 0 || l > 1)
            {
                ErrorEnterConstructorParametrs();
                return;
            }
            this.l = l;
        }

        //ПЕРЕОПРЕДЕЛЯЮЩИЕСЯ МЕТОДЫ
        //Вычисление интенсивности деформаций
        protected virtual double ei(int x, int y, int z)
        {
            return 0;
        }
        public virtual void StaticDecisionMethodVariationIteration(Load F, int BorderExitIterationMethods,
    int BorderExitMethodABS, int BorderExitPhisicalIteration, int TypeBorder1, int TypeBorder2, int TypeBorder3, int TypeBorder4)
        {
        }
        public virtual void StaticDecisionMethodKantorovichVlasovOneLevel(Load F, int BorderExitIterationMethods,
            int BorderExitMethodABS, int BorderExitPhisicalIteration, int TypeBorder1, int TypeBorder2, int TypeBorder3, int TypeBorder4)
        {
        }
        public virtual void StaticDecisionMethodVariationIterartionTwoLevel(Load F, int BorderExitIterationMethods,
            int BorderExitMethodABS, int BorderExitPhisicalIteration, int TypeBorder1, int TypeBorder2, int TypeBorder3, int TypeBorder4)
        {
        }
        public virtual void StaticDecisionMethodKantorovichVlasovTwoLevel(Load F, int BorderExitIterationMethods,
            int BorderExitMethodABS, int BorderExitPhisicalIteration, int TypeBorder1, int TypeBorder2, int TypeBorder3, int TypeBorder4)
        {
        }
        public virtual void StaticDecisionMethodBubnovGalercin(Load F, int DegreeOfDifficFunctionMethodBubnovaGalercin,
            int BorderExitMethodABS, int BorderExitPhisicalIteration, int TypeBorder1, int TypeBorder2, int TypeBorder3, int TypeBorder4)
        {
        }
        public virtual void DinamiDecisionMethodRugeKyta4(Load F, double t0, double t1, double Deltat, int BorderExitPhisicalIteration,
            int TypeBorder1, int TypeBorder2, int TypeBorder3, int TypeBorder4, int WriteShag)
        {
        }
        public virtual void DinamiDecisionMethodVariationIteration(Load F, double t0, double t1, double Deltat, int BorderExitPhisicalIteration,
            int TypeBorder1, int TypeBorder2, int TypeBorder3, int TypeBorder4, int WriteShag)
        {
        }


        //Базовый конструктор
        protected Plane(int N, int M, int P, double n, double m, double p)
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
            W = new double[N + 2, M + 2];
            PrintW = new double[N, M];
            Eii = new double[N, M, this.P];
            T = new double[N, M, this.P];
            Temp = new double[N, M, this.P];
            DZ = new double[N, M];
            alf = new double[N, M, this.P];
            C = new double[N, M, this.P];
            beta = new double[N, M, this.P];
            K1 = new double[N, M, this.P];
            K2 = new double[N, M, this.P];
            SigmaT1 = new double[N, M, this.P];
            SigmaT2 = new double[N, M, this.P];
            G000 = new double[N, M, this.P];
            G001 = new double[N, M, this.P];


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
            this.Gk2 = 0;
            // учёт нано
            this.l = 0;

            //cвязанность с температурой
            this.FlagDeformMapT = false;
            this.FlagTemperatureProblem = false;
            this.FlagTemperatureWithENUG0 = false;

            this.Temp0 = 0;
            this.Temp1 = 0;
            this.ei0 = 0;
            this.ei1 = 0;
            //учёты распределения градиента
            this.FlagGradientProblem = false;
            this.ConstRaspredGradient = 0;
            this.GradientE = new double [N,M,this.P];
            this.GradientNU = new double[N, M, this.P];
            this.GradientAlf = new double[N, M, this.P];

            mW = 0; //максимальный прогиб


            //Временной параметр
            this.FlagTimeProblem = false;
            this.Dissipation = 0;
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

        //базовое заполнение распределение разбиения по оси z (введено для криволиненых поверхностей)
        protected void LoadDZ(double Dz)
        {
            for (int x = 0; x < N; x++)
                for (int y = 0; y < M; y++)
                    DZ[x, y] = Dz;
        }
       
        // зависимости напряжений от деформаций
        protected double SigmaEii(double eii)
        {
            double ans = 0;
            if (eii < 0)
            {
                if (FlagLinearyAndExponentModelSigmaEpsilon2)
                    ans = (double)(-SigmS1 * (1 - Math.Exp(eii / eS1)));
                else
                    ans = (double)(3 * Gk2 * eii);
            }
            else
            {
                if (FlagLinearyAndExponentModelSigmaEpsilon1)
                { 
                    ans = (double)(SigmS * (1 - Math.Exp(-eii / eS)));
                }
                else
                {
                    if (FlagLinearyHard1)
                    {
                        if (eii <= eS)
                            ans = 3 * eii;
                        else
                            ans = (eS * (G0 - G01) + G01 * eii) * 3 / G0;
                    }
                    else
                        ans = (double)(3 * eii);
                }
            }
            return ans;
        }
        protected double SigmaEii(double eii,int x,int y,int z)
        {
            if(eii>=0)
                return SigmaT1[x, y, z] * (1 - Math.Exp(-eii / ((SigmaT1[x, y, z]) /( 3*G000[x,y,z]))));
            else
                return -SigmaT2[x, y, z] * (1 - Math.Exp(eii / ((SigmaT2[x, y, z]) / (3 * G001[x, y, z]))));
        }

        // функция для оперделения пористости
        protected double PoristFunctionToHeigth(double x, double y, double z)
        {
            switch (ModelPoriststb)
            {
                //U-P
                case 1:
                    return 1 - ConstPoriststb;
                case 2:
                    if (FlagСurvilinearPlane) return 1 - (Math.Abs(z)) * ConstPoriststb / (1 + FunctionСurvilinearPlane(x, y));
                    return 1 - ConstPoriststb * Math.Abs(z);
                case 3:
                    if (FlagСurvilinearPlane) return 1 - ((1 + FunctionСurvilinearPlane(x, y)) - Math.Abs(z)) * ConstPoriststb / (1 + FunctionСurvilinearPlane(x, y));
                    else return 1 - ConstPoriststb * (1 - Math.Abs(z));
                case 4:
                    return 1 - ConstPoriststb * Math.Cos(z * Math.PI / 2);
                default:
                    return 0;
            }
        }
        //учёт пористости
        protected void Porisoty(double[,,] Mas)
        {
            for (int x = 0; x < N; x++)
                for (int y = 0; y < M; y++)
                    for (int z = 0; z < P; z++)
                        Mas[x, y, z] = Mas[x, y, z] * PoristFunctionToHeigth(dx * x, dy * y, DZ[x, y] * (z - (P - 1) / 2));
        }
        protected double FunctionGradient(int x,int y,int z,double Constgradients)
        {
            if (z == 0) return 0;
            return Math.Pow(1.0 / 2 + (DZ[x, y] * (z - ((P - 1) / 2)) / (2 * (1 + FunctionСurvilinearPlane(dx * x, dy * y)))), Constgradients);
        }
        protected double FunctionGradientT(int x, int y, int z, double Constgradients)
        {
            if (z == 0) return 0;
            return Math.Pow(1.0 / 2 + ((2 * DZ[x, y] * (z - ((P - 1) / 2)) + 1) / (2 * (1 + FunctionСurvilinearPlane(dx * x, dy * y)))), Constgradients);
        }
        //учёт градиента
        protected void Gradient(double[,,] Mas, double[,,] GradientMaterial)
        {
            for (int x = 0; x < N; x++)
                for (int y = 0; y < M; y++)
                    for (int z = 0; z < P; z++)
                    {
                        if (z == 0) Mas[x, y, z] = GradientMaterial[x, y, z];
                        else if (z == P - 1) Mas[x, y, z] = Mas[x, y, z];
                        else Mas[x, y, z] = GradientMaterial[x, y, z] + (Mas[x, y, z] - GradientMaterial[x, y, z]) * FunctionGradient(x, y, z,ConstRaspredGradient);
                    }
                        
        }

        // Расчёт распределения температуры при её учёте в зависимости
        //от граничных условий с точностью до 6 знака
        public void BezrazmernostbT()
        {
            double MinTemp = T[0,0,0];
            for (int x = 0; x < N; x++)
                for (int y = 0; y < M; y++)
                    for (int z = 0; z < P; z++)
                        if (MinTemp > T[x, y, z]) MinTemp = T[x, y, z];
            for (int x = 0; x < N; x++)
                for (int y = 0; y < M; y++)
                    for (int z = 0; z < P; z++)
                    {
                        //T[x, y, z] = T[x, y, z] * alf[x, y, z] * lam1 * lam1;
                        T[x, y, z] = (T[x, y, z] - MinTemp) * alf[x, y, z] * lam1 * lam1;
                    }
        }
        public void BezrazmernostbC()
        {
            for (int x = 0; x < N; x++)
                for (int y = 0; y < M; y++)
                    for (int z = 0; z < P; z++)
                        C[x, y, z] = C[x, y, z] * beta[x, y, z] * lam1 * lam1;
        }
        private void LoadDeformMapT(String PathDeformMap)
        {
            double r = lam1 * lam1 / G0;

            System.Threading.Thread.CurrentThread.CurrentCulture = new System.Globalization.CultureInfo("en-US");// смена точки на запятую

            //Считываем данные с распределение напряжений
            FileStream file2 = new FileStream(PathDeformMap, FileMode.Open); //создаем файловый поток
            StreamReader reader2 = new StreamReader(file2);//создаем «потоковый писатель» и связываем его с файловым потоком 
                                                           //определение размерности массива зависимости
            int e = 0;
            int t = 0;
            string s = "";
            s = "" + reader2.ReadLine();
            foreach (var number in s.Split())
                if (number != "")
                    t++;
            while (true)
            {
                e++;
                s = "" + reader2.ReadLine();
                if (s == "") break;
            }
            DeformMapT = new double[t + 2, e + 2];
            //возвращение чтения в начальную позицию
            reader2.BaseStream.Position = 0;
            //Забитие безразмерных данных в массив


            for (int i = 0; i < e; i++)
            {
                int j = 0;
                s = "" + reader2.ReadLine();
                j = 0;
                foreach (var number in s.Split())
                    if (number != "")
                    {
                        DeformMapT[j + 1, i + 1] = Convert.ToDouble(number) * r;
                        j++;
                    }
            }
            //Установление за граничных значений по линейной зависимости
            for (int i = 0; i < e; i++)
            {
                DeformMapT[0, i + 1] = 2 * DeformMapT[1, i + 1] - DeformMapT[2, i + 1];
                DeformMapT[t + 1, i + 1] = 2 * DeformMapT[t, i + 1] - DeformMapT[t - 1, i + 1];
            }
            for (int j = 0; j < t; j++)
            {
                DeformMapT[j + 1, 0] = 2 * DeformMapT[j + 1, 1] - DeformMapT[j + 1, 2];
                DeformMapT[j + 1, e + 1] = 2 * DeformMapT[j + 1, e] - DeformMapT[j + 1, e - 1];
            }
            //угловые точки
            DeformMapT[0, 0] = DeformMapT[0, 1] - 0.5 * DeformMapT[0, 2] + DeformMapT[1, 0] - 0.5 * DeformMapT[2, 0];
            DeformMapT[0, e + 1] = DeformMapT[0, e] - 0.5 * DeformMapT[0, e - 1] + DeformMapT[1, e + 1] - 0.5 * DeformMapT[2, e + 1];
            DeformMapT[t + 1, 0] = DeformMapT[t + 1, 1] - 0.5 * DeformMapT[t + 1, 2] + DeformMapT[t, 0] - 0.5 * DeformMapT[t - 1, 0];
            DeformMapT[t + 1, e + 1] = DeformMapT[t + 1, e] - 0.5 * DeformMapT[t + 1, e - 1] + DeformMapT[t, e + 1] - 0.5 * DeformMapT[t - 1, e + 1];

            reader2.Close(); //закрываем поток. Не закрыв поток, в файл ничего не запишется так как первичная запись идёт в ОЗУ
        }
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
        private void LoadFileMassiv(String PathTemp,double[,,] Massive)
        {
            System.Threading.Thread.CurrentThread.CurrentCulture = new System.Globalization.CultureInfo("en-US");// смена точки на запятую
            FileStream file1 = new FileStream(PathTemp, FileMode.Open); //создаем файловый поток
            StreamReader reader = new StreamReader(file1);//создаем «потоковый писатель» и связываем его с файловым потоком 
            string s = "";
            while (s != "% Data")
                s = "" + reader.ReadLine();
            int x;
            for (int z = 0; z < P; z++)

                for (int y = 0; y < M; y++)
                {
                    s = "" + reader.ReadLine();
                    x = 0;
                    foreach (var number in s.Split())
                    {
                        if (number != "")
                        {
                            Massive[x, y, z] = Convert.ToDouble(number);
                            x++;
                        }
                    }
                }
            reader.Close(); //закрываем поток. Не закрыв поток, в файл ничего не запишется так как первичная запись идёт в ОЗУ
        }

        private void LoadTFULL(double ConstTemp)
        {
            for (int x = 0; x < N; x++)
                for (int y = 0; y < M; y++)
                    for (int z = 0; z < P; z++)
                        T[x, y, z] = ConstTemp ;
        }
        private void LoadTFULL(String PathTemp, double a)
        {
            LoadFileMassiv(PathTemp, T);
            for (int i = 0; i < T.GetLength(0); i++)
            {
                for (int j = 0; j < T.GetLength(1); j++)
                {
                    for (int k = 0; k < T.GetLength(2); k++)
                    {
                        T[i, j, k] = T[i, j, k] * a;
                    }
                }
            }
        }
        private void LoadTFULL(double ConstTemp1, double ConstTemp2, double gamma)
        {
            for (int x = 0; x < N; x++)
                for (int y = 0; y < M; y++)
                    for (int z = 0; z < P; z++)
                    {
                        T[x, y, z] = ConstTemp1 + (ConstTemp2 - ConstTemp1) * FunctionGradientT(x, y, z, gamma);
                        //if (z < P / 2 + 1)
                        //    T[x, y, z] = ConstTemp1 + (ConstTemp2 - ConstTemp1) * FunctionGradientT(x, y, z, gamma);
                        //else
                        //    T[x, y, z] = T[x, y, P - 1 - z];
                    }
        }

        private void LoadAlfFULL(double ConstTemp)
        {
            for (int x = 0; x < N; x++)
                for (int y = 0; y < M; y++)
                    for (int z = 0; z < P; z++)
                        alf[x, y, z] = ConstTemp;
        }
        //граничные условия распределение температуры (учитываются в LoadT())
        private void BorderT(double left, double right, double top, double down, double front, double back)
        {
            double r = lam1 * lam1;
            //левая грань
            for (int y = 0; y < M; y++)
                for (int z = 0; z < P; z++)
                    T[0, y, z] = left * alf[0, y, z] * r;
            //правая грань
            for (int y = 0; y < M; y++)
                for (int z = 0; z < P; z++)
                    T[N - 1, y, z] = right * alf[N - 1, y, z] * r;
            //верхняя грань
            for (int y = 0; y < M; y++)
                for (int x = 0; x < N; x++)
                    T[x, y, P - 1] = top * alf[x, y, P - 1] * r;
            //нижняя грань
            for (int y = 0; y < M; y++)
                for (int x = 0; x < N; x++)
                    T[x, y, 0] = down * alf[x, y, 0] * r;
            //фронт 
            for (int x = 0; x < N; x++)
                for (int z = 0; z < P; z++)
                    T[x, 0, z] = front * alf[x, 0, z] * r;
            //задняя грань
            for (int x = 0; x < N; x++)
                for (int z = 0; z < P; z++)
                    T[x, M - 1, z] = back * alf[x, M - 1, z] * r;
        }
        // уравнение для расчёта точек температуры 7 точечного шаблона
        private double YravnenieForT(double[,,] A, int x, int y, int z)
        {
            double del = 2 / (dx * dx) + 2 * lam * lam / (dy * dy) + 2 * lam1 * lam1 / (dz * dz);
            return ((A[x + 1, y, z] + A[x - 1, y, z]) / (dx * dx) +
                lam * lam * (A[x, y + 1, z] + A[x, y - 1, z]) / (dy * dy) +
                lam1 * lam1 * (A[x, y, z + 1] + A[x, y, z - 1]) / (dz * dz) + Temp[x, y, z]) / del;
        }
        private void LoadCFULL(double ConstTemp)
        {
            for (int x = 0; x < N; x++)
                for (int y = 0; y < M; y++)
                    for (int z = 0; z < P; z++)
                        C[x, y, z] = ConstTemp;
        }
        private void LoadCFULL(double ConstCemp1, double ConstCemp2, double gamma)
        {
            for (int x = 0; x < N; x++)
                for (int y = 0; y < M; y++)
                    for (int z = 0; z < P; z++)
                        if (z < P / 2 + 1)
                            C[x, y, z] = ConstCemp1 + (ConstCemp2 - ConstCemp1) * FunctionGradientT(x, y, z, gamma);
                        else
                            C[x, y, z] = C[x, y, P - 1 - z];
        }
        private void LoadCFULL(String PathC, double a)
        {
            LoadFileMassiv(PathC, C);
            for (int i = 0; i < C.GetLength(0); i++)
            {
                for (int j = 0; j < C.GetLength(1); j++)
                {
                    for (int k = 0; k < C.GetLength(2); k++)
                    {
                        C[i, j, k] = C[i, j, k] * a;
                    }
                }
            }
        }
        private void LoadBetaFULL(double ConstTemp)
        {
            for (int x = 0; x < N; x++)
                for (int y = 0; y < M; y++)
                    for (int z = 0; z < P; z++)
                    {
                        beta[x, y, z] = ConstTemp;
                    }
        }

        protected double FunctionСurvilinearPlane(double x, double y)
        {
            //максимум 0.125
            switch (ModelFunctionСurvilinearPlane)
            {
                case 1:
                    return -(Math.Pow(Math.Cos(x * Math.PI) * Math.Cos(y * Math.PI), 2) - 1) / 8;
                case 2:
                    return Math.Pow(Math.Cos(x * Math.PI) * Math.Cos(y * Math.PI), 2) / 8;
                case 3:
                    return -(Math.Pow(Math.Sin((x - y) * Math.PI) * Math.Sin((x + y) * Math.PI), 2) - 1) / 8;
                default:
                    return 0;

            }
        }
       
        protected void LoadDzForСurvilinearPlane()
        {
            for (int x = 0; x < N; x++)
                for (int y = 0; y < M; y++)
                    DZ[x, y] = dz * (FunctionСurvilinearPlane(dx * x, dy * y) + 1);
        }
        protected double SigmaEiiTIntorpaletion(double eii, double Tem)
        {
            //установка минимальной и максимальной температуры в файле распределения
            //температура в кельвинах
            double x = Tem;
            double y = eii;

            //точки определяющие целочисленное положение на массиве s(e,t)
            double DeltTemp = (Temp1 - Temp0) / (DeformMapT.GetLength(0) - 3);
            double DeltEi = (ei1 - ei0) / (DeformMapT.GetLength(1) - 3);

            //ццент относительно которого считается кубическая интерполяция
            int Centri0 = (int)Math.Truncate((Tem - Temp0) / DeltTemp) + 1;
            int Centrj0 = (int)Math.Truncate((eii - ei0) / DeltEi) + 1;

            if (Centri0 > DeformMapT.GetLength(0) - 3) Centri0 = DeformMapT.GetLength(0) - 3;
            if (Centrj0 > DeformMapT.GetLength(1) - 3) Centrj0 = DeformMapT.GetLength(1) - 3;
            //условие для отрицатеьных деформаций и температур
            if (Centri0 < 1) Centri0 = 1;
            if (Centrj0 < 1) Centrj0 = 1;
            //модификация eii и Т для вычисления бикубической инеропляции
            x = (x - (Centri0 - 1) * DeltTemp - Temp0) / DeltTemp;
            y = (y - (Centrj0 - 1) * DeltEi) / DeltEi;
            //ограничение значений на вылет из массива

            if (x > 1) x = 1;
            if (x < 0) x = 0;
            if (y > 1) y = 1;

            //Значения точек в массивах
            double f00 = DeformMapT[Centri0 - 1, Centrj0 - 1];
            double f01 = DeformMapT[Centri0 - 1, Centrj0];
            double f10 = DeformMapT[Centri0, Centrj0 - 1];
            double f11 = DeformMapT[Centri0, Centrj0];
            double f02 = DeformMapT[Centri0 - 1, Centrj0 + 1];
            double f20 = DeformMapT[Centri0 + 1, Centrj0 - 1];
            double f12 = DeformMapT[Centri0, Centrj0 + 1];
            double f21 = DeformMapT[Centri0 + 1, Centrj0];
            double f22 = DeformMapT[Centri0 + 1, Centrj0 + 1];
            double f03 = DeformMapT[Centri0 - 1, Centrj0 + 2];
            double f30 = DeformMapT[Centri0 + 2, Centrj0 - 1];
            double f13 = DeformMapT[Centri0, Centrj0 + 2];
            double f31 = DeformMapT[Centri0 + 2, Centrj0];
            double f23 = DeformMapT[Centri0 + 1, Centrj0 + 2];
            double f32 = DeformMapT[Centri0 + 2, Centrj0 + 1];
            double f33 = DeformMapT[Centri0 + 2, Centrj0 + 2];

            //члены бикубической интерполяции
            double a00 = f11;
            double a01 = -0.5 * f10 + 0.5 * f12;
            double a02 = f10 - 2.5 * f11 + 2 * f12 - 0.5 * f13;
            double a03 = -0.5 * f10 + 1.5 * f11 - 1.5 * f12 + 0.5 * f13;
            double a10 = -0.5 * f01 + 0.5 * f21;
            double a11 = 0.25 * f00 - 0.25 * f02 - 0.25 * f20 + 0.25 * f22;
            double a12 = -0.5 * f00 + 1.25 * f01 - f02 + 0.25 * f03 + 0.5 * f20 - 1.25 * f21 + f22 - 0.25 * f23;
            double a13 = 0.25 * f00 - 0.75 * f01 + 0.75 * f02 - 0.25 * f03 - 0.25 * f20 + 0.75 * f21 - 0.75 * f22 + 0.25 * f23;
            double a20 = f01 - 2.5 * f11 + 2 * f21 - 0.5 * f31;
            double a21 = -0.5 * f00 + 0.5 * f02 + 1.25 * f10 - 1.25 * f12 - f20 + f22 + 0.25 * f30 - 0.25 * f32;
            double a22 = f00 - 2.5 * f01 + 2 * f02 - 0.5 * f03 - 2.5 * f10 + 6.25 * f11 - 5 * f12 + 1.25 * f13
                + 2 * f20 - 5 * f21 + 4 * f22 - f23 - 0.5 * f30 + 1.25 * f31 - f32 + 0.25 * f33;
            double a23 = -0.5 * f00 + 1.5 * f01 - 1.5 * f02 + 0.5 * f03 + 1.25 * f10 - 3.75 * f11 + 3.75 * f12
                - 1.25 * f13 - f20 + 3 * f21 - 3 * f22 + f23 + 0.25 * f30 - 0.75 * f31 + 0.75 * f32 - 0.25 * f33;
            double a30 = -0.5 * f01 + 1.5 * f11 - 1.5 * f21 + 0.5 * f31;
            double a31 = 0.25 * f00 - 0.25 * f02 - 0.75 * f10 + 0.75 * f12 + 0.75 * f20 - 0.75 * f22 - 0.25 * f30 + 0.25 * f32;
            double a32 = -0.5 * f00 + 1.25 * f01 - f02 + 0.25 * f03 + 1.5 * f10 - 3.75 * f11 + 3 * f12
                - 0.75 * f13 - 1.5 * f20 + 3.75 * f21 - 3 * f22 + 0.75 * f23 + 0.5 * f30 - 1.25 * f31 + f32 - 0.25 * f33;
            double a33 = 0.25 * f00 - 0.75 * f01 + 0.75 * f02 - 0.25 * f03 - 0.75 * f10 + 2.25 * f11 - 2.25 * f12 + 0.75 * f13
                + 0.75 * f20 - 2.25 * f21 + 2.25 * f22 - 0.75 * f23 - 0.25 * f30 + 0.75 * f31 - 0.75 * f32 + 0.25 * f33;
            double sigmaI = (a00 + a01 * y + a02 * y * y + a03 * y * y * y) +
                            (a10 + a11 * y + a12 * y * y + a13 * y * y * y) * x +
                            (a20 + a21 * y + a22 * y * y + a23 * y * y * y) * x * x +
                            (a30 + a31 * y + a32 * y * y + a33 * y * y * y) * x * x * x;

            return sigmaI;

        }
        // подсчёт G
        protected double G(double eii)
        {
            if (eii == 0)
                return 1;
            else
                return SigmaEii(eii) / (3 * eii);
        }
        protected double G(double eii, double Tempature)
        {
            if (eii == 0)
                return 1;
            else
                return SigmaEiiTIntorpaletion(eii, Tempature) / (3 * eii);
        }
        protected double G(double eii, int x,int y,int z)
        {
            if (eii == 0)
                return G000[x,y,z];
            else
                return SigmaEii(eii, x,y,z) / (3 * eii);
        }

        protected void LoadEii()
        {
            for (int y = 0; y < M; y++)
            {
                for (int x = 0; x < N; x++)
                {
                    for (int z = 0; z < P; z++)
                        Eii[x, y, z] = ei(x, y, z);
                }
            }
        }
        // расчёт массива коэффициентов Пуссона
        protected void LoadNU()
        {
            double Gg = 0;
            for (int x = 0; x < N; x++)
            {
                for (int y = 0; y < M; y++)
                {
                    for (int z = 0; z < P; z++)
                    {
                        Gg = G(Eii[x, y, z]);
                        if (FlagDeformMapT)
                            Gg = G(Eii[x, y, z], T[x, y, z]);
                        if (FlagTemperatureWithENUG0)
                            Gg = G(Eii[x, y, z], x, y, z);
                        if (Eii[x, y, z] < 0)
                            NU[x, y, z] = (3 * K2[x,y,z] - 2 * Gg) / (6 * K2[x, y, z] + 2 * Gg);
                        else
                            NU[x, y, z] = (3 * K1[x, y, z] - 2 * Gg) / (6 * K1[x, y, z] + 2 * Gg);
                    }
                }
            }
        }
        // расчёт массива Модулей Юнга
        protected void LoadE()
        {
            for (int x = 0; x < N; x++)
            {
                for (int y = 0; y < M; y++)
                {
                    for (int z = 0; z < P; z++)
                    {
                        double Gg = G(Eii[x, y, z]);
                        if (FlagDeformMapT)
                            Gg = G(Eii[x, y, z], T[x, y, z]);
                        if (FlagTemperatureWithENUG0)
                            Gg = G(Eii[x, y, z], x, y, z);
                        if (Eii[x, y, z] < 0)
                            E[x, y, z] = 9 * K2[x, y, z] * Gg / (3 * K2[x, y, z] + Gg);
                        else
                            E[x, y, z] = 9 * K1[x, y, z] * Gg / (3 * K1[x, y, z] + Gg);
                    }
                }
            }
        }
    }
}
