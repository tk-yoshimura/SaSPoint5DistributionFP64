// Author and Approximation Formula Coefficient Generator: T.Yoshimura
// Github: https://github.com/tk-yoshimura
// Original Code: https://github.com/tk-yoshimura/SaSPoint5DistributionFP64

using SaSPoint5DistributionFP64.InternalUtils;
using SaSPoint5DistributionFP64.RandomGeneration;
using System.Collections.ObjectModel;
using System.Diagnostics;
using static System.Double;

namespace SaSPoint5DistributionFP64 {
    [DebuggerDisplay("{ToString(),nq}")]
    public class SaSPoint5Distribution {

        public double Mu { get; }

        public double C { get; }

        private readonly double c_inv;

        private static readonly double entropy_base = 3.63992444568030649573;


        public SaSPoint5Distribution() : this(mu: 0d, c: 1d) { }

        public SaSPoint5Distribution(double c) : this(mu: 0d, c: c) { }

        public SaSPoint5Distribution(double mu, double c) {
            if (!IsFinite(mu)) {
                throw new ArgumentOutOfRangeException(nameof(mu), "Invalid location parameter.");
            }
            if (!(c > 0 && IsFinite(c))) {
                throw new ArgumentOutOfRangeException(nameof(c), "Invalid scale parameter.");
            }

            Mu = mu;
            C = c;

            c_inv = 1d / c;
        }

        public double PDF(double x) {
            double u = (x - Mu) * c_inv;

            if (IsNaN(u)) {
                return NaN;
            }
            if (IsInfinity(u)) {
                return 0d;
            }

            double pdf = PDFPade.Value(u) * c_inv;

            return pdf;
        }

        public double CDF(double x, Interval interval = Interval.Lower) {
            double u = (x - Mu) * c_inv;

            if (IsNaN(u)) {
                return NaN;
            }

            double cdf = (interval == Interval.Lower) ? CDFPade.Value(-u) : CDFPade.Value(u);

            return cdf;
        }

        public double Quantile(double p, Interval interval = Interval.Lower) {
            if (!(p >= 0d && p <= 1d)) {
                return NaN;
            }

            if (interval == Interval.Lower) {
                double x = Mu - C * QuantilePade.Value(p);

                return x;
            }
            else {
                double x = Mu + C * QuantilePade.Value(p);

                return x;
            }
        }

        public double Sample(Random random) {
            double u = random.NextUniformOpenInterval01() - 0.5d;
            double w = random.NextUniformOpenInterval01();

            double cu = CosPi(u);

            double r = SinPi(u * 0.5d) * CosPi(u * 0.5) / (Log(w) * cu * cu);

            double v = r * C + Mu;

            return v;
        }

        public bool Symmetric => true;

        public double Median => Mu;

        public double Mode => Mu;

        public double Mean => NaN;

        public double Variance => NaN;

        public double Skewness => NaN;

        public double Kurtosis => NaN;

        public double Entropy => entropy_base + Log(C);

        public double Alpha => 0.5d;

        public double Beta => 0d;

        public static SaSPoint5Distribution operator +(SaSPoint5Distribution dist1, SaSPoint5Distribution dist2) {
            return new(dist1.Mu + dist2.Mu, ExMath.Square(Sqrt(dist1.C) + Sqrt(dist2.C)));
        }

        public static SaSPoint5Distribution operator -(SaSPoint5Distribution dist1, SaSPoint5Distribution dist2) {
            return new(dist1.Mu - dist2.Mu, ExMath.Square(Sqrt(dist1.C) + Sqrt(dist2.C)));
        }

        public static SaSPoint5Distribution operator +(SaSPoint5Distribution dist, double s) {
            return new(dist.Mu + s, dist.C);
        }

        public static SaSPoint5Distribution operator -(SaSPoint5Distribution dist, double s) {
            return new(dist.Mu - s, dist.C);
        }

        public static SaSPoint5Distribution operator *(SaSPoint5Distribution dist, double k) {
            return new(dist.Mu * k, dist.C * k);
        }

        public static SaSPoint5Distribution operator /(SaSPoint5Distribution dist, double k) {
            return new(dist.Mu / k, dist.C / k);
        }

        public override string ToString() {
            return $"{typeof(SaSPoint5Distribution).Name}[mu={Mu},c={C}]";
        }

        public string Formula => "p(x; mu, c) := stable_distribution(x; alpha = 1/2, beta = 0, mu, c)";

        private static class PDFPade {
            private static readonly (ReadOnlyCollection<double> numer, ReadOnlyCollection<double> denom) pade_plus_0_0p125 = new(
                new ReadOnlyCollection<double>([
                    6.36619772367581343076e-1,
                    2.17275699713513462507e2,
                    3.49063163361344578910e4,
                    3.40332906932698464252e6,
                    2.19485577044357440949e8,
                    9.66086435948730562464e9,
                    2.90571833690383003932e11,
                    5.83089315593106044683e12,
                    7.37911022713775715766e13,
                    5.26757196603002476852e14,
                    1.75780353683063527570e15,
                    1.85883041942144306222e15,
                    4.19828222275972713819e14,
                ]),
                new ReadOnlyCollection<double>([
                    1.00000000000000000000e0,
                    3.41295871011779138155e2,
                    5.48907134827349102297e4,
                    5.36641455324410261980e6,
                    3.48045461004960397915e8,
                    1.54920747349701741537e10,
                    4.76490595358644532404e11,
                    1.00104823128402735005e13,
                    1.39703522470411802507e14,
                    1.23724881334160220266e15,
                    6.47437580921138359461e15,
                    1.77627318260037604066e16,
                    2.04792815832538146160e16,
                    7.45102534638640681964e15,
                    3.68496090049571174527e14,
                ])
            );

            private static readonly (ReadOnlyCollection<double> numer, ReadOnlyCollection<double> denom) pade_plus_0p125_0p25 = new(
                new ReadOnlyCollection<double>([
                    4.35668401768623200524e-1,
                    7.12477357389655327116e0,
                    4.02466317948738993787e1,
                    9.04888497628205955839e1,
                    7.56175387288619211460e1,
                    1.26950253999694502457e1,
                    -6.59304802132933325219e-1,
                ]),
                new ReadOnlyCollection<double>([
                    1.00000000000000000000e0,
                    1.98623818041545101115e1,
                    1.52856383017632616759e2,
                    5.70706902111659740041e2,
                    1.06454927680197927878e3,
                    9.13160352749764887791e2,
                    2.58872466837209126618e2,
                ])
            );

            private static readonly (ReadOnlyCollection<double> numer, ReadOnlyCollection<double> denom) pade_plus_0p25_0p5 = new(
                new ReadOnlyCollection<double>([
                    2.95645445681747568732e-1,
                    2.23779537590791610124e0,
                    5.01302198171248036052e0,
                    2.76363131116340641935e0,
                    1.18134858311074670327e-1,
                    2.00287083462139382715e-2,
                    -7.53979800555375661516e-3,
                    1.37294648777729527395e-3,
                ]),
                new ReadOnlyCollection<double>([
                    1.00000000000000000000e0,
                    1.02879626214781666701e1,
                    3.85125274509784615691e1,
                    6.18474367367800231625e1,
                    3.77100050087302476029e1,
                    5.41866360740066443656e0,
                ])
            );

            private static readonly (ReadOnlyCollection<double> numer, ReadOnlyCollection<double> denom) pade_plus_0p5_1 = new(
                new ReadOnlyCollection<double>([
                    1.70762401725206223811e-1,
                    8.43343631021918972436e-1,
                    1.39703819152564365627e0,
                    8.75843324574692085009e-1,
                    1.86199552443747562584e-1,
                    7.35858280181579907616e-3,
                    -1.03693607694266081126e-4,
                ]),
                new ReadOnlyCollection<double>([
                    1.00000000000000000000e0,
                    6.73363440952557318819e0,
                    1.74288966619209299976e1,
                    2.15943268035083671893e1,
                    1.29818726981381859879e1,
                    3.40707211426946022041e0,
                    2.80229012541729457678e-1,
                ])
            );

            private static readonly (ReadOnlyCollection<double> numer, ReadOnlyCollection<double> denom) pade_plus_1_2 = new(
                new ReadOnlyCollection<double>([
                    8.61071469126041183247e-2,
                    1.69689585946245345838e-1,
                    1.09494833291892212033e-1,
                    2.76619622453130604637e-2,
                    2.44972748006913061509e-3,
                    4.09853605772288438003e-5,
                    -2.63561415158954865283e-7,
                ]),
                new ReadOnlyCollection<double>([
                    1.00000000000000000000e0,
                    3.04082856018856244947e0,
                    3.52558663323956252986e0,
                    1.94795523079701426332e0,
                    5.23956733400745421623e-1,
                    6.19453597593998871667e-2,
                    2.31061984192347753499e-3,
                ])
            );

            private static readonly (ReadOnlyCollection<double> numer, ReadOnlyCollection<double> denom) pade_plus_2_4 = new(
                new ReadOnlyCollection<double>([
                    3.91428580496513429479e-2,
                    4.07162484034780126757e-2,
                    1.43342733342753081931e-2,
                    2.01622178115394696215e-3,
                    1.00648013467757737201e-4,
                    9.51545046750892356441e-7,
                    -3.56598940936439037087e-9,
                ]),
                new ReadOnlyCollection<double>([
                    1.00000000000000000000e0,
                    1.63904431617187026619e0,
                    1.03812003196677309121e0,
                    3.18144310790210668797e-1,
                    4.81930155615666517263e-2,
                    3.25435391589941361778e-3,
                    7.01626957128181647457e-5,
                ])
            );

            private static readonly (ReadOnlyCollection<double> numer, ReadOnlyCollection<double> denom) pade_plus_4_8 = new(
                new ReadOnlyCollection<double>([
                    1.65057384221262866484e-2,
                    8.05429762031495873704e-3,
                    1.35249234647852784985e-3,
                    9.18685252682786794440e-5,
                    2.23447790937806602674e-6,
                    1.03176916111395079569e-8,
                    -1.94913182592441292094e-11,
                ]),
                new ReadOnlyCollection<double>([
                    1.00000000000000000000e0,
                    8.10113554189626079232e-1,
                    2.54175325409968367580e-1,
                    3.87119072807894983910e-2,
                    2.92520770162792443587e-3,
                    9.89094130526684467420e-5,
                    1.07148513311070719488e-6,
                ])
            );

            private static readonly (ReadOnlyCollection<double> numer, ReadOnlyCollection<double> denom) pade_plus_8_16 = new(
                new ReadOnlyCollection<double>([
                    6.60044810497290557553e-3,
                    1.59342644994950292031e-3,
                    1.32429706922966110874e-4,
                    4.45378136978435909660e-6,
                    5.36409958111394628239e-8,
                    1.22293787679910067873e-10,
                    -1.16300443044165216564e-13,
                ]),
                new ReadOnlyCollection<double>([
                    1.00000000000000000000e0,
                    4.10446485803039594111e-1,
                    6.51887342399859289520e-2,
                    5.02151225308643905366e-3,
                    1.91741179639551137839e-4,
                    3.27316600311598190022e-6,
                    1.78840301213102212857e-8,
                ])
            );

            private static readonly (ReadOnlyCollection<double> numer, ReadOnlyCollection<double> denom) pade_plus_16_32 = new(
                new ReadOnlyCollection<double>([
                    2.54339461777955741686e-3,
                    3.10069525357852579756e-4,
                    1.30082682796085732756e-5,
                    2.20715868479255585050e-7,
                    1.33996659756026452288e-9,
                    1.53505360463827994365e-12,
                    -7.42649416356965421308e-16,
                ]),
                new ReadOnlyCollection<double>([
                    1.00000000000000000000e0,
                    2.09203384450859785642e-1,
                    1.69422626897631306130e-2,
                    6.65649059670689720386e-4,
                    1.29654785666009849481e-5,
                    1.12886139474560969619e-7,
                    3.14420104899170413840e-10,
                ])
            );

            private static readonly (ReadOnlyCollection<double> numer, ReadOnlyCollection<double> denom) pade_plus_32_64 = new(
                new ReadOnlyCollection<double>([
                    9.55085695067883584460e-4,
                    5.86125496733202756668e-5,
                    1.23753971325810931282e-6,
                    1.05643819745933041408e-8,
                    3.22502949410095015524e-11,
                    1.85366144680157942079e-14,
                    -4.53975807317403152058e-18,
                ]),
                new ReadOnlyCollection<double>([
                    1.00000000000000000000e0,
                    1.05980850386474826374e-1,
                    4.34966042652000070674e-3,
                    8.66341538387446465700e-5,
                    8.55608082202236124363e-7,
                    3.77719968378509293354e-9,
                    5.33287361559571716670e-12,
                ])
            );

            private static readonly (ReadOnlyCollection<double> numer, ReadOnlyCollection<double> denom) pade_plus_limit = new(
                new ReadOnlyCollection<double>([
                    1.99471140200716338970e-1,
                    -1.93310094131437487158e-2,
                    -8.44282614309073196195e-3,
                    3.47296024282356038069e-3,
                    -4.05398011689821941383e-4,
                ]),
                new ReadOnlyCollection<double>([
                    1.00000000000000000000e0,
                    7.00973251258577238892e-1,
                    2.66969681258835723157e-1,
                    5.51785147503612200456e-2,
                    6.50130030979966274341e-3,
                ])
            );

            public static double Value(double x) {
                x = Abs(x);

                double y;
                if (x <= 0.125d) {
                    Debug.WriteLine("pade minimum segment passed");

                    y = ApproxUtil.Pade(x, pade_plus_0_0p125);
                }
                else if (x <= 0.25d) {
                    y = ApproxUtil.Pade(x - 0.125d, pade_plus_0p125_0p25);
                }
                else if (x <= 0.5d) {
                    y = ApproxUtil.Pade(x - 0.25d, pade_plus_0p25_0p5);
                }
                else if (x <= 1d) {
                    y = ApproxUtil.Pade(x - 0.5d, pade_plus_0p5_1);
                }
                else if (x <= 2d) {
                    y = ApproxUtil.Pade(x - 1d, pade_plus_1_2);
                }
                else if (x <= 4d) {
                    y = ApproxUtil.Pade(x - 2d, pade_plus_2_4);
                }
                else if (x <= 8d) {
                    y = ApproxUtil.Pade(x - 4d, pade_plus_4_8);
                }
                else if (x <= 16d) {
                    y = ApproxUtil.Pade(x - 8d, pade_plus_8_16);
                }
                else if (x <= 32d) {
                    y = ApproxUtil.Pade(x - 16d, pade_plus_16_32);
                }
                else if (x <= 64d) {
                    y = ApproxUtil.Pade(x - 32d, pade_plus_32_64);
                }
                else {
                    double v = Sqrt(x);
                    double u = 1d / v;

                    y = ApproxUtil.Pade(u, pade_plus_limit) * (u * u * u);
                }

                return y;
            }
        }

        private static class CDFPade {
            private static readonly (ReadOnlyCollection<double> numer, ReadOnlyCollection<double> denom) pade_plus_0_0p5 = new(
                new ReadOnlyCollection<double>([
                    5.00000000000000000000e-1,
                    1.11530082549581486148e2,
                    1.18564167533523512811e4,
                    7.51503793077701705413e5,
                    3.05648233678438482191e7,
                    8.12176734530090957088e8,
                    1.39533182836234507573e10,
                    1.50394359286077974212e11,
                    9.79057903542935575811e11,
                    3.73800992855150140014e12,
                    8.12697090329432868343e12,
                    9.63154058643818290870e12,
                    5.77714904017642642181e12,
                    1.53321958252091815685e12,
                    1.36220966258718212359e11,
                    1.70766655065405022702e9,
                ]),
                new ReadOnlyCollection<double>([
                    1.00000000000000000000e0,
                    2.24333404643898143947e2,
                    2.39984636687021023600e4,
                    1.53353791432086858132e6,
                    6.30764952479861776476e7,
                    1.70405769169309597488e9,
                    3.00381227010195289341e10,
                    3.37519046677507392667e11,
                    2.35001610518109063314e12,
                    9.90961948200767679416e12,
                    2.47066673978544828258e13,
                    3.51442593932882610556e13,
                    2.68891431106117733130e13,
                    9.99723484253582494535e12,
                    1.49190229409236772612e12,
                    5.68752980146893975323e10,
                ])
            );

            private static readonly (ReadOnlyCollection<double> numer, ReadOnlyCollection<double> denom) pade_plus_0p5_1 = new(
                new ReadOnlyCollection<double>([
                    3.31309550000758082456e-1,
                    1.63012162307622129396e0,
                    2.97763161467248770571e0,
                    2.49277948739575294031e0,
                    9.49619262302649586821e-1,
                    1.38360148984087584165e-1,
                    4.00812864075652334798e-3,
                    -4.82051978765960490940e-5,
                ]),
                new ReadOnlyCollection<double>([
                    1.00000000000000000000e0,
                    5.43565383128046471592e0,
                    1.13265160672130133152e1,
                    1.13352316246726435292e1,
                    5.56671465170409694873e0,
                    1.21011708389501479550e0,
                    8.34618282872428849500e-2,
                ])
            );

            private static readonly (ReadOnlyCollection<double> numer, ReadOnlyCollection<double> denom) pade_plus_1_2 = new(
                new ReadOnlyCollection<double>([
                    2.71280312689343248819e-1,
                    7.44610837974139249205e-1,
                    7.17844128359406982825e-1,
                    2.98789060945288850507e-1,
                    5.22747411439102272576e-2,
                    3.06447984437786430265e-3,
                    2.60407071021044908690e-5,
                ]),
                new ReadOnlyCollection<double>([
                    1.00000000000000000000e0,
                    3.06221257507188300824e0,
                    3.44827372231472308047e0,
                    1.78166113338930668519e0,
                    4.25580478492907232687e-1,
                    4.09983847731128510426e-2,
                    1.04343172183467651240e-3,
                ])
            );

            private static readonly (ReadOnlyCollection<double> numer, ReadOnlyCollection<double> denom) pade_plus_2_4 = new(
                new ReadOnlyCollection<double>([
                    2.13928162275383716645e-1,
                    2.35139109235828185307e-1,
                    9.35967515134932733243e-2,
                    1.64310489592753858417e-2,
                    1.23186728989215889119e-3,
                    3.13500969261032539402e-5,
                    1.17021346758965979212e-7,
                ]),
                new ReadOnlyCollection<double>([
                    1.00000000000000000000e0,
                    1.28212183177829510267e0,
                    6.17321009406850420793e-1,
                    1.38400318019319970893e-1,
                    1.44994794535896837497e-2,
                    6.17774446282546623636e-4,
                    7.00521050169239269819e-6,
                ])
            );

            private static readonly (ReadOnlyCollection<double> numer, ReadOnlyCollection<double> denom) pade_plus_4_8 = new(
                new ReadOnlyCollection<double>([
                    1.63772802979087193656e-1,
                    9.69009603942214234119e-2,
                    2.08261725719828138744e-2,
                    1.97965182693146960970e-3,
                    8.05499273532204276894e-5,
                    1.11401971145777879684e-6,
                    2.25932082770588727842e-9,
                ]),
                new ReadOnlyCollection<double>([
                    1.00000000000000000000e0,
                    6.92463563872865541733e-1,
                    1.80720987166755982366e-1,
                    2.20416647324531054557e-2,
                    1.26052070140663063778e-3,
                    2.93967534265875431639e-5,
                    1.82706995042259549615e-7,
                ])
            );

            private static readonly (ReadOnlyCollection<double> numer, ReadOnlyCollection<double> denom) pade_plus_8_16 = new(
                new ReadOnlyCollection<double>([
                    1.22610122564874280532e-1,
                    3.70273222121572231593e-2,
                    4.06083618461789591121e-3,
                    1.96898134215932126299e-4,
                    4.08421066512186972853e-6,
                    2.87707419853226244584e-8,
                    2.96850126180387702894e-11,
                ]),
                new ReadOnlyCollection<double>([
                    1.00000000000000000000e0,
                    3.55825191301363023576e-1,
                    4.77251766176046719729e-2,
                    2.99136605131226103925e-3,
                    8.78895785432321899939e-5,
                    1.05235770624006494709e-6,
                    3.35423877769913468556e-9,
                ])
            );

            private static readonly (ReadOnlyCollection<double> numer, ReadOnlyCollection<double> denom) pade_plus_16_32 = new(
                new ReadOnlyCollection<double>([
                    9.03056141356415077080e-2,
                    1.37568904417652631821e-2,
                    7.60947271383247418831e-4,
                    1.86048302967560067128e-5,
                    1.94537860496575427218e-7,
                    6.90524093915996283104e-10,
                    3.58808434477817122371e-13,
                ]),
                new ReadOnlyCollection<double>([
                    1.00000000000000000000e0,
                    1.80501347735272292079e-1,
                    1.22807958286146936376e-2,
                    3.90421541115275676253e-4,
                    5.81669449234915057779e-6,
                    3.53005415676201803667e-8,
                    5.69883025435873921433e-11,
                ])
            );

            private static readonly (ReadOnlyCollection<double> numer, ReadOnlyCollection<double> denom) pade_plus_32_64 = new(
                new ReadOnlyCollection<double>([
                    6.57333571766941474226e-2,
                    5.02795551798163084224e-3,
                    1.39633616037997111325e-4,
                    1.71386564634533872559e-6,
                    8.99508156357247137439e-9,
                    1.60229460572297160486e-11,
                    4.17711709622960498456e-15,
                ]),
                new ReadOnlyCollection<double>([
                    1.00000000000000000000e0,
                    9.10198637347368265508e-2,
                    3.12263472357578263712e-3,
                    5.00524795130325614005e-5,
                    3.75913188747149725195e-7,
                    1.14970132098893394023e-9,
                    9.34957119271300093120e-13,
                ])
            );

            private static readonly (ReadOnlyCollection<double> numer, ReadOnlyCollection<double> denom) pade_plus_limit = new(
                new ReadOnlyCollection<double>([
                    3.98942280401432677940e-1,
                    8.12222388783621449146e-2,
                    1.68515703707271703934e-2,
                    2.19801627205374824460e-3,
                    -5.63321705854968264807e-5,
                ]),
                new ReadOnlyCollection<double>([
                    1.00000000000000000000e0,
                    6.02536240902768558315e-1,
                    1.99284471400121092380e-1,
                    3.48012577961755452113e-2,
                    3.38545004473058881799e-3,
                ])
            );

            public static double Value(double x) {
                if (IsNegative(x)) {
                    return 1d - Value(-x);
                }

                double y;
                if (x <= 0.5d) {
                    Debug.WriteLine("pade minimum segment passed");

                    y = ApproxUtil.Pade(x, pade_plus_0_0p5);
                }
                else if (x <= 1d) {
                    y = ApproxUtil.Pade(x - 0.5d, pade_plus_0p5_1);
                }
                else if (x <= 2d) {
                    y = ApproxUtil.Pade(x - 1d, pade_plus_1_2);
                }
                else if (x <= 4d) {
                    y = ApproxUtil.Pade(x - 2d, pade_plus_2_4);
                }
                else if (x <= 8d) {
                    y = ApproxUtil.Pade(x - 4d, pade_plus_4_8);
                }
                else if (x <= 16d) {
                    y = ApproxUtil.Pade(x - 8d, pade_plus_8_16);
                }
                else if (x <= 32d) {
                    y = ApproxUtil.Pade(x - 16d, pade_plus_16_32);
                }
                else if (x <= 64d) {
                    y = ApproxUtil.Pade(x - 32d, pade_plus_32_64);
                }
                else {
                    double v = Sqrt(x);
                    double u = 1d / v;

                    y = ApproxUtil.Pade(u, pade_plus_limit) * u;
                }

                return y;
            }
        }

        private static class QuantilePade {
            private static readonly (ReadOnlyCollection<double> numer, ReadOnlyCollection<double> denom) pade_plus_expm1_1p125 = new(
                new ReadOnlyCollection<double>([
                    0.00000000000000000000e0,
                    1.36099130643975127045e-1,
                    2.19634434498311523885e1,
                    1.70276954848343179287e3,
                    8.02187341786354339306e4,
                    2.48750112198456813443e6,
                    5.20617858300443231437e7,
                    7.31202030685167303439e8,
                    6.66061403138355591915e9,
                    3.65687892725590813998e10,
                    1.06061776220305595494e11,
                    1.23930642673461465346e11,
                    1.49986408149520127078e10,
                    -6.17325587219357123900e8,
                ]),
                new ReadOnlyCollection<double>([
                    1.00000000000000000000e0,
                    1.63111146753825227716e2,
                    1.27864461509685444043e4,
                    6.10371533241799228037e5,
                    1.92422115963507708309e7,
                    4.11544185502250709497e8,
                    5.95343302992055062258e9,
                    5.65615858889758369947e10,
                    3.30833154992293143503e11,
                    1.06032392136054207216e12,
                    1.50071282012095447931e12,
                    5.43552396263989180433e11,
                    9.57434915768660935004e10,
                ])
            );

            private static readonly (ReadOnlyCollection<double> numer, ReadOnlyCollection<double> denom) pade_plus_expm1p125_1p25 = new(
                new ReadOnlyCollection<double>([
                    1.46698650748920243698e-2,
                    3.58380131788385557227e-1,
                    3.39153750029553194566e0,
                    1.55457424873957272207e1,
                    3.44403897039657057261e1,
                    3.01881531964962975320e1,
                    2.77679052294606319767e0,
                    -7.76665288232972435969e-2,
                ]),
                new ReadOnlyCollection<double>([
                    1.00000000000000000000e0,
                    1.72584280323876188464e1,
                    1.11983518800147654866e2,
                    3.25969893054048132145e2,
                    3.91978809680672051666e2,
                    1.29874252720714897530e2,
                    2.08740114519610102248e1,
                ])
            );

            private static readonly (ReadOnlyCollection<double> numer, ReadOnlyCollection<double> denom) pade_plus_expm1p25_1p5 = new(
                new ReadOnlyCollection<double>([
                    2.69627866689346445458e-2,
                    3.23091180507445216811e-1,
                    1.42164019533549860681e0,
                    2.74613170828120023406e0,
                    2.07865023346180997996e0,
                    2.53267176863740856907e-1,
                    -2.55816250186301841152e-2,
                    3.02683750470398342224e-3,
                ]),
                new ReadOnlyCollection<double>([
                    1.00000000000000000000e0,
                    8.55049920135376003042e0,
                    2.48726119139047911316e1,
                    2.79519589592198994574e1,
                    9.88212916161823866098e0,
                    1.39749417956251951564e0,
                ])
            );

            private static readonly (ReadOnlyCollection<double> numer, ReadOnlyCollection<double> denom) pade_plus_expm1p5_2 = new(
                new ReadOnlyCollection<double>([
                    4.79518653373241051274e-2,
                    3.81837125793765918564e-1,
                    1.13370353708146321188e0,
                    1.55218145762186846509e0,
                    9.60938271141036509605e-1,
                    2.11811755464425606950e-1,
                    8.84533960603915742831e-3,
                    1.73314614571009160225e-3,
                    -3.63491208733876986098e-5,
                ]),
                new ReadOnlyCollection<double>([
                    1.00000000000000000000e0,
                    6.36954463000253710936e0,
                    1.40601897306833147611e1,
                    1.33838075106916667084e1,
                    5.60958095533108032859e0,
                    1.11796035623375210182e0,
                    1.12508482637488861060e-1,
                    5.18503975949799718538e-3,
                ])
            );

            private static readonly (ReadOnlyCollection<double> numer, ReadOnlyCollection<double> denom) pade_plus_expm2_4 = new(
                new ReadOnlyCollection<double>([
                    8.02395484493329835881e-2,
                    2.46132933068351274622e-1,
                    2.81820176867119231101e-1,
                    1.47754061028371025893e-1,
                    3.54638964490281023406e-2,
                    3.99998730093393774294e-3,
                    3.81581928434827040262e-4,
                    1.82520920154354221101e-5,
                    -2.06151396745690348445e-7,
                    6.77986548138011345849e-9,
                ]),
                new ReadOnlyCollection<double>([
                    1.00000000000000000000e0,
                    2.39244329037830026691e0,
                    2.12683465416376620896e0,
                    9.02612272334554457823e-1,
                    2.06667959191488815314e-1,
                    2.79328968525257867541e-2,
                    2.28216286216537879937e-3,
                    1.04195690531437767679e-4,
                ])
            );

            private static readonly (ReadOnlyCollection<double> numer, ReadOnlyCollection<double> denom) pade_plus_expm4_8 = new(
                new ReadOnlyCollection<double>([
                    1.39293493266195561875e-1,
                    1.26741380938661691592e-1,
                    4.31117040307200265931e-2,
                    7.50528269269498076949e-3,
                    8.63100497178570310436e-4,
                    6.75686286034521991703e-5,
                    3.11102625473120771882e-6,
                    9.63513655399980075083e-8,
                    -6.40223609013005302318e-11,
                ]),
                new ReadOnlyCollection<double>([
                    1.00000000000000000000e0,
                    8.11234548272888947555e-1,
                    2.63525516991753831892e-1,
                    4.77118226533147280522e-2,
                    5.46090741266888954909e-3,
                    4.15325425646862026425e-4,
                    2.02377681998442384863e-5,
                    5.79823311154876056655e-7,
                ])
            );

            private static readonly (ReadOnlyCollection<double> numer, ReadOnlyCollection<double> denom) pade_plus_expm8_16 = new(
                new ReadOnlyCollection<double>([
                    1.57911660613037760235e-1,
                    5.59740955695099219682e-2,
                    8.92895854008560399142e-3,
                    8.88795299273855801726e-4,
                    5.66358335596607738071e-5,
                    2.46733195253941569922e-6,
                    6.44829870181825872501e-8,
                    7.62193242864380357931e-10,
                    -7.82035413331699873450e-14,
                ]),
                new ReadOnlyCollection<double>([
                    1.00000000000000000000e0,
                    3.49007782566002620811e-1,
                    5.65303702876260444572e-2,
                    5.54316442661801299351e-3,
                    3.58498995501703237922e-4,
                    1.53872913968336341278e-5,
                    4.08512152326482573624e-7,
                    4.72959615756470826429e-9,
                ])
            );

            private static readonly (ReadOnlyCollection<double> numer, ReadOnlyCollection<double> denom) pade_plus_expm16_32 = new(
                new ReadOnlyCollection<double>([
                    1.59150086070234563099e-1,
                    6.07144002506911115092e-2,
                    1.10026443723891740392e-2,
                    1.24892739209332398698e-3,
                    9.82922518655171276487e-5,
                    5.58366837526347222893e-6,
                    2.29005408647580194007e-7,
                    6.44325718317518336404e-9,
                    1.05110361316230054467e-10,
                    1.48083450629432857655e-18,
                ]),
                new ReadOnlyCollection<double>([
                    1.00000000000000000000e0,
                    3.81470315977341203351e-1,
                    6.91330250512167919573e-2,
                    7.84712209182587717077e-3,
                    6.17595479676821181012e-4,
                    3.50829361179041199953e-5,
                    1.43889153071571504712e-6,
                    4.04840254888235877998e-8,
                    6.60429636407045050112e-10,
                ])
            );

            private static readonly (ReadOnlyCollection<double> numer, ReadOnlyCollection<double> denom) pade_plus_expm32_64 = new(
                new ReadOnlyCollection<double>([
                    1.59154943017783026201e-1,
                    6.91506515614472069475e-2,
                    1.44590186111155933843e-2,
                    1.92616138327724025421e-3,
                    1.79640147906775699469e-4,
                    1.30852535070639833809e-5,
                    5.55259657884038297268e-7,
                    3.50107118687544980820e-8,
                    -1.47102592933729597720e-22,
                ]),
                new ReadOnlyCollection<double>([
                    1.00000000000000000000e0,
                    4.34486357752330500669e-1,
                    9.08486933075320995164e-2,
                    1.21024289017243304241e-2,
                    1.12871233794777525784e-3,
                    8.22170725751776749123e-5,
                    3.48879932410650101194e-6,
                    2.19978790407451988423e-7,
                ])
            );

            public static double Value(double x) {
                if (x > 0.5) {
                    return -Value(1d - x);
                }

                double v;
                int exponent = ILogB(x);

                if (exponent >= -2) {
                    double u = -Log2(ScaleB(x, 1));

                    if (u <= 0.125) {
                        v = ApproxUtil.Pade(u, pade_plus_expm1_1p125);
                    }
                    else if (u <= 0.25) {
                        v = ApproxUtil.Pade(u - 0.125, pade_plus_expm1p125_1p25);
                    }
                    else if (u <= 0.5) {
                        v = ApproxUtil.Pade(u - 0.25, pade_plus_expm1p25_1p5);
                    }
                    else {
                        v = ApproxUtil.Pade(u - 0.5, pade_plus_expm1p5_2);
                    }
                }
                else if (exponent >= -4) {
                    v = ApproxUtil.Pade(-Log2(ScaleB(x, 2)), pade_plus_expm2_4);
                }
                else if (exponent >= -8) {
                    v = ApproxUtil.Pade(-Log2(ScaleB(x, 4)), pade_plus_expm4_8);
                }
                else if (exponent >= -16) {
                    v = ApproxUtil.Pade(-Log2(ScaleB(x, 8)), pade_plus_expm8_16);
                }
                else if (exponent >= -32) {
                    v = ApproxUtil.Pade(-Log2(ScaleB(x, 16)), pade_plus_expm16_32);
                }
                else if (exponent >= -64) {
                    v = ApproxUtil.Pade(-Log2(ScaleB(x, 32)), pade_plus_expm32_64);
                }
                else {
                    v = ScaleB(1d / Pi, -1);
                }

                double y = v / (x * x);

                return y;
            }
        }
    }
}