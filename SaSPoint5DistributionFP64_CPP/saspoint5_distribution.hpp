// Author and Approximation Formula Coefficient Generator: T.Yoshimura
// Github: https://github.com/tk-yoshimura
// Original Code: https://github.com/tk-yoshimura/SaSPoint5DistributionFP64
// C++20 implement

#include <vector>
#include <cmath>
#include <cassert>
#include <numbers>
#include <limits>

using namespace std;
using namespace std::numbers;

double poly(double x, vector<double> coef) {
    double s = coef[coef.size() - 1];

    for (int i = coef.size() - 2; i >= 0; i--) {
        s = s * x + coef[i];
    }

    return s;
}

double pade(double x, vector<double> numer, vector<double> denom) {
    double sc = poly(x, numer), sd = poly(x, denom);

    assert(sd >= 0.5);

    return sc / sd;
}

double saspoint5_pdf(double x) {
    static const vector<double> pade_plus_0_0p25_numer = {
        6.36619772367581343076e-1,
        2.08695290618996857630e2,
        3.27138231260799880672e4,
        3.12635260622309707764e6,
        1.98048081337216511381e8,
        8.56483679198669196462e9,
        2.52897673008030765179e11,
        4.97565108033799014273e12,
        6.17809824993213970941e13,
        4.39359884488046050026e14,
        1.62407973388917606167e15,
        3.05602866995962977259e15,
        4.16283837663175466754e15,
        4.33932052046202879986e15,
        1.98037038829086082130e15,
        2.11673879981348951641e14,
    };
    static const vector<double> pade_plus_0_0p25_denom = {
        1.00000000000000000000e0,
        3.27817795923715088045e2,
        5.14467532018582990645e4,
        4.93053225788651363719e6,
        3.14164883876153657006e8,
        1.37444895081586710404e10,
        4.15331400451571813658e11,
        8.56868997354451418025e12,
        1.17651029229388237648e14,
        1.03622647547110519740e15,
        5.60383117050645616112e15,
        1.78609307498162312531e16,
        3.40114300916463634855e16,
        4.63563483570838873155e16,
        4.77894054692871071343e16,
        2.65130046251435019059e16,
        5.00269030614708036673e15,
        1.28890058960331124978e14,
    };
    static const vector<double> pade_plus_0p25_0p5_numer = {
        2.95645445681747568732e-1,
        1.97411200998490091540e0,
        2.96639538497923278703e0,
        -2.11314593078651532185e0,
        -3.24647624911685089557e0,
        -4.16259541605596944447e-1,
        1.40473115910450455571e-2,
    };
    static const vector<double> pade_plus_0p25_0p5_denom = {
        1.00000000000000000000e0,
        9.39607211623133600737e0,
        2.91651053913427485821e1,
        2.56594161726268984972e1,
        -2.43932432079171316468e1,
        -3.87779015677750197072e1,
        -9.84153122176909305802e0,
    };
    static const vector<double> pade_plus_0p5_1_numer = {
        1.70762401725206223811e-1,
        8.43343631021918972436e-1,
        1.39703819152564365627e0,
        8.75843324574692085009e-1,
        1.86199552443747562584e-1,
        7.35858280181579907616e-3,
        -1.03693607694266081126e-4,
    };
    static const vector<double> pade_plus_0p5_1_denom = {
        1.00000000000000000000e0,
        6.73363440952557318819e0,
        1.74288966619209299976e1,
        2.15943268035083671893e1,
        1.29818726981381859879e1,
        3.40707211426946022041e0,
        2.80229012541729457678e-1,
    };
    static const vector<double> pade_plus_1_2_numer = {
        8.61071469126041183247e-2,
        1.69689585946245345838e-1,
        1.09494833291892212033e-1,
        2.76619622453130604637e-2,
        2.44972748006913061509e-3,
        4.09853605772288438003e-5,
        -2.63561415158954865283e-7,
    };
    static const vector<double> pade_plus_1_2_denom = {
        1.00000000000000000000e0,
        3.04082856018856244947e0,
        3.52558663323956252986e0,
        1.94795523079701426332e0,
        5.23956733400745421623e-1,
        6.19453597593998871667e-2,
        2.31061984192347753499e-3,
    };
    static const vector<double> pade_plus_2_4_numer = {
        3.91428580496513429479e-2,
        4.07162484034780126757e-2,
        1.43342733342753081931e-2,
        2.01622178115394696215e-3,
        1.00648013467757737201e-4,
        9.51545046750892356441e-7,
        -3.56598940936439037087e-9,
    };
    static const vector<double> pade_plus_2_4_denom = {
        1.00000000000000000000e0,
        1.63904431617187026619e0,
        1.03812003196677309121e0,
        3.18144310790210668797e-1,
        4.81930155615666517263e-2,
        3.25435391589941361778e-3,
        7.01626957128181647457e-5,
    };
    static const vector<double> pade_plus_4_8_numer = {
        1.65057384221262866484e-2,
        8.05429762031495873704e-3,
        1.35249234647852784985e-3,
        9.18685252682786794440e-5,
        2.23447790937806602674e-6,
        1.03176916111395079569e-8,
        -1.94913182592441292094e-11,
    };
    static const vector<double> pade_plus_4_8_denom = {
        1.00000000000000000000e0,
        8.10113554189626079232e-1,
        2.54175325409968367580e-1,
        3.87119072807894983910e-2,
        2.92520770162792443587e-3,
        9.89094130526684467420e-5,
        1.07148513311070719488e-6,
    };
    static const vector<double> pade_plus_8_16_numer = {
        6.60044810497290557553e-3,
        1.59342644994950292031e-3,
        1.32429706922966110874e-4,
        4.45378136978435909660e-6,
        5.36409958111394628239e-8,
        1.22293787679910067873e-10,
        -1.16300443044165216564e-13,
    };
    static const vector<double> pade_plus_8_16_denom = {
        1.00000000000000000000e0,
        4.10446485803039594111e-1,
        6.51887342399859289520e-2,
        5.02151225308643905366e-3,
        1.91741179639551137839e-4,
        3.27316600311598190022e-6,
        1.78840301213102212857e-8,
    };
    static const vector<double> pade_plus_16_32_numer = {
        2.54339461777955741686e-3,
        3.10069525357852579756e-4,
        1.30082682796085732756e-5,
        2.20715868479255585050e-7,
        1.33996659756026452288e-9,
        1.53505360463827994365e-12,
        -7.42649416356965421308e-16,
    };
    static const vector<double> pade_plus_16_32_denom = {
        1.00000000000000000000e0,
        2.09203384450859785642e-1,
        1.69422626897631306130e-2,
        6.65649059670689720386e-4,
        1.29654785666009849481e-5,
        1.12886139474560969619e-7,
        3.14420104899170413840e-10,
    };
    static const vector<double> pade_plus_32_64_numer = {
        9.55085695067883584460e-4,
        5.86125496733202756668e-5,
        1.23753971325810931282e-6,
        1.05643819745933041408e-8,
        3.22502949410095015524e-11,
        1.85366144680157942079e-14,
        -4.53975807317403152058e-18,
    };
    static const vector<double> pade_plus_32_64_denom = {
        1.00000000000000000000e0,
        1.05980850386474826374e-1,
        4.34966042652000070674e-3,
        8.66341538387446465700e-5,
        8.55608082202236124363e-7,
        3.77719968378509293354e-9,
        5.33287361559571716670e-12,
    };
    static const vector<double> pade_plus_limit_numer = {
        1.99471140200716338970e-1,
        -1.93310094131437487158e-2,
        -8.44282614309073196195e-3,
        3.47296024282356038069e-3,
        -4.05398011689821941383e-4,
    };
    static const vector<double> pade_plus_limit_denom = {
        1.00000000000000000000e0,
        7.00973251258577238892e-1,
        2.66969681258835723157e-1,
        5.51785147503612200456e-2,
        6.50130030979966274341e-3,
    };

    x = abs(x);

    double y;
    if (x <= 0.25) {
        y = pade(x, pade_plus_0_0p25_numer, pade_plus_0_0p25_denom);
    }
    else if (x <= 0.5) {
        y = pade(x - 0.25, pade_plus_0p25_0p5_numer, pade_plus_0p25_0p5_denom);
    }
    else if (x <= 1) {
        y = pade(x - 0.5, pade_plus_0p5_1_numer, pade_plus_0p5_1_denom);
    }
    else if (x <= 2) {
        y = pade(x - 1, pade_plus_1_2_numer, pade_plus_1_2_denom);
    }
    else if (x <= 4) {
        y = pade(x - 2, pade_plus_2_4_numer, pade_plus_2_4_denom);
    }
    else if (x <= 8) {
        y = pade(x - 4, pade_plus_4_8_numer, pade_plus_4_8_denom);
    }
    else if (x <= 16) {
        y = pade(x - 8, pade_plus_8_16_numer, pade_plus_8_16_denom);
    }
    else if (x <= 32) {
        y = pade(x - 16, pade_plus_16_32_numer, pade_plus_16_32_denom);
    }
    else if (x <= 64) {
        y = pade(x - 32, pade_plus_32_64_numer, pade_plus_32_64_denom);
    }
    else {
        double v = sqrt(x);
        double u = 1 / v;

        y = pade(u, pade_plus_limit_numer, pade_plus_limit_denom) * (u * u * u);
    }

    return y;
}

double saspoint5_cdf(double x, bool complementary = false) {
    static const vector<double> pade_plus_0_0p5_numer = {
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
    };
    static const vector<double> pade_plus_0_0p5_denom = {
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
    };
    static const vector<double> pade_plus_0p5_1_numer = {
        3.31309550000758082456e-1,
        1.63012162307622129396e0,
        2.97763161467248770571e0,
        2.49277948739575294031e0,
        9.49619262302649586821e-1,
        1.38360148984087584165e-1,
        4.00812864075652334798e-3,
        -4.82051978765960490940e-5,
    };
    static const vector<double> pade_plus_0p5_1_denom = {
        1.00000000000000000000e0,
        5.43565383128046471592e0,
        1.13265160672130133152e1,
        1.13352316246726435292e1,
        5.56671465170409694873e0,
        1.21011708389501479550e0,
        8.34618282872428849500e-2,
    };
    static const vector<double> pade_plus_1_2_numer = {
        2.71280312689343248819e-1,
        7.44610837974139249205e-1,
        7.17844128359406982825e-1,
        2.98789060945288850507e-1,
        5.22747411439102272576e-2,
        3.06447984437786430265e-3,
        2.60407071021044908690e-5,
    };
    static const vector<double> pade_plus_1_2_denom = {
        1.00000000000000000000e0,
        3.06221257507188300824e0,
        3.44827372231472308047e0,
        1.78166113338930668519e0,
        4.25580478492907232687e-1,
        4.09983847731128510426e-2,
        1.04343172183467651240e-3,
    };
    static const vector<double> pade_plus_2_4_numer = {
        2.13928162275383716645e-1,
        2.35139109235828185307e-1,
        9.35967515134932733243e-2,
        1.64310489592753858417e-2,
        1.23186728989215889119e-3,
        3.13500969261032539402e-5,
        1.17021346758965979212e-7,
    };
    static const vector<double> pade_plus_2_4_denom = {
        1.00000000000000000000e0,
        1.28212183177829510267e0,
        6.17321009406850420793e-1,
        1.38400318019319970893e-1,
        1.44994794535896837497e-2,
        6.17774446282546623636e-4,
        7.00521050169239269819e-6,
    };
    static const vector<double> pade_plus_4_8_numer = {
        1.63772802979087193656e-1,
        9.69009603942214234119e-2,
        2.08261725719828138744e-2,
        1.97965182693146960970e-3,
        8.05499273532204276894e-5,
        1.11401971145777879684e-6,
        2.25932082770588727842e-9,
    };
    static const vector<double> pade_plus_4_8_denom = {
        1.00000000000000000000e0,
        6.92463563872865541733e-1,
        1.80720987166755982366e-1,
        2.20416647324531054557e-2,
        1.26052070140663063778e-3,
        2.93967534265875431639e-5,
        1.82706995042259549615e-7,
    };
    static const vector<double> pade_plus_8_16_numer = {
        1.22610122564874280532e-1,
        3.70273222121572231593e-2,
        4.06083618461789591121e-3,
        1.96898134215932126299e-4,
        4.08421066512186972853e-6,
        2.87707419853226244584e-8,
        2.96850126180387702894e-11,
    };
    static const vector<double> pade_plus_8_16_denom = {
        1.00000000000000000000e0,
        3.55825191301363023576e-1,
        4.77251766176046719729e-2,
        2.99136605131226103925e-3,
        8.78895785432321899939e-5,
        1.05235770624006494709e-6,
        3.35423877769913468556e-9,
    };
    static const vector<double> pade_plus_16_32_numer = {
        9.03056141356415077080e-2,
        1.37568904417652631821e-2,
        7.60947271383247418831e-4,
        1.86048302967560067128e-5,
        1.94537860496575427218e-7,
        6.90524093915996283104e-10,
        3.58808434477817122371e-13,
    };
    static const vector<double> pade_plus_16_32_denom = {
        1.00000000000000000000e0,
        1.80501347735272292079e-1,
        1.22807958286146936376e-2,
        3.90421541115275676253e-4,
        5.81669449234915057779e-6,
        3.53005415676201803667e-8,
        5.69883025435873921433e-11,
    };
    static const vector<double> pade_plus_32_64_numer = {
        6.57333571766941474226e-2,
        5.02795551798163084224e-3,
        1.39633616037997111325e-4,
        1.71386564634533872559e-6,
        8.99508156357247137439e-9,
        1.60229460572297160486e-11,
        4.17711709622960498456e-15,
    };
    static const vector<double> pade_plus_32_64_denom = {
        1.00000000000000000000e0,
        9.10198637347368265508e-2,
        3.12263472357578263712e-3,
        5.00524795130325614005e-5,
        3.75913188747149725195e-7,
        1.14970132098893394023e-9,
        9.34957119271300093120e-13,
    };
    static const vector<double> pade_plus_limit_numer = {
        3.98942280401432677940e-1,
        8.12222388783621449146e-2,
        1.68515703707271703934e-2,
        2.19801627205374824460e-3,
        -5.63321705854968264807e-5,
    };
    static const vector<double> pade_plus_limit_denom = {
        1.00000000000000000000e0,
        6.02536240902768558315e-1,
        1.99284471400121092380e-1,
        3.48012577961755452113e-2,
        3.38545004473058881799e-3,
    };

    bool inversion = (x <= 0) ^ complementary;

    x = abs(x);

    double y;
    if (x <= 0.5) {
        y = pade(x, pade_plus_0_0p5_numer, pade_plus_0_0p5_denom);
    }
    else if (x <= 1) {
        y = pade(x - 0.5, pade_plus_0p5_1_numer, pade_plus_0p5_1_denom);
    }
    else if (x <= 2) {
        y = pade(x - 1, pade_plus_1_2_numer, pade_plus_1_2_denom);
    }
    else if (x <= 4) {
        y = pade(x - 2, pade_plus_2_4_numer, pade_plus_2_4_denom);
    }
    else if (x <= 8) {
        y = pade(x - 4, pade_plus_4_8_numer, pade_plus_4_8_denom);
    }
    else if (x <= 16) {
        y = pade(x - 8, pade_plus_8_16_numer, pade_plus_8_16_denom);
    }
    else if (x <= 32) {
        y = pade(x - 16, pade_plus_16_32_numer, pade_plus_16_32_denom);
    }
    else if (x <= 64) {
        y = pade(x - 32, pade_plus_32_64_numer, pade_plus_32_64_denom);
    }
    else {
        double v = sqrt(x);
        double u = 1 / v;

        y = pade(u, pade_plus_limit_numer, pade_plus_limit_denom) * u;
    }

    y = inversion ? y : 1 - y;

    return y;
}

double saspoint5_quantile(double x, bool complementary = false) {
    static const vector<double> pade_plus_0p125_0p15625_numer = {
    4.49147047287704191198e0,
    3.71414591873162733815e1,
    1.03480125109733625638e2,
    1.00005220470908823354e2,
    2.74730973141070583590e0,
    -1.49189870645073451633e1,
    };
    static const vector<double> pade_plus_0p125_0p15625_denom = {
        1.00000000000000000000e0,
        -7.40327929970295513340e0,
        -2.07000208788561899835e1,
        1.26466884410752108685e2,
        2.98398241252531847810e2,
        2.44953613390623685184e1,
        -3.01956129294066642747e1,
    };
    static const vector<double> pade_plus_0p15625_0p1875_numer = {
        2.83944503273842706084e0,
        1.67404492991634154469e1,
        2.61760101641248660785e1,
        4.92318752899232235761e0,
        -1.33297300307054215288e0,
    };
    static const vector<double> pade_plus_0p15625_0p1875_denom = {
        1.00000000000000000000e0,
        -7.91824101512135661001e0,
        -4.36560728869362624099e-1,
        7.35318086832885976497e1,
        1.23880777723503405014e1,
    };
    static const vector<double> pade_plus_0p1875_0p25_numer = {
        1.28383277518932773555e0,
        1.54688305214846834739e1,
        6.21226143476123511484e1,
        8.79420018191270817296e1,
        1.79710864893179400702e1,
        -4.36296483289971333525e0,
    };
    static const vector<double> pade_plus_0p1875_0p25_denom = {
        1.00000000000000000000e0,
        1.53476487482968293892e-1,
        -3.27465714746398966728e1,
        1.45811276737475915785e0,
        2.40881382059354704440e2,
        4.56347287861047693971e1,
    };
    static const vector<double> pade_plus_0p25_0p3125_numer = {
        6.21997830440617022638e-1,
        8.11360490622525648490e0,
        3.54777843088616889615e1,
        5.57432092610261643096e1,
        1.52522387203702199687e1,
        -2.61310790608187412340e0,
    };
    static const vector<double> pade_plus_0p25_0p3125_denom = {
        1.00000000000000000000e0,
        1.51582630412651833313e0,
        -2.40350415645727043604e1,
        -1.99617786212255367547e1,
        1.43144810168931016975e2,
        3.94332982945258486603e1,
    };
    static const vector<double> pade_plus_0p3125_0p375_numer = {
        2.93390030300024489254e-1,
        5.16683967106977943664e0,
        2.97430566765511734325e1,
        6.04452611872185179046e1,
        2.33091461217488287282e1,
        -2.53187049065204826355e0,
    };
    static const vector<double> pade_plus_0p3125_0p375_denom = {
        1.00000000000000000000e0,
        4.70163738422098119446e0,
        -2.04200239512608550672e1,
        -5.86169085045103369448e1,
        1.38804873716322353304e2,
        6.12970262611389719488e1,
    };
    static const vector<double> pade_plus_0p375_0p4375_numer = {
        1.12403044179836065451e-1,
        9.66281090097717968770e0,
        2.82936597616711818703e2,
        3.70468489885038565990e3,
        2.20859955075149305099e4,
        4.88019507624858356395e4,
        4.74293142423451556386e3,
        -4.39007128051960125520e3,
    };
    static const vector<double> pade_plus_0p375_0p4375_denom = {
        1.00000000000000000000e0,
        6.64372788710680763873e1,
        1.14299363729551464845e3,
        5.21928413581303380919e3,
        -1.48642173255156825356e4,
        -7.79821178543143197493e4,
        1.58132628754840899695e5,
    };
    static const vector<double> pade_plus_0p4375_0p46875_numer = {
        5.11633396058323290276e-2,
        4.63682274468470911206e0,
        1.60248860748557199402e2,
        2.63720558501783652599e3,
        2.04317281067803852047e4,
        5.93233409946824883602e4,
        8.81331971324281446276e3,
    };
    static const vector<double> pade_plus_0p4375_0p46875_denom = {
        1.00000000000000000000e0,
        5.62583376348831271328e1,
        1.09471802023964334935e3,
        7.32578868098044330527e3,
        -5.79750588904401014870e3,
        -1.22300195620521909334e5,
        1.89789689002325527907e5,
    };
    static const vector<double> pade_plus_0p46875_0p484375_numer = {
        2.48262870167331472810e-2,
        3.15160321160735006235e0,
        1.35873615593591715436e2,
        2.57741528351772386498e3,
        1.84962907958895484391e4,
        2.77498760256500183962e4,
        4.66650326473025064328e4,
    };
    static const vector<double> pade_plus_0p46875_0p484375_denom = {
        1.00000000000000000000e0,
        6.15450908942569260141e1,
        1.31938240470054998868e3,
        7.41488495176428683684e3,
        -2.89904373344249493042e4,
    };
    static const vector<double> pade_plus_0p484375_0p5_numer = {
        0.00000000000000000000e0,
        1.57079632679489658420e0,
        3.68677664732235740582e2,
        5.12991233169657765833e4,
        4.15925967174184357395e6,
        2.20306977390050161287e8,
        7.02631670287223054378e9,
        1.27090189596731191250e11,
        8.38808550969735622972e11,
        -4.95994728465451538598e11,
    };
    static const vector<double> pade_plus_0p484375_0p5_denom = {
        1.00000000000000000000e0,
        2.34707490998847188599e2,
        3.26086881866984488126e4,
        2.63628459501019657057e6,
        1.38653708172292751211e8,
        4.34560305046086344625e9,
        7.44152155659186741721e10,
        3.45841494519572088548e11,
        -2.83594324749756534784e12,
        4.11908049498574834947e12,
    };
    static const vector<double> pade_plus_expm3_4_numer = {
        1.19511170764405535039e-1,
        1.37927972427045248330e-1,
        5.24650188881053184583e-2,
        8.03959641438786924480e-3,
        7.65196262623027097818e-4,
        4.48726026525832428753e-5,
        -2.39900494956368998141e-7,
    };
    static const vector<double> pade_plus_expm3_4_denom = {
        1.00000000000000000000e0,
        9.25160082350049994731e-1,
        3.05761874775646930171e-1,
        5.22137081244781154561e-2,
        4.97356364271483616303e-3,
        2.50262885907540348010e-4,
    };
    static const vector<double> pade_plus_expm4_8_numer = {
        1.39293493266195561875e-1,
        1.26741380938661691592e-1,
        4.31117040307200265931e-2,
        7.50528269269498076949e-3,
        8.63100497178570310436e-4,
        6.75686286034521991703e-5,
        3.11102625473120771882e-6,
        9.63513655399980075083e-8,
        -6.40223609013005302318e-11,
    };
    static const vector<double> pade_plus_expm4_8_denom = {
        1.00000000000000000000e0,
        8.11234548272888947555e-1,
        2.63525516991753831892e-1,
        4.77118226533147280522e-2,
        5.46090741266888954909e-3,
        4.15325425646862026425e-4,
        2.02377681998442384863e-5,
        5.79823311154876056655e-7,
    };
    static const vector<double> pade_plus_expm8_16_numer = {
        1.57911660613037760235e-1,
        5.59740955695099219682e-2,
        8.92895854008560399142e-3,
        8.88795299273855801726e-4,
        5.66358335596607738071e-5,
        2.46733195253941569922e-6,
        6.44829870181825872501e-8,
        7.62193242864380357931e-10,
        -7.82035413331699873450e-14,
    };
    static const vector<double> pade_plus_expm8_16_denom = {
        1.00000000000000000000e0,
        3.49007782566002620811e-1,
        5.65303702876260444572e-2,
        5.54316442661801299351e-3,
        3.58498995501703237922e-4,
        1.53872913968336341278e-5,
        4.08512152326482573624e-7,
        4.72959615756470826429e-9,
    };
    static const vector<double> pade_plus_expm16_32_numer = {
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
    };
    static const vector<double> pade_plus_expm16_32_denom = {
        1.00000000000000000000e0,
        3.81470315977341203351e-1,
        6.91330250512167919573e-2,
        7.84712209182587717077e-3,
        6.17595479676821181012e-4,
        3.50829361179041199953e-5,
        1.43889153071571504712e-6,
        4.04840254888235877998e-8,
        6.60429636407045050112e-10,
    };
    static const vector<double> pade_plus_expm32_64_numer = {
        1.59154943017783026201e-1,
        6.91506515614472069475e-2,
        1.44590186111155933843e-2,
        1.92616138327724025421e-3,
        1.79640147906775699469e-4,
        1.30852535070639833809e-5,
        5.55259657884038297268e-7,
        3.50107118687544980820e-8,
        -1.47102592933729597720e-22,
    };
    static const vector<double> pade_plus_expm32_64_denom = {
        1.00000000000000000000e0,
        4.34486357752330500669e-1,
        9.08486933075320995164e-2,
        1.21024289017243304241e-2,
        1.12871233794777525784e-3,
        8.22170725751776749123e-5,
        3.48879932410650101194e-6,
        2.19978790407451988423e-7,
    };
    static const vector<double> pade_plus_expm64_96_numer = {
        1.59154943091895335754e-1,
        1.01884113307066544176e-1,
    };
    static const vector<double> pade_plus_expm64_96_denom = {
        1.00000000000000000000e0,
        6.40156763765980698583e-1,
    };

    if (x > 0.5) {
        return -saspoint5_quantile(1 - x, complementary);
    }

    double y;
    if (x >= 0.125) {
        if (x <= 0.15625) {
            y = pade(0.15625 - x, pade_plus_0p125_0p15625_numer, pade_plus_0p125_0p15625_denom);
        }
        else if (x <= 0.1875) {
            y = pade(0.1875 - x, pade_plus_0p15625_0p1875_numer, pade_plus_0p15625_0p1875_denom);
        }
        else if (x <= 0.25) {
            y = pade(0.25 - x, pade_plus_0p1875_0p25_numer, pade_plus_0p1875_0p25_denom);
        }
        else if (x <= 0.3125) {
            y = pade(0.3125 - x, pade_plus_0p25_0p3125_numer, pade_plus_0p25_0p3125_denom);
        }
        else if (x <= 0.375) {
            y = pade(0.375 - x, pade_plus_0p3125_0p375_numer, pade_plus_0p3125_0p375_denom);
        }
        else if (x <= 0.4375) {
            y = pade(0.4375 - x, pade_plus_0p375_0p4375_numer, pade_plus_0p375_0p4375_denom);
        }
        else if (x <= 0.46875) {
            y = pade(0.46875 - x, pade_plus_0p4375_0p46875_numer, pade_plus_0p4375_0p46875_denom);
        }
        else if (x <= 0.484375) {
            y = pade(0.484375 - x, pade_plus_0p46875_0p484375_numer, pade_plus_0p46875_0p484375_denom);
        }
        else {
            y = pade(0.5 - x, pade_plus_0p484375_0p5_numer, pade_plus_0p484375_0p5_denom);
        }
    }
    else {
        double v;
        int exponent = ilogb(x);

        if (exponent >= -4) {
            v = pade(-log2(ldexp(x, 3)), pade_plus_expm3_4_numer, pade_plus_expm3_4_denom);
        }
        else if (exponent >= -8) {
            v = pade(-log2(ldexp(x, 4)), pade_plus_expm4_8_numer, pade_plus_expm4_8_denom);
        }
        else if (exponent >= -16) {
            v = pade(-log2(ldexp(x, 8)), pade_plus_expm8_16_numer, pade_plus_expm8_16_denom);
        }
        else if (exponent >= -32) {
            v = pade(-log2(ldexp(x, 16)), pade_plus_expm16_32_numer, pade_plus_expm16_32_denom);
        }
        else if (exponent >= -64) {
            v = pade(-log2(ldexp(x, 32)), pade_plus_expm32_64_numer, pade_plus_expm32_64_denom);
        }
        else if (exponent >= -96) {
            v = pade(-log2(ldexp(x, 64)), pade_plus_expm64_96_numer, pade_plus_expm64_96_denom);
        }
        else {
            v = ldexp(1 / pi, -1);
        }

        y = v / (x * x);
    }

    y = complementary ? y : -y;

    return y;
}