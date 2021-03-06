// --------------------------------------------------------------------------------------------------------------------------------

#include "GreensFunction3DAbsSym_test.hpp"
#include "GreensFunction3DAbsSym.hpp"

// --------------------------------------------------------------------------------------------------------------------------------

int GreensFunction3DAbsSym_drawTime_small_aD()
{
   auto D = 1e-12;
   auto a = 1e-7;
   auto gf = GreensFunction3DAbsSym(D, a);
   double rnd_and_time[25][2] = { { 0.315811079434995963262, 0.001866099146214802628 },
                                 { 0.745372736243210645656, 0.000935070166686555458 },
                                 { 0.677021235593504536851, 0.001051735968802082478 },
                                 { 0.806800763907906506931, 0.000830918098511605136 },
                                 { 0.852424862903069594895, 0.000751060369004093120 },
                                 { 0.795064986644005454127, 0.000850969145106009118 },
                                 { 0.953419927844552691524, 0.000537771057223989287 },
                                 { 0.640251668263613308548, 0.001116380730818125579 },
                                 { 0.196608708442741987555, 0.002349368535322374736 },
                                 { 0.003473064802086145951, 0.006439837670894756810 },
                                 { 0.707237250769455252549, 0.000999766304582363884 },
                                 { 0.214137137310948565563, 0.002262555864859691926 },
                                 { 0.258209103153735122826, 0.002071982235742249236 },
                                 { 0.456664658091753616602, 0.001483877157722694958 },
                                 { 0.699517037905452441347, 0.001012967148705036000 },
                                 { 0.787737713343543266272, 0.000863430043287744672 },
                                 { 0.934610618775201494913, 0.000585925896667676615 },
                                 { 0.290369748136950251586, 0.001952108581255179079 },
                                 { 0.178035668385141932446, 0.002450160421808604828 },
                                 { 0.983605578227976891861, 0.000430736023747279411 },
                                 { 0.285096921937765078872, 0.001970845987562706600 },
                                 { 0.498599704289611787933, 0.001390843206697094586 },
                                 { 0.342637646021285668608, 0.001782351739653680852 },
                                 { 0.041373657129610551438, 0.003929488023784270642 },
                                 { 0.331982320943767041474, 0.001814837451512876149 } };

   int ndatapoint = sizeof(rnd_and_time) / sizeof(rnd_and_time[0]);
   for (int index = 0; index < ndatapoint; index++)
   {
      auto time = gf.drawTime(rnd_and_time[index][0]);
      TINYTEST_ALMOST_EQUALS(rnd_and_time[index][1], time);
   }

   return 1;
}

// --------------------------------------------------------------------------------------------------------------------------------

int GreensFunction3DAbsSym_drawTime_large_aD()
{
   auto D = 1e-10;
   auto a = 1e-5;
   auto gf = GreensFunction3DAbsSym(D, a);
   double rnd_and_time[25][2] = { { 0.962716520673271531810, 0.051046489213631109716 },
                                 { 0.748677775299869754975, 0.092948711835637703152 },
                                 { 0.675912593070708313147, 0.105365999075643173474 },
                                 { 0.279714141552093204901, 0.199032523886163768268 },
                                 { 0.210463703690972640572, 0.228015134307036474970 },
                                 { 0.713119821289867938237, 0.098973705255855917295 },
                                 { 0.238038457792912989564, 0.215487283238214162583 },
                                 { 0.593561928506718141462, 0.120159065288583226464 },
                                 { 0.449621580117532932417, 0.150022325652018641677 },
                                 { 0.310360711569103987765, 0.188394669659152199302 },
                                 { 0.165435400409754749267, 0.252467522984175599583 },
                                 { 0.969271356224258461776, 0.048901518437516923975 },
                                 { 0.677704401894629216097, 0.105055103344210919347 },
                                 { 0.246820845780703327091, 0.211796504912558462672 },
                                 { 0.582376507096108170691, 0.122264948877969920841 },
                                 { 0.727524504629876840639, 0.096526715257150930889 },
                                 { 0.881664482690083542078, 0.069700215265919309050 },
                                 { 0.187639837080878714296, 0.239680314167070247817 },
                                 { 0.932150730339122969806, 0.059173095542875944372 },
                                 { 0.315446232755978443773, 0.186728455314237000396 },
                                 { 0.869044562939920289461, 0.072072684918660385113 },
                                 { 0.555890689049327306809, 0.127369052988740635475 },
                                 { 0.865845970056310553586, 0.072663746166881606096 },
                                 { 0.274064908240819005913, 0.201116566050357020205 },
                                 { 0.590142030787426677404, 0.120800002238968141363 } };

   int ndatapoint = sizeof(rnd_and_time) / sizeof(rnd_and_time[0]);
   for (int index = 0; index < ndatapoint; index++)
   {
      auto time = gf.drawTime(rnd_and_time[index][0]);
      TINYTEST_ALMOST_EQUALS(rnd_and_time[index][1], time);
   }

   return 1;
}

// --------------------------------------------------------------------------------------------------------------------------------

int GreensFunction3DAbsSym_drawR_prCuttOffRegion_small_aDt()
{
   auto D = 1e-12;
   auto a = 1e-8;
   auto time = 1e-5;
   auto gf = GreensFunction3DAbsSym(D, a);

   double rnd_and_r[10][2] = { { 0.025888941676681059247, 1.866244264996706452603e-9 },
                              { 0.179480275071015114571, 3.750103912827415496206e-9 },
                              { 0.282370436151638222395, 4.493858332626677474982e-9 },
                              { 0.382221453526523005345, 5.116501067825286874521e-9 },
                              { 0.440786578726513916131, 5.459283298556091837304e-9 },
                              { 0.555538743934046246823, 6.111910415293572856232e-9 },
                              { 0.675856382630283861321, 6.804416958100045087907e-9 },
                              { 0.720552564985632027712, 7.073449451042937295445e-9 },
                              { 0.898213535409815805495, 8.332126801770337646157e-9 },
                              { 0.937039267671582088959, 8.708016637881057504890e-9 } };

   int ndatapoint = sizeof(rnd_and_r) / sizeof(rnd_and_r[0]);
   for (int index = 0; index < ndatapoint; index++)
   {
      auto r = gf.drawR(rnd_and_r[index][0], time);
      TINYTEST_ALMOST_EQUALS(rnd_and_r[index][1], r);
   }

   return 1;
}

// --------------------------------------------------------------------------------------------------------------------------------

int GreensFunction3DAbsSym_drawR_prCuttOffRegion_large_aDt()
{
   auto D = 1e-10;
   auto a = 1e-5;
   auto time = 0.5;
   auto gf = GreensFunction3DAbsSym(D, a);

   double rnd_and_r[10][2] = { { 0.011363384872635650770, 1.523241873280303539606e-6 },
                              { 0.169505874124747309065, 3.918343692226884178983e-6 },
                              { 0.250220523244313384369, 4.544314360787700781926e-6 },
                              { 0.352712214404456419789, 5.214573536179250546283e-6 },
                              { 0.481352243014477146801, 5.962082249528399129215e-6 },
                              { 0.573342075873749849532, 6.470497578347329349624e-6 },
                              { 0.658822613299087747655, 6.943108119763216269618e-6 },
                              { 0.755732097180417307635, 7.501937284102392259794e-6 },
                              { 0.834921673883016806791, 8.004132284544126182648e-6 },
                              { 0.967608269166751911656, 9.163874330319150278776e-6 } };

   int ndatapoint = sizeof(rnd_and_r) / sizeof(rnd_and_r[0]);
   for (int index = 0; index < ndatapoint; index++)
   {
      auto r = gf.drawR(rnd_and_r[index][0], time);
      TINYTEST_ALMOST_EQUALS(rnd_and_r[index][1], r);
   }

   return 1;
}

// --------------------------------------------------------------------------------------------------------------------------------

int GreensFunction3DAbsSym_drawR_pFreeRegion_large_aDt()
{
   auto D = 1e-12;
   auto a = 1e-5;
   auto time = 0.02;
   auto gf = GreensFunction3DAbsSym(D, a);

   double rnd_and_r[25][2] = { { 0.995687669379400700523, 7.254072777131840817540e-7 },
                              { 0.036306784946310520869, 1.058833350602473724249e-7 },
                              { 0.001320129400853341281, 3.421626703304860055866e-8 },
                              { 0.557882030149868634895, 3.279569388953348274740e-7 },
                              { 0.366189389606645409438, 2.618448703251907762065e-7 },
                              { 0.016534055728839962291, 8.051739899256218766683e-8 },
                              { 0.659398309201323681158, 3.661169526347550271103e-7 },
                              { 0.423118393328652720889, 2.813181459069758020946e-7 },
                              { 0.316598183792522345836, 2.445485256264671423999e-7 },
                              { 0.199253645971683120013, 2.002088431218724537828e-7 },
                              { 0.389561485267561618532, 2.698674078609035080961e-7 },
                              { 0.279556539480342358456, 2.312292288191360358704e-7 },
                              { 0.848177718303045369203, 4.599535025410062892951e-7 },
                              { 0.334101943244114140189, 2.507088419970611429919e-7 },
                              { 0.847172798797933644423, 4.592847086775160627780e-7 },
                              { 0.963042564823356208644, 5.826326254754003488306e-7 },
                              { 0.029527911458207471081, 9.846647519128216513995e-8 },
                              { 0.827196806895627848799, 4.466004230352254120243e-7 },
                              { 0.406422502035251381673, 2.756269674209372276570e-7 },
                              { 0.204753473617299102683, 2.024667962441998526799e-7 },
                              { 0.419129267285419939498, 2.799589353141872509002e-7 },
                              { 0.320066781067059329361, 2.457750799582370673658e-7 },
                              { 0.495969182183347888014, 3.062416645418982145505e-7 },
                              { 0.642817143489710992443, 3.595696482102742795381e-7 },
                              { 0.730858390360877387667, 3.964759691369815057391e-7 } };

   int ndatapoint = sizeof(rnd_and_r) / sizeof(rnd_and_r[0]);
   for (int index = 0; index < ndatapoint; index++)
   {
      auto r = gf.drawR(rnd_and_r[index][0], time);
      TINYTEST_ALMOST_EQUALS(rnd_and_r[index][1], r);
   }

   return 1;
}

// --------------------------------------------------------------------------------------------------------------------------------

int GreensFunction3DAbsSym_drawR_pFreeRegion_small_aDt()
{
   auto D = 1e-14;
   auto a = 1e-8;
   auto time = 5e-12;
   auto gf = GreensFunction3DAbsSym(D, a);

   double rnd_and_r[25][2] = { { 0.476496687573469052758, 4.736168411869674033309e-13 },
                              { 0.976471327132704386863, 9.737314409076234641515e-13 },
                              { 0.624549022175279472748, 5.573861919719451714650e-13 },
                              { 0.049192886040670318347, 1.864864218985457086476e-13 },
                              { 0.801908578223165663325, 6.829599426999018466993e-13 },
                              { 0.988828733216949346158, 1.053804242932285534030e-12 },
                              { 0.225583787718236207623, 3.333186017952819208834e-13 },
                              { 0.427031562651046409180, 4.469110056185909348214e-13 },
                              { 0.138795975643732303247, 2.740004451081458012335e-13 },
                              { 0.950976013771376230261, 8.864931838864092011542e-13 },
                              { 0.536256669426121101178, 5.064110656218978834594e-13 },
                              { 0.421151852204202912807, 4.437436306213792208549e-13 },
                              { 0.776006968677008780106, 6.612131816364206679842e-13 },
                              { 0.822219649423775661940, 7.013982598691030073386e-13 },
                              { 0.405959080278658853477, 4.355545332218840404967e-13 },
                              { 0.566229184566119278040, 5.232802764795937074178e-13 },
                              { 0.067153813052748146403, 2.086513325872588876711e-13 },
                              { 0.375570620517567484742, 4.191154642670313058662e-13 },
                              { 0.890865191383022252569, 7.779092016976743468796e-13 },
                              { 0.214136389054974780374, 3.261313679263674009086e-13 },
                              { 0.274379745095616228175, 3.626025540483144061816e-13 },
                              { 0.408870927437179011467, 4.371249382088549445448e-13 },
                              { 0.340942622861992999675, 4.001832535677679414077e-13 },
                              { 0.230266874204663362899, 3.362190330213599116887e-13 },
                              { 0.461243249780868616046, 4.653595849985797557088e-13 } };

   int ndatapoint = sizeof(rnd_and_r) / sizeof(rnd_and_r[0]);
   for (int index = 0; index < ndatapoint; index++)
   {
      auto r = gf.drawR(rnd_and_r[index][0], time);
      TINYTEST_ALMOST_EQUALS(rnd_and_r[index][1], r);
   }

   return 1;
}

// --------------------------------------------------------------------------------------------------------------------------------

class UnitTestGreensFunction3DAbsSym
{
public:
   static int test_p_int_r_is_p_int_r_free_with_large_shell()
   {
      double D = 1e-12;
      double a = 1e-6;
      auto gf = GreensFunction3DAbsSym(D, a);

      double r = 1e-9;
      double t = 1e-6;

      double p = gf.p_int_r(r, t);
      double p_free = gf.p_int_r_free(r, t);
      TINYTEST_ALMOST_EQUAL(p, p_free, 2.5E-12);
      return 1;
   }

   static int test_p_int_r_at_a_is_p_survival()
   {
      double D = 1e-12;
      double a = 1e-8;
      auto gf = GreensFunction3DAbsSym(D, a);

      double t = 1e-5;
      double p1 = gf.p_int_r(a, t);
      double psurv1 = gf.p_survival(t);
      TINYTEST_ALMOST_EQUALS(p1, psurv1);

      t = 1e-3;
      double p2 = gf.p_int_r(a, t);
      double psurv2 = gf.p_survival(t);
      TINYTEST_ALMOST_EQUALS(p2, psurv2);
      return 1;
   }

   //static int test_p_int_r_and_p_r_fourier()
   //{
   //   double D = 1e-12;
   //   double a = 1e-6;
   //   auto gf = GreensFunction3DAbsSym(D, a);

   //   double r = 1e-9;
   //   double t = 1e-6;

   //   double pi = gf.p_int_r(r, t);
   //   double pf = gf.p_r_fourier(r, t);
   //   TINYTEST_ALMOST_EQUALS(pi, pf);
   //   return 1;
   //}

   static int test_ellipticTheta4Zero()
   {
      auto gf = GreensFunction3DAbsSym(0, 0);

      TINYTEST_ALMOST_EQUALS(0.12112420800258050246, gf.ellipticTheta4Zero(0.5));     // et4z(1/2) ~= 0.121124208002580502460849293181867505809858246820960597233�
      TINYTEST_ALMOST_EQUALS(2e-15, (1.0 - gf.ellipticTheta4Zero(1e-15)));      // et4z(1e-15) ~= 1 - 2e-15
      TINYTEST_ALMOST_EQUALS(2e-16, (1.0 - gf.ellipticTheta4Zero(1e-16)));      // et4z(1e-16) ~= 1 - 2e-16
      TINYTEST_ALMOST_EQUALS(2e-17, (1.0 - gf.ellipticTheta4Zero(1e-17)));      // et4z(1e-17) ~= 1 - 2e-17

      return 1;
   }
};

// --------------------------------------------------------------------------------------------------------------------------------

int GreensFunction3DAbsSym_p_int_r_is_p_int_r_free_with_large_shell()
{
   return UnitTestGreensFunction3DAbsSym::test_p_int_r_is_p_int_r_free_with_large_shell();
}

int GreensFunction3DAbsSym_p_int_r_at_a_is_p_survival()
{
   return UnitTestGreensFunction3DAbsSym::test_p_int_r_at_a_is_p_survival();
}

//int GreensFunction3DAbsSym_p_int_r_and_p_r_fourier()
//{
//   return UnitTestGreensFunction3DAbsSym::test_p_int_r_and_p_r_fourier();
//}

int GreensFunction3DAbsSym_ellipticTheta4Zero()
{
   return UnitTestGreensFunction3DAbsSym::test_ellipticTheta4Zero();
}

// --------------------------------------------------------------------------------------------------------------------------------

TINYTEST_START_SUITE(GreensFunction3DAbsSym);
TINYTEST_ADD_TEST(GreensFunction3DAbsSym_drawTime_small_aD);
TINYTEST_ADD_TEST(GreensFunction3DAbsSym_drawTime_large_aD);
TINYTEST_ADD_TEST(GreensFunction3DAbsSym_drawR_prCuttOffRegion_small_aDt);
TINYTEST_ADD_TEST(GreensFunction3DAbsSym_drawR_prCuttOffRegion_large_aDt);
TINYTEST_ADD_TEST(GreensFunction3DAbsSym_drawR_pFreeRegion_large_aDt);
TINYTEST_ADD_TEST(GreensFunction3DAbsSym_drawR_pFreeRegion_small_aDt);
TINYTEST_ADD_TEST(GreensFunction3DAbsSym_p_int_r_is_p_int_r_free_with_large_shell)
TINYTEST_ADD_TEST(GreensFunction3DAbsSym_p_int_r_at_a_is_p_survival)
//TINYTEST_ADD_TEST(GreensFunction3DAbsSym_p_int_r_and_p_r_fourier)
TINYTEST_ADD_TEST(GreensFunction3DAbsSym_ellipticTheta4Zero)
TINYTEST_END_SUITE();

// --------------------------------------------------------------------------------------------------------------------------------
