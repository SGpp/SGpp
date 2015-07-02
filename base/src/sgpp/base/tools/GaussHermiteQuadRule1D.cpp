// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include "GaussHermiteQuadRule1D.hpp"
#include <sgpp/base/exception/factory_exception.hpp>

using namespace SGPP::base;

namespace SGPP {
namespace base {

GaussHermiteQuadRule1D::GaussHermiteQuadRule1D() {
    //------------------------------------------------------------
    // n = 1
    coordinates[0] = new DataVector(1);
    coordinates[0]->set(0, 0.000000000000000);
    //------------------------------------------------------------
    weights[0] = new DataVector(1);
    weights[0]->set(0, 1.000000000000000);
    //------------------------------------------------------------
    // n = 2
    coordinates[1] = new DataVector(2);
    coordinates[1]->set(0, -1.000000000000000);
    coordinates[1]->set(1, 1.000000000000000);
    //------------------------------------------------------------
    weights[1] = new DataVector(2);
    weights[1]->set(0, 0.500000000000000);
    weights[1]->set(1, 0.500000000000000);
    //------------------------------------------------------------
    // n = 3
    coordinates[2] = new DataVector(3);
    coordinates[2]->set(0, -1.732050807568877);
    coordinates[2]->set(1, 0.000000000000000);
    coordinates[2]->set(2, 1.732050807568877);
    //------------------------------------------------------------
    weights[2] = new DataVector(3);
    weights[2]->set(0, 0.166666666666667);
    weights[2]->set(1, 0.666666666666667);
    weights[2]->set(2, 0.166666666666667);
    //------------------------------------------------------------
    // n = 4
    coordinates[3] = new DataVector(4);
    coordinates[3]->set(0, -2.334414218338978);
    coordinates[3]->set(1, -0.741963784302726);
    coordinates[3]->set(2, 0.741963784302726);
    coordinates[3]->set(3, 2.334414218338978);
    //------------------------------------------------------------
    weights[3] = new DataVector(4);
    weights[3]->set(0, 0.045875854768068);
    weights[3]->set(1, 0.454124145231932);
    weights[3]->set(2, 0.454124145231932);
    weights[3]->set(3, 0.045875854768068);
    //------------------------------------------------------------
    // n = 5
    coordinates[4] = new DataVector(5);
    coordinates[4]->set(0, -2.856970013872806);
    coordinates[4]->set(1, -1.355626179974266);
    coordinates[4]->set(2, 0.000000000000000);
    coordinates[4]->set(3, 1.355626179974266);
    coordinates[4]->set(4, 2.856970013872806);
    //------------------------------------------------------------
    weights[4] = new DataVector(5);
    weights[4]->set(0, 0.011257411327721);
    weights[4]->set(1, 0.222075922005613);
    weights[4]->set(2, 0.533333333333333);
    weights[4]->set(3, 0.222075922005613);
    weights[4]->set(4, 0.011257411327721);
    //------------------------------------------------------------
    // n = 6
    coordinates[5] = new DataVector(6);
    coordinates[5]->set(0, -3.324257433552119);
    coordinates[5]->set(1, -1.889175877753711);
    coordinates[5]->set(2, -0.616706590192594);
    coordinates[5]->set(3, 0.616706590192594);
    coordinates[5]->set(4, 1.889175877753711);
    coordinates[5]->set(5, 3.324257433552119);
    //------------------------------------------------------------
    weights[5] = new DataVector(6);
    weights[5]->set(0, 0.002555784402056);
    weights[5]->set(1, 0.088615746041915);
    weights[5]->set(2, 0.408828469556029);
    weights[5]->set(3, 0.408828469556029);
    weights[5]->set(4, 0.088615746041915);
    weights[5]->set(5, 0.002555784402056);
    //------------------------------------------------------------
    // n = 7
    coordinates[6] = new DataVector(7);
    coordinates[6]->set(0, -3.750439717725742);
    coordinates[6]->set(1, -2.366759410734542);
    coordinates[6]->set(2, -1.154405394739968);
    coordinates[6]->set(3, 0.000000000000000);
    coordinates[6]->set(4, 1.154405394739968);
    coordinates[6]->set(5, 2.366759410734542);
    coordinates[6]->set(6, 3.750439717725742);
    //------------------------------------------------------------
    weights[6] = new DataVector(7);
    weights[6]->set(0, 0.000548268855972);
    weights[6]->set(1, 0.030757123967586);
    weights[6]->set(2, 0.240123178605013);
    weights[6]->set(3, 0.457142857142857);
    weights[6]->set(4, 0.240123178605013);
    weights[6]->set(5, 0.030757123967586);
    weights[6]->set(6, 0.000548268855972);
    //------------------------------------------------------------
    // n = 8
    coordinates[7] = new DataVector(8);
    coordinates[7]->set(0, -4.144547186125894);
    coordinates[7]->set(1, -2.802485861287542);
    coordinates[7]->set(2, -1.636519042435108);
    coordinates[7]->set(3, -0.539079811351375);
    coordinates[7]->set(4, 0.539079811351375);
    coordinates[7]->set(5, 1.636519042435108);
    coordinates[7]->set(6, 2.802485861287542);
    coordinates[7]->set(7, 4.144547186125894);
    //------------------------------------------------------------
    weights[7] = new DataVector(8);
    weights[7]->set(0, 0.000112614538375);
    weights[7]->set(1, 0.009635220120788);
    weights[7]->set(2, 0.117239907661759);
    weights[7]->set(3, 0.373012257679077);
    weights[7]->set(4, 0.373012257679077);
    weights[7]->set(5, 0.117239907661759);
    weights[7]->set(6, 0.009635220120788);
    weights[7]->set(7, 0.000112614538375);
    //------------------------------------------------------------
    // n = 9
    coordinates[8] = new DataVector(9);
    coordinates[8]->set(0, -4.512745863399783);
    coordinates[8]->set(1, -3.205429002856470);
    coordinates[8]->set(2, -2.076847978677830);
    coordinates[8]->set(3, -1.023255663789133);
    coordinates[8]->set(4, 0.000000000000000);
    coordinates[8]->set(5, 1.023255663789133);
    coordinates[8]->set(6, 2.076847978677830);
    coordinates[8]->set(7, 3.205429002856470);
    coordinates[8]->set(8, 4.512745863399783);
    //------------------------------------------------------------
    weights[8] = new DataVector(9);
    weights[8]->set(0, 0.000022345844008);
    weights[8]->set(1, 0.002789141321232);
    weights[8]->set(2, 0.049916406765218);
    weights[8]->set(3, 0.244097502894939);
    weights[8]->set(4, 0.406349206349206);
    weights[8]->set(5, 0.244097502894939);
    weights[8]->set(6, 0.049916406765218);
    weights[8]->set(7, 0.002789141321232);
    weights[8]->set(8, 0.000022345844008);
    //------------------------------------------------------------
    // n = 10
    coordinates[9] = new DataVector(10);
    coordinates[9]->set(0, -4.859462828332312);
    coordinates[9]->set(1, -3.581823483551927);
    coordinates[9]->set(2, -2.484325841638955);
    coordinates[9]->set(3, -1.465989094391158);
    coordinates[9]->set(4, -0.484935707515498);
    coordinates[9]->set(5, 0.484935707515498);
    coordinates[9]->set(6, 1.465989094391158);
    coordinates[9]->set(7, 2.484325841638955);
    coordinates[9]->set(8, 3.581823483551927);
    coordinates[9]->set(9, 4.859462828332312);
    //------------------------------------------------------------
    weights[9] = new DataVector(10);
    weights[9]->set(0, 0.000004310652631);
    weights[9]->set(1, 0.000758070934312);
    weights[9]->set(2, 0.019111580500770);
    weights[9]->set(3, 0.135483702980268);
    weights[9]->set(4, 0.344642334932019);
    weights[9]->set(5, 0.344642334932019);
    weights[9]->set(6, 0.135483702980268);
    weights[9]->set(7, 0.019111580500770);
    weights[9]->set(8, 0.000758070934312);
    weights[9]->set(9, 0.000004310652631);
    //------------------------------------------------------------
    // n = 11
    coordinates[10] = new DataVector(11);
    coordinates[10]->set(0, -5.188001224374871);
    coordinates[10]->set(1, -3.936166607129977);
    coordinates[10]->set(2, -2.865123160643646);
    coordinates[10]->set(3, -1.876035020154846);
    coordinates[10]->set(4, -0.928868997381064);
    coordinates[10]->set(5, 0.000000000000000);
    coordinates[10]->set(6, 0.928868997381064);
    coordinates[10]->set(7, 1.876035020154846);
    coordinates[10]->set(8, 2.865123160643646);
    coordinates[10]->set(9, 3.936166607129977);
    coordinates[10]->set(10, 5.188001224374871);
    //------------------------------------------------------------
    weights[10] = new DataVector(11);
    weights[10]->set(0, 0.000000812184979);
    weights[10]->set(1, 0.000195671930271);
    weights[10]->set(2, 0.006720285235537);
    weights[10]->set(3, 0.066138746071058);
    weights[10]->set(4, 0.242240299873970);
    weights[10]->set(5, 0.369408369408370);
    weights[10]->set(6, 0.242240299873970);
    weights[10]->set(7, 0.066138746071058);
    weights[10]->set(8, 0.006720285235537);
    weights[10]->set(9, 0.000195671930271);
    weights[10]->set(10, 0.000000812184979);
    //------------------------------------------------------------
    // n = 12
    coordinates[11] = new DataVector(12);
    coordinates[11]->set(0, -5.500901704467749);
    coordinates[11]->set(1, -4.271825847932282);
    coordinates[11]->set(2, -3.223709828770097);
    coordinates[11]->set(3, -2.259464451000799);
    coordinates[11]->set(4, -1.340375197151617);
    coordinates[11]->set(5, -0.444403001944139);
    coordinates[11]->set(6, 0.444403001944139);
    coordinates[11]->set(7, 1.340375197151617);
    coordinates[11]->set(8, 2.259464451000799);
    coordinates[11]->set(9, 3.223709828770097);
    coordinates[11]->set(10, 4.271825847932282);
    coordinates[11]->set(11, 5.500901704467749);
    //------------------------------------------------------------
    weights[11] = new DataVector(12);
    weights[11]->set(0, 0.000000149992717);
    weights[11]->set(1, 0.000048371849226);
    weights[11]->set(2, 0.002203380687533);
    weights[11]->set(3, 0.029116687912364);
    weights[11]->set(4, 0.146967048045330);
    weights[11]->set(5, 0.321664361512830);
    weights[11]->set(6, 0.321664361512830);
    weights[11]->set(7, 0.146967048045330);
    weights[11]->set(8, 0.029116687912364);
    weights[11]->set(9, 0.002203380687533);
    weights[11]->set(10, 0.000048371849226);
    weights[11]->set(11, 0.000000149992717);
    //------------------------------------------------------------
    // n = 13
    coordinates[12] = new DataVector(13);
    coordinates[12]->set(0, -5.800167252386501);
    coordinates[12]->set(1, -4.591398448936521);
    coordinates[12]->set(2, -3.563444380281635);
    coordinates[12]->set(3, -2.620689973432215);
    coordinates[12]->set(4, -1.725418379588239);
    coordinates[12]->set(5, -0.856679493519450);
    coordinates[12]->set(6, 0.000000000000000);
    coordinates[12]->set(7, 0.856679493519450);
    coordinates[12]->set(8, 1.725418379588239);
    coordinates[12]->set(9, 2.620689973432215);
    coordinates[12]->set(10, 3.563444380281635);
    coordinates[12]->set(11, 4.591398448936521);
    coordinates[12]->set(12, 5.800167252386501);
    //------------------------------------------------------------
    weights[12] = new DataVector(13);
    weights[12]->set(0, 0.000000027226276);
    weights[12]->set(1, 0.000011526596527);
    weights[12]->set(2, 0.000681236350443);
    weights[12]->set(3, 0.011770560505997);
    weights[12]->set(4, 0.079168955860450);
    weights[12]->set(5, 0.237871522964136);
    weights[12]->set(6, 0.340992340992341);
    weights[12]->set(7, 0.237871522964136);
    weights[12]->set(8, 0.079168955860450);
    weights[12]->set(9, 0.011770560505997);
    weights[12]->set(10, 0.000681236350443);
    weights[12]->set(11, 0.000011526596527);
    weights[12]->set(12, 0.000000027226276);
    //------------------------------------------------------------
    // n = 14
    coordinates[13] = new DataVector(14);
    coordinates[13]->set(0, -6.087409546901291);
    coordinates[13]->set(1, -4.896936397345565);
    coordinates[13]->set(2, -3.886924575059770);
    coordinates[13]->set(3, -2.963036579838668);
    coordinates[13]->set(4, -2.088344745701944);
    coordinates[13]->set(5, -1.242688955485464);
    coordinates[13]->set(6, -0.412590457954602);
    coordinates[13]->set(7, 0.412590457954602);
    coordinates[13]->set(8, 1.242688955485464);
    coordinates[13]->set(9, 2.088344745701944);
    coordinates[13]->set(10, 2.963036579838668);
    coordinates[13]->set(11, 3.886924575059770);
    coordinates[13]->set(12, 4.896936397345565);
    coordinates[13]->set(13, 6.087409546901291);
    //------------------------------------------------------------
    weights[13] = new DataVector(14);
    weights[13]->set(0, 0.000000004868161);
    weights[13]->set(1, 0.000002660991344);
    weights[13]->set(2, 0.000200339553761);
    weights[13]->set(3, 0.004428919106947);
    weights[13]->set(4, 0.038650108824253);
    weights[13]->set(5, 0.154083339842514);
    weights[13]->set(6, 0.302634626813020);
    weights[13]->set(7, 0.302634626813020);
    weights[13]->set(8, 0.154083339842514);
    weights[13]->set(9, 0.038650108824253);
    weights[13]->set(10, 0.004428919106947);
    weights[13]->set(11, 0.000200339553761);
    weights[13]->set(12, 0.000002660991344);
    weights[13]->set(13, 0.000000004868161);
    //------------------------------------------------------------
    // n = 15
    coordinates[14] = new DataVector(15);
    coordinates[14]->set(0, -6.363947888829840);
    coordinates[14]->set(1, -5.190093591304782);
    coordinates[14]->set(2, -4.196207711269015);
    coordinates[14]->set(3, -3.289082424398767);
    coordinates[14]->set(4, -2.432436827009758);
    coordinates[14]->set(5, -1.606710069028730);
    coordinates[14]->set(6, -0.799129068324548);
    coordinates[14]->set(7, 0.000000000000000);
    coordinates[14]->set(8, 0.799129068324548);
    coordinates[14]->set(9, 1.606710069028730);
    coordinates[14]->set(10, 2.432436827009758);
    coordinates[14]->set(11, 3.289082424398767);
    coordinates[14]->set(12, 4.196207711269015);
    coordinates[14]->set(13, 5.190093591304782);
    coordinates[14]->set(14, 6.363947888829840);
    //------------------------------------------------------------
    weights[14] = new DataVector(15);
    weights[14]->set(0, 0.000000000858965);
    weights[14]->set(1, 0.000000597541960);
    weights[14]->set(2, 0.000056421464052);
    weights[14]->set(3, 0.001567357503550);
    weights[14]->set(4, 0.017365774492138);
    weights[14]->set(5, 0.089417795399844);
    weights[14]->set(6, 0.232462293609732);
    weights[14]->set(7, 0.318259518259518);
    weights[14]->set(8, 0.232462293609732);
    weights[14]->set(9, 0.089417795399844);
    weights[14]->set(10, 0.017365774492138);
    weights[14]->set(11, 0.001567357503550);
    weights[14]->set(12, 0.000056421464052);
    weights[14]->set(13, 0.000000597541960);
    weights[14]->set(14, 0.000000000858965);
    //------------------------------------------------------------
    // n = 16
    coordinates[15] = new DataVector(16);
    coordinates[15]->set(0, -6.630878198393129);
    coordinates[15]->set(1, -5.472225705949343);
    coordinates[15]->set(2, -4.492955302520012);
    coordinates[15]->set(3, -3.600873624171549);
    coordinates[15]->set(4, -2.760245047630702);
    coordinates[15]->set(5, -1.951980345716334);
    coordinates[15]->set(6, -1.163829100554965);
    coordinates[15]->set(7, -0.386760604500557);
    coordinates[15]->set(8, 0.386760604500557);
    coordinates[15]->set(9, 1.163829100554965);
    coordinates[15]->set(10, 1.951980345716334);
    coordinates[15]->set(11, 2.760245047630702);
    coordinates[15]->set(12, 3.600873624171549);
    coordinates[15]->set(13, 4.492955302520012);
    coordinates[15]->set(14, 5.472225705949343);
    coordinates[15]->set(15, 6.630878198393129);
    //------------------------------------------------------------
    weights[15] = new DataVector(16);
    weights[15]->set(0, 0.000000000149781);
    weights[15]->set(1, 0.000000130947322);
    weights[15]->set(2, 0.000015300032162);
    weights[15]->set(3, 0.000525984926574);
    weights[15]->set(4, 0.007266937601185);
    weights[15]->set(5, 0.047284752354014);
    weights[15]->set(6, 0.158338372750950);
    weights[15]->set(7, 0.286568521238013);
    weights[15]->set(8, 0.286568521238013);
    weights[15]->set(9, 0.158338372750950);
    weights[15]->set(10, 0.047284752354014);
    weights[15]->set(11, 0.007266937601185);
    weights[15]->set(12, 0.000525984926574);
    weights[15]->set(13, 0.000015300032162);
    weights[15]->set(14, 0.000000130947322);
    weights[15]->set(15, 0.000000000149781);
    //------------------------------------------------------------
    // n = 17
    coordinates[16] = new DataVector(17);
    coordinates[16]->set(0, -6.889122439895333);
    coordinates[16]->set(1, -5.744460078659407);
    coordinates[16]->set(2, -4.778531589629984);
    coordinates[16]->set(3, -3.900065717198010);
    coordinates[16]->set(4, -3.073797175328194);
    coordinates[16]->set(5, -2.281019440252989);
    coordinates[16]->set(6, -1.509883307796741);
    coordinates[16]->set(7, -0.751842600703896);
    coordinates[16]->set(8, 0.000000000000000);
    coordinates[16]->set(9, 0.751842600703896);
    coordinates[16]->set(10, 1.509883307796741);
    coordinates[16]->set(11, 2.281019440252989);
    coordinates[16]->set(12, 3.073797175328194);
    coordinates[16]->set(13, 3.900065717198010);
    coordinates[16]->set(14, 4.778531589629984);
    coordinates[16]->set(15, 5.744460078659407);
    coordinates[16]->set(16, 6.889122439895333);
    //------------------------------------------------------------
    weights[16] = new DataVector(17);
    weights[16]->set(0, 0.000000000025843);
    weights[16]->set(1, 0.000000028080161);
    weights[16]->set(2, 0.000004012679448);
    weights[16]->set(3, 0.000168491431551);
    weights[16]->set(4, 0.002858946062285);
    weights[16]->set(5, 0.023086657025711);
    weights[16]->set(6, 0.097406371162721);
    weights[16]->set(7, 0.226706308468977);
    weights[16]->set(8, 0.299538370126606);
    weights[16]->set(9, 0.226706308468977);
    weights[16]->set(10, 0.097406371162721);
    weights[16]->set(11, 0.023086657025711);
    weights[16]->set(12, 0.002858946062285);
    weights[16]->set(13, 0.000168491431551);
    weights[16]->set(14, 0.000004012679448);
    weights[16]->set(15, 0.000000028080161);
    weights[16]->set(16, 0.000000000025843);
    //------------------------------------------------------------
    // n = 18
    coordinates[17] = new DataVector(18);
    coordinates[17]->set(0, -7.139464849146480);
    coordinates[17]->set(1, -6.007745911359598);
    coordinates[17]->set(2, -5.054072685442740);
    coordinates[17]->set(3, -4.188020231629404);
    coordinates[17]->set(4, -3.374736535778092);
    coordinates[17]->set(5, -2.595833688911241);
    coordinates[17]->set(6, -1.839779921508645);
    coordinates[17]->set(7, -1.098395518091501);
    coordinates[17]->set(8, -0.365245755507698);
    coordinates[17]->set(9, 0.365245755507698);
    coordinates[17]->set(10, 1.098395518091501);
    coordinates[17]->set(11, 1.839779921508645);
    coordinates[17]->set(12, 2.595833688911241);
    coordinates[17]->set(13, 3.374736535778092);
    coordinates[17]->set(14, 4.188020231629404);
    coordinates[17]->set(15, 5.054072685442740);
    coordinates[17]->set(16, 6.007745911359598);
    coordinates[17]->set(17, 7.139464849146480);
    //------------------------------------------------------------
    weights[17] = new DataVector(18);
    weights[17]->set(0, 0.000000000004417);
    weights[17]->set(1, 0.000000005905488);
    weights[17]->set(2, 0.000001021552398);
    weights[17]->set(3, 0.000051798961441);
    weights[17]->set(4, 0.001065484796292);
    weights[17]->set(5, 0.010516517751941);
    weights[17]->set(6, 0.054896632480222);
    weights[17]->set(7, 0.160685303893513);
    weights[17]->set(8, 0.272783234654288);
    weights[17]->set(9, 0.272783234654288);
    weights[17]->set(10, 0.160685303893513);
    weights[17]->set(11, 0.054896632480222);
    weights[17]->set(12, 0.010516517751941);
    weights[17]->set(13, 0.001065484796292);
    weights[17]->set(14, 0.000051798961441);
    weights[17]->set(15, 0.000001021552398);
    weights[17]->set(16, 0.000000005905488);
    weights[17]->set(17, 0.000000000004417);
    //------------------------------------------------------------
    // n = 19
    coordinates[18] = new DataVector(19);
    coordinates[18]->set(0, -7.382579024030432);
    coordinates[18]->set(1, -6.262891156513253);
    coordinates[18]->set(2, -5.320536377336039);
    coordinates[18]->set(3, -4.465872626831032);
    coordinates[18]->set(4, -3.664416547450639);
    coordinates[18]->set(5, -2.898051276515754);
    coordinates[18]->set(6, -2.155502761316936);
    coordinates[18]->set(7, -1.428876676078373);
    coordinates[18]->set(8, -0.712085044042380);
    coordinates[18]->set(9, 0.000000000000000);
    coordinates[18]->set(10, 0.712085044042380);
    coordinates[18]->set(11, 1.428876676078373);
    coordinates[18]->set(12, 2.155502761316936);
    coordinates[18]->set(13, 2.898051276515754);
    coordinates[18]->set(14, 3.664416547450639);
    coordinates[18]->set(15, 4.465872626831032);
    coordinates[18]->set(16, 5.320536377336039);
    coordinates[18]->set(17, 6.262891156513253);
    coordinates[18]->set(18, 7.382579024030432);
    //------------------------------------------------------------
    weights[18] = new DataVector(19);
    weights[18]->set(0, 0.000000000000748);
    weights[18]->set(1, 0.000000001220371);
    weights[18]->set(2, 0.000000253222003);
    weights[18]->set(3, 0.000015351145955);
    weights[18]->set(4, 0.000378502109414);
    weights[18]->set(5, 0.004507235420342);
    weights[18]->set(6, 0.028666691030118);
    weights[18]->set(7, 0.103603657276144);
    weights[18]->set(8, 0.220941712199144);
    weights[18]->set(9, 0.283773192751521);
    weights[18]->set(10, 0.220941712199144);
    weights[18]->set(11, 0.103603657276144);
    weights[18]->set(12, 0.028666691030118);
    weights[18]->set(13, 0.004507235420342);
    weights[18]->set(14, 0.000378502109414);
    weights[18]->set(15, 0.000015351145955);
    weights[18]->set(16, 0.000000253222003);
    weights[18]->set(17, 0.000000001220371);
    weights[18]->set(18, 0.000000000000748);
    //------------------------------------------------------------
    // n = 20
    coordinates[19] = new DataVector(20);
    coordinates[19]->set(0, -7.619048541679759);
    coordinates[19]->set(1, -6.510590157013655);
    coordinates[19]->set(2, -5.578738805893201);
    coordinates[19]->set(3, -4.734581334046055);
    coordinates[19]->set(4, -3.943967350657316);
    coordinates[19]->set(5, -3.189014816553390);
    coordinates[19]->set(6, -2.458663611172368);
    coordinates[19]->set(7, -1.745247320814127);
    coordinates[19]->set(8, -1.042945348802751);
    coordinates[19]->set(9, -0.346964157081356);
    coordinates[19]->set(10, 0.346964157081356);
    coordinates[19]->set(11, 1.042945348802751);
    coordinates[19]->set(12, 1.745247320814127);
    coordinates[19]->set(13, 2.458663611172368);
    coordinates[19]->set(14, 3.189014816553390);
    coordinates[19]->set(15, 3.943967350657316);
    coordinates[19]->set(16, 4.734581334046055);
    coordinates[19]->set(17, 5.578738805893201);
    coordinates[19]->set(18, 6.510590157013655);
    coordinates[19]->set(19, 7.619048541679759);
    //------------------------------------------------------------
    weights[19] = new DataVector(20);
    weights[19]->set(0, 0.000000000000126);
    weights[19]->set(1, 0.000000000248206);
    weights[19]->set(2, 0.000000061274903);
    weights[19]->set(3, 0.000004402121090);
    weights[19]->set(4, 0.000128826279962);
    weights[19]->set(5, 0.001830103131080);
    weights[19]->set(6, 0.013997837447101);
    weights[19]->set(7, 0.061506372063977);
    weights[19]->set(8, 0.161739333984000);
    weights[19]->set(9, 0.260793063449555);
    weights[19]->set(10, 0.260793063449555);
    weights[19]->set(11, 0.161739333984000);
    weights[19]->set(12, 0.061506372063977);
    weights[19]->set(13, 0.013997837447101);
    weights[19]->set(14, 0.001830103131080);
    weights[19]->set(15, 0.000128826279962);
    weights[19]->set(16, 0.000004402121090);
    weights[19]->set(17, 0.000000061274903);
    weights[19]->set(18, 0.000000000248206);
    weights[19]->set(19, 0.000000000000126);
}

GaussHermiteQuadRule1D::~GaussHermiteQuadRule1D() {
}

void GaussHermiteQuadRule1D::getLevelPointsAndWeightsNormalized(size_t level,
        DataVector& pcoordinates, DataVector& pweights, float_t mean,
        float_t stdd) {
    if (level < 1 || level > maxSupportedLevel) {
        throw factory_exception(
                "GaussHermiteQuadRule1D::getLevelPointsAndWeightsNormalized : order of gauss quadrature has to be within {1, ..., 20}");
    }

    getLevelPointsAndWeights(level, pcoordinates, pweights);

    // scale coordiantes
    for (size_t i = 0; i < level; i++) {
        pcoordinates[i] = (pcoordinates[i] + mean) * stdd;
    }
}

} /* namespace base */
} /* namespace SGPP */
