// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include "GaussLegendreQuadRule1D.hpp"
#include <sgpp/base/exception/factory_exception.hpp>

#include <iostream>

using namespace SGPP::base;

namespace SGPP {
namespace base {

GaussLegendreQuadRule1D::GaussLegendreQuadRule1D() {
    //------------------------------------------------------------
    // n = 1
    coordinates[0] = new DataVector(1);
    coordinates[0]->set(0, 0.000000000000000);
    //------------------------------------------------------------
    weights[0] = new DataVector(1);
    weights[0]->set(0, 2.000000000000000);
    //------------------------------------------------------------
    // n = 2
    coordinates[1] = new DataVector(2);
    coordinates[1]->set(0, -0.577350269189626);
    coordinates[1]->set(1, 0.577350269189626);
    //------------------------------------------------------------
    weights[1] = new DataVector(2);
    weights[1]->set(0, 1.000000000000000);
    weights[1]->set(1, 1.000000000000000);
    //------------------------------------------------------------
    // n = 3
    coordinates[2] = new DataVector(3);
    coordinates[2]->set(0, -0.774596669241483);
    coordinates[2]->set(1, 0.000000000000000);
    coordinates[2]->set(2, 0.774596669241483);
    //------------------------------------------------------------
    weights[2] = new DataVector(3);
    weights[2]->set(0, 0.555555555555555);
    weights[2]->set(1, 0.888888888888889);
    weights[2]->set(2, 0.555555555555555);
    //------------------------------------------------------------
    // n = 4
    coordinates[3] = new DataVector(4);
    coordinates[3]->set(0, -0.861136311594053);
    coordinates[3]->set(1, -0.339981043584856);
    coordinates[3]->set(2, 0.339981043584856);
    coordinates[3]->set(3, 0.861136311594053);
    //------------------------------------------------------------
    weights[3] = new DataVector(4);
    weights[3]->set(0, 0.347854845137454);
    weights[3]->set(1, 0.652145154862546);
    weights[3]->set(2, 0.652145154862546);
    weights[3]->set(3, 0.347854845137454);
    //------------------------------------------------------------
    // n = 5
    coordinates[4] = new DataVector(5);
    coordinates[4]->set(0, -0.906179845938664);
    coordinates[4]->set(1, -0.538469310105683);
    coordinates[4]->set(2, 0.000000000000000);
    coordinates[4]->set(3, 0.538469310105683);
    coordinates[4]->set(4, 0.906179845938664);
    //------------------------------------------------------------
    weights[4] = new DataVector(5);
    weights[4]->set(0, 0.236926885056189);
    weights[4]->set(1, 0.478628670499367);
    weights[4]->set(2, 0.568888888888889);
    weights[4]->set(3, 0.478628670499367);
    weights[4]->set(4, 0.236926885056189);
    //------------------------------------------------------------
    // n = 6
    coordinates[5] = new DataVector(6);
    coordinates[5]->set(0, -0.932469514203152);
    coordinates[5]->set(1, -0.661209386466264);
    coordinates[5]->set(2, -0.238619186083197);
    coordinates[5]->set(3, 0.238619186083197);
    coordinates[5]->set(4, 0.661209386466264);
    coordinates[5]->set(5, 0.932469514203152);
    //------------------------------------------------------------
    weights[5] = new DataVector(6);
    weights[5]->set(0, 0.171324492379170);
    weights[5]->set(1, 0.360761573048139);
    weights[5]->set(2, 0.467913934572691);
    weights[5]->set(3, 0.467913934572691);
    weights[5]->set(4, 0.360761573048139);
    weights[5]->set(5, 0.171324492379170);
    //------------------------------------------------------------
    // n = 7
    coordinates[6] = new DataVector(7);
    coordinates[6]->set(0, -0.949107912342759);
    coordinates[6]->set(1, -0.741531185599394);
    coordinates[6]->set(2, -0.405845151377397);
    coordinates[6]->set(3, 0.000000000000000);
    coordinates[6]->set(4, 0.405845151377397);
    coordinates[6]->set(5, 0.741531185599394);
    coordinates[6]->set(6, 0.949107912342759);
    //------------------------------------------------------------
    weights[6] = new DataVector(7);
    weights[6]->set(0, 0.129484966168869);
    weights[6]->set(1, 0.279705391489277);
    weights[6]->set(2, 0.381830050505119);
    weights[6]->set(3, 0.417959183673470);
    weights[6]->set(4, 0.381830050505119);
    weights[6]->set(5, 0.279705391489277);
    weights[6]->set(6, 0.129484966168869);
    //------------------------------------------------------------
    // n = 8
    coordinates[7] = new DataVector(8);
    coordinates[7]->set(0, -0.960289856497536);
    coordinates[7]->set(1, -0.796666477413627);
    coordinates[7]->set(2, -0.525532409916329);
    coordinates[7]->set(3, -0.183434642495650);
    coordinates[7]->set(4, 0.183434642495650);
    coordinates[7]->set(5, 0.525532409916329);
    coordinates[7]->set(6, 0.796666477413627);
    coordinates[7]->set(7, 0.960289856497536);
    //------------------------------------------------------------
    weights[7] = new DataVector(8);
    weights[7]->set(0, 0.101228536290376);
    weights[7]->set(1, 0.222381034453374);
    weights[7]->set(2, 0.313706645877888);
    weights[7]->set(3, 0.362683783378362);
    weights[7]->set(4, 0.362683783378362);
    weights[7]->set(5, 0.313706645877888);
    weights[7]->set(6, 0.222381034453374);
    weights[7]->set(7, 0.101228536290376);
    //------------------------------------------------------------
    // n = 9
    coordinates[8] = new DataVector(9);
    coordinates[8]->set(0, -0.968160239507626);
    coordinates[8]->set(1, -0.836031107326636);
    coordinates[8]->set(2, -0.613371432700590);
    coordinates[8]->set(3, -0.324253423403809);
    coordinates[8]->set(4, 0.000000000000000);
    coordinates[8]->set(5, 0.324253423403809);
    coordinates[8]->set(6, 0.613371432700590);
    coordinates[8]->set(7, 0.836031107326636);
    coordinates[8]->set(8, 0.968160239507626);
    //------------------------------------------------------------
    weights[8] = new DataVector(9);
    weights[8]->set(0, 0.081274388361575);
    weights[8]->set(1, 0.180648160694857);
    weights[8]->set(2, 0.260610696402935);
    weights[8]->set(3, 0.312347077040002);
    weights[8]->set(4, 0.330239355001259);
    weights[8]->set(5, 0.312347077040002);
    weights[8]->set(6, 0.260610696402935);
    weights[8]->set(7, 0.180648160694857);
    weights[8]->set(8, 0.081274388361575);
    //------------------------------------------------------------
    // n = 10
    coordinates[9] = new DataVector(10);
    coordinates[9]->set(0, -0.973906528517172);
    coordinates[9]->set(1, -0.865063366688985);
    coordinates[9]->set(2, -0.679409568299024);
    coordinates[9]->set(3, -0.433395394129247);
    coordinates[9]->set(4, -0.148874338981631);
    coordinates[9]->set(5, 0.148874338981631);
    coordinates[9]->set(6, 0.433395394129247);
    coordinates[9]->set(7, 0.679409568299024);
    coordinates[9]->set(8, 0.865063366688985);
    coordinates[9]->set(9, 0.973906528517172);
    //------------------------------------------------------------
    weights[9] = new DataVector(10);
    weights[9]->set(0, 0.066671344308687);
    weights[9]->set(1, 0.149451349150581);
    weights[9]->set(2, 0.219086362515983);
    weights[9]->set(3, 0.269266719309997);
    weights[9]->set(4, 0.295524224714753);
    weights[9]->set(5, 0.295524224714753);
    weights[9]->set(6, 0.269266719309997);
    weights[9]->set(7, 0.219086362515983);
    weights[9]->set(8, 0.149451349150581);
    weights[9]->set(9, 0.066671344308687);
    //------------------------------------------------------------
    // n = 11
    coordinates[10] = new DataVector(11);
    coordinates[10]->set(0, -0.978228658146057);
    coordinates[10]->set(1, -0.887062599768095);
    coordinates[10]->set(2, -0.730152005574049);
    coordinates[10]->set(3, -0.519096129206812);
    coordinates[10]->set(4, -0.269543155952345);
    coordinates[10]->set(5, 0.000000000000000);
    coordinates[10]->set(6, 0.269543155952345);
    coordinates[10]->set(7, 0.519096129206812);
    coordinates[10]->set(8, 0.730152005574049);
    coordinates[10]->set(9, 0.887062599768095);
    coordinates[10]->set(10, 0.978228658146057);
    //------------------------------------------------------------
    weights[10] = new DataVector(11);
    weights[10]->set(0, 0.055668567116173);
    weights[10]->set(1, 0.125580369464905);
    weights[10]->set(2, 0.186290210927734);
    weights[10]->set(3, 0.233193764591990);
    weights[10]->set(4, 0.262804544510247);
    weights[10]->set(5, 0.272925086777901);
    weights[10]->set(6, 0.262804544510247);
    weights[10]->set(7, 0.233193764591990);
    weights[10]->set(8, 0.186290210927734);
    weights[10]->set(9, 0.125580369464905);
    weights[10]->set(10, 0.055668567116173);
    //------------------------------------------------------------
    // n = 12
    coordinates[11] = new DataVector(12);
    coordinates[11]->set(0, -0.981560634246719);
    coordinates[11]->set(1, -0.904117256370475);
    coordinates[11]->set(2, -0.769902674194305);
    coordinates[11]->set(3, -0.587317954286617);
    coordinates[11]->set(4, -0.367831498998180);
    coordinates[11]->set(5, -0.125233408511469);
    coordinates[11]->set(6, 0.125233408511469);
    coordinates[11]->set(7, 0.367831498998180);
    coordinates[11]->set(8, 0.587317954286617);
    coordinates[11]->set(9, 0.769902674194305);
    coordinates[11]->set(10, 0.904117256370475);
    coordinates[11]->set(11, 0.981560634246719);
    //------------------------------------------------------------
    weights[11] = new DataVector(12);
    weights[11]->set(0, 0.047175336386511);
    weights[11]->set(1, 0.106939325995319);
    weights[11]->set(2, 0.160078328543346);
    weights[11]->set(3, 0.203167426723066);
    weights[11]->set(4, 0.233492536538355);
    weights[11]->set(5, 0.249147045813403);
    weights[11]->set(6, 0.249147045813403);
    weights[11]->set(7, 0.233492536538355);
    weights[11]->set(8, 0.203167426723066);
    weights[11]->set(9, 0.160078328543346);
    weights[11]->set(10, 0.106939325995319);
    weights[11]->set(11, 0.047175336386511);
    //------------------------------------------------------------
    // n = 13
    coordinates[12] = new DataVector(13);
    coordinates[12]->set(0, -0.984183054718588);
    coordinates[12]->set(1, -0.917598399222978);
    coordinates[12]->set(2, -0.801578090733310);
    coordinates[12]->set(3, -0.642349339440340);
    coordinates[12]->set(4, -0.448492751036447);
    coordinates[12]->set(5, -0.230458315955135);
    coordinates[12]->set(6, 0.000000000000000);
    coordinates[12]->set(7, 0.230458315955135);
    coordinates[12]->set(8, 0.448492751036447);
    coordinates[12]->set(9, 0.642349339440340);
    coordinates[12]->set(10, 0.801578090733310);
    coordinates[12]->set(11, 0.917598399222978);
    coordinates[12]->set(12, 0.984183054718588);
    //------------------------------------------------------------
    weights[12] = new DataVector(13);
    weights[12]->set(0, 0.040484004765313);
    weights[12]->set(1, 0.092121499837729);
    weights[12]->set(2, 0.138873510219788);
    weights[12]->set(3, 0.178145980761946);
    weights[12]->set(4, 0.207816047536889);
    weights[12]->set(5, 0.226283180262898);
    weights[12]->set(6, 0.232551553230875);
    weights[12]->set(7, 0.226283180262898);
    weights[12]->set(8, 0.207816047536889);
    weights[12]->set(9, 0.178145980761946);
    weights[12]->set(10, 0.138873510219788);
    weights[12]->set(11, 0.092121499837729);
    weights[12]->set(12, 0.040484004765313);
    //------------------------------------------------------------
    // n = 14
    coordinates[13] = new DataVector(14);
    coordinates[13]->set(0, -0.986283808696812);
    coordinates[13]->set(1, -0.928434883663574);
    coordinates[13]->set(2, -0.827201315069765);
    coordinates[13]->set(3, -0.687292904811685);
    coordinates[13]->set(4, -0.515248636358154);
    coordinates[13]->set(5, -0.319112368927890);
    coordinates[13]->set(6, -0.108054948707344);
    coordinates[13]->set(7, 0.108054948707344);
    coordinates[13]->set(8, 0.319112368927890);
    coordinates[13]->set(9, 0.515248636358154);
    coordinates[13]->set(10, 0.687292904811685);
    coordinates[13]->set(11, 0.827201315069765);
    coordinates[13]->set(12, 0.928434883663574);
    coordinates[13]->set(13, 0.986283808696812);
    //------------------------------------------------------------
    weights[13] = new DataVector(14);
    weights[13]->set(0, 0.035119460331756);
    weights[13]->set(1, 0.080158087159761);
    weights[13]->set(2, 0.121518570687902);
    weights[13]->set(3, 0.157203167158193);
    weights[13]->set(4, 0.185538397477937);
    weights[13]->set(5, 0.205198463721295);
    weights[13]->set(6, 0.215263853463157);
    weights[13]->set(7, 0.215263853463157);
    weights[13]->set(8, 0.205198463721295);
    weights[13]->set(9, 0.185538397477937);
    weights[13]->set(10, 0.157203167158193);
    weights[13]->set(11, 0.121518570687902);
    weights[13]->set(12, 0.080158087159761);
    weights[13]->set(13, 0.035119460331756);
    //------------------------------------------------------------
    // n = 15
    coordinates[14] = new DataVector(15);
    coordinates[14]->set(0, -0.987992518020485);
    coordinates[14]->set(1, -0.937273392400706);
    coordinates[14]->set(2, -0.848206583410427);
    coordinates[14]->set(3, -0.724417731360170);
    coordinates[14]->set(4, -0.570972172608539);
    coordinates[14]->set(5, -0.394151347077563);
    coordinates[14]->set(6, -0.201194093997435);
    coordinates[14]->set(7, 0.000000000000000);
    coordinates[14]->set(8, 0.201194093997435);
    coordinates[14]->set(9, 0.394151347077563);
    coordinates[14]->set(10, 0.570972172608539);
    coordinates[14]->set(11, 0.724417731360170);
    coordinates[14]->set(12, 0.848206583410427);
    coordinates[14]->set(13, 0.937273392400706);
    coordinates[14]->set(14, 0.987992518020485);
    //------------------------------------------------------------
    weights[14] = new DataVector(15);
    weights[14]->set(0, 0.030753241996118);
    weights[14]->set(1, 0.070366047488108);
    weights[14]->set(2, 0.107159220467172);
    weights[14]->set(3, 0.139570677926154);
    weights[14]->set(4, 0.166269205816994);
    weights[14]->set(5, 0.186161000015562);
    weights[14]->set(6, 0.198431485327112);
    weights[14]->set(7, 0.202578241925561);
    weights[14]->set(8, 0.198431485327112);
    weights[14]->set(9, 0.186161000015562);
    weights[14]->set(10, 0.166269205816994);
    weights[14]->set(11, 0.139570677926154);
    weights[14]->set(12, 0.107159220467172);
    weights[14]->set(13, 0.070366047488108);
    weights[14]->set(14, 0.030753241996118);
    //------------------------------------------------------------
    // n = 16
    coordinates[15] = new DataVector(16);
    coordinates[15]->set(0, -0.989400934991650);
    coordinates[15]->set(1, -0.944575023073233);
    coordinates[15]->set(2, -0.865631202387832);
    coordinates[15]->set(3, -0.755404408355003);
    coordinates[15]->set(4, -0.617876244402644);
    coordinates[15]->set(5, -0.458016777657227);
    coordinates[15]->set(6, -0.281603550779259);
    coordinates[15]->set(7, -0.095012509837637);
    coordinates[15]->set(8, 0.095012509837637);
    coordinates[15]->set(9, 0.281603550779259);
    coordinates[15]->set(10, 0.458016777657227);
    coordinates[15]->set(11, 0.617876244402644);
    coordinates[15]->set(12, 0.755404408355003);
    coordinates[15]->set(13, 0.865631202387832);
    coordinates[15]->set(14, 0.944575023073233);
    coordinates[15]->set(15, 0.989400934991650);
    //------------------------------------------------------------
    weights[15] = new DataVector(16);
    weights[15]->set(0, 0.027152459411758);
    weights[15]->set(1, 0.062253523938649);
    weights[15]->set(2, 0.095158511682492);
    weights[15]->set(3, 0.124628971255533);
    weights[15]->set(4, 0.149595988816576);
    weights[15]->set(5, 0.169156519395002);
    weights[15]->set(6, 0.182603415044923);
    weights[15]->set(7, 0.189450610455067);
    weights[15]->set(8, 0.189450610455067);
    weights[15]->set(9, 0.182603415044923);
    weights[15]->set(10, 0.169156519395002);
    weights[15]->set(11, 0.149595988816576);
    weights[15]->set(12, 0.124628971255533);
    weights[15]->set(13, 0.095158511682492);
    weights[15]->set(14, 0.062253523938649);
    weights[15]->set(15, 0.027152459411758);
    //------------------------------------------------------------
    // n = 17
    coordinates[16] = new DataVector(17);
    coordinates[16]->set(0, -0.990575475314417);
    coordinates[16]->set(1, -0.950675521768768);
    coordinates[16]->set(2, -0.880239153726986);
    coordinates[16]->set(3, -0.781514003896801);
    coordinates[16]->set(4, -0.657671159216691);
    coordinates[16]->set(5, -0.512690537086477);
    coordinates[16]->set(6, -0.351231763453876);
    coordinates[16]->set(7, -0.178484181495848);
    coordinates[16]->set(8, 0.000000000000000);
    coordinates[16]->set(9, 0.178484181495848);
    coordinates[16]->set(10, 0.351231763453876);
    coordinates[16]->set(11, 0.512690537086477);
    coordinates[16]->set(12, 0.657671159216691);
    coordinates[16]->set(13, 0.781514003896801);
    coordinates[16]->set(14, 0.880239153726986);
    coordinates[16]->set(15, 0.950675521768768);
    coordinates[16]->set(16, 0.990575475314417);
    //------------------------------------------------------------
    weights[16] = new DataVector(17);
    weights[16]->set(0, 0.024148302868548);
    weights[16]->set(1, 0.055459529373986);
    weights[16]->set(2, 0.085036148317179);
    weights[16]->set(3, 0.111883847193404);
    weights[16]->set(4, 0.135136368468526);
    weights[16]->set(5, 0.154045761076811);
    weights[16]->set(6, 0.168004102156450);
    weights[16]->set(7, 0.176562705366993);
    weights[16]->set(8, 0.179446470356207);
    weights[16]->set(9, 0.176562705366993);
    weights[16]->set(10, 0.168004102156450);
    weights[16]->set(11, 0.154045761076811);
    weights[16]->set(12, 0.135136368468526);
    weights[16]->set(13, 0.111883847193404);
    weights[16]->set(14, 0.085036148317179);
    weights[16]->set(15, 0.055459529373986);
    weights[16]->set(16, 0.024148302868548);
    //------------------------------------------------------------
    // n = 18
    coordinates[17] = new DataVector(18);
    coordinates[17]->set(0, -0.991565168420931);
    coordinates[17]->set(1, -0.955823949571398);
    coordinates[17]->set(2, -0.892602466497556);
    coordinates[17]->set(3, -0.803704958972523);
    coordinates[17]->set(4, -0.691687043060353);
    coordinates[17]->set(5, -0.559770831073948);
    coordinates[17]->set(6, -0.411751161462843);
    coordinates[17]->set(7, -0.251886225691505);
    coordinates[17]->set(8, -0.084775013041735);
    coordinates[17]->set(9, 0.084775013041735);
    coordinates[17]->set(10, 0.251886225691505);
    coordinates[17]->set(11, 0.411751161462843);
    coordinates[17]->set(12, 0.559770831073948);
    coordinates[17]->set(13, 0.691687043060353);
    coordinates[17]->set(14, 0.803704958972523);
    coordinates[17]->set(15, 0.892602466497556);
    coordinates[17]->set(16, 0.955823949571398);
    coordinates[17]->set(17, 0.991565168420931);
    //------------------------------------------------------------
    weights[17] = new DataVector(18);
    weights[17]->set(0, 0.021616013526484);
    weights[17]->set(1, 0.049714548894969);
    weights[17]->set(2, 0.076425730254890);
    weights[17]->set(3, 0.100942044106287);
    weights[17]->set(4, 0.122555206711478);
    weights[17]->set(5, 0.140642914670651);
    weights[17]->set(6, 0.154684675126265);
    weights[17]->set(7, 0.164276483745833);
    weights[17]->set(8, 0.169142382963144);
    weights[17]->set(9, 0.169142382963144);
    weights[17]->set(10, 0.164276483745833);
    weights[17]->set(11, 0.154684675126265);
    weights[17]->set(12, 0.140642914670651);
    weights[17]->set(13, 0.122555206711478);
    weights[17]->set(14, 0.100942044106287);
    weights[17]->set(15, 0.076425730254890);
    weights[17]->set(16, 0.049714548894969);
    weights[17]->set(17, 0.021616013526484);
    //------------------------------------------------------------
    // n = 19
    coordinates[18] = new DataVector(19);
    coordinates[18]->set(0, -0.992406843843584);
    coordinates[18]->set(1, -0.960208152134830);
    coordinates[18]->set(2, -0.903155903614818);
    coordinates[18]->set(3, -0.822714656537143);
    coordinates[18]->set(4, -0.720966177335229);
    coordinates[18]->set(5, -0.600545304661681);
    coordinates[18]->set(6, -0.464570741375961);
    coordinates[18]->set(7, -0.316564099963630);
    coordinates[18]->set(8, -0.160358645640225);
    coordinates[18]->set(9, 0.000000000000000);
    coordinates[18]->set(10, 0.160358645640225);
    coordinates[18]->set(11, 0.316564099963630);
    coordinates[18]->set(12, 0.464570741375961);
    coordinates[18]->set(13, 0.600545304661681);
    coordinates[18]->set(14, 0.720966177335229);
    coordinates[18]->set(15, 0.822714656537143);
    coordinates[18]->set(16, 0.903155903614818);
    coordinates[18]->set(17, 0.960208152134830);
    coordinates[18]->set(18, 0.992406843843584);
    //------------------------------------------------------------
    weights[18] = new DataVector(19);
    weights[18]->set(0, 0.019461788229726);
    weights[18]->set(1, 0.044814226765701);
    weights[18]->set(2, 0.069044542737641);
    weights[18]->set(3, 0.091490021622450);
    weights[18]->set(4, 0.111566645547334);
    weights[18]->set(5, 0.128753962539336);
    weights[18]->set(6, 0.142606702173607);
    weights[18]->set(7, 0.152766042065860);
    weights[18]->set(8, 0.158968843393954);
    weights[18]->set(9, 0.161054449848784);
    weights[18]->set(10, 0.158968843393954);
    weights[18]->set(11, 0.152766042065860);
    weights[18]->set(12, 0.142606702173607);
    weights[18]->set(13, 0.128753962539336);
    weights[18]->set(14, 0.111566645547334);
    weights[18]->set(15, 0.091490021622450);
    weights[18]->set(16, 0.069044542737641);
    weights[18]->set(17, 0.044814226765701);
    weights[18]->set(18, 0.019461788229726);
    //------------------------------------------------------------
    // n = 20
    coordinates[19] = new DataVector(20);
    coordinates[19]->set(0, -0.993128599185095);
    coordinates[19]->set(1, -0.963971927277914);
    coordinates[19]->set(2, -0.912234428251326);
    coordinates[19]->set(3, -0.839116971822219);
    coordinates[19]->set(4, -0.746331906460151);
    coordinates[19]->set(5, -0.636053680726515);
    coordinates[19]->set(6, -0.510867001950827);
    coordinates[19]->set(7, -0.373706088715420);
    coordinates[19]->set(8, -0.227785851141645);
    coordinates[19]->set(9, -0.076526521133497);
    coordinates[19]->set(10, 0.076526521133497);
    coordinates[19]->set(11, 0.227785851141645);
    coordinates[19]->set(12, 0.373706088715420);
    coordinates[19]->set(13, 0.510867001950827);
    coordinates[19]->set(14, 0.636053680726515);
    coordinates[19]->set(15, 0.746331906460151);
    coordinates[19]->set(16, 0.839116971822219);
    coordinates[19]->set(17, 0.912234428251326);
    coordinates[19]->set(18, 0.963971927277914);
    coordinates[19]->set(19, 0.993128599185095);
    //------------------------------------------------------------
    weights[19] = new DataVector(20);
    weights[19]->set(0, 0.017614007139150);
    weights[19]->set(1, 0.040601429800387);
    weights[19]->set(2, 0.062672048334110);
    weights[19]->set(3, 0.083276741576705);
    weights[19]->set(4, 0.101930119817240);
    weights[19]->set(5, 0.118194531961519);
    weights[19]->set(6, 0.131688638449177);
    weights[19]->set(7, 0.142096109318382);
    weights[19]->set(8, 0.149172986472604);
    weights[19]->set(9, 0.152753387130726);
    weights[19]->set(10, 0.152753387130726);
    weights[19]->set(11, 0.149172986472604);
    weights[19]->set(12, 0.142096109318382);
    weights[19]->set(13, 0.131688638449177);
    weights[19]->set(14, 0.118194531961519);
    weights[19]->set(15, 0.101930119817240);
    weights[19]->set(16, 0.083276741576705);
    weights[19]->set(17, 0.062672048334110);
    weights[19]->set(18, 0.040601429800387);
    weights[19]->set(19, 0.017614007139150);
}

GaussLegendreQuadRule1D::~GaussLegendreQuadRule1D() {
}

void GaussLegendreQuadRule1D::getLevelPointsAndWeightsNormalized(size_t level,
        DataVector& pcoordinates, DataVector& pweights) {
    if (level < 1 || level > maxSupportedLevel) {
        throw factory_exception(
                "GaussLegendreQuadRule1D::getLevelPointsAndWeightsNormalized : order of gauss quadrature has to be within {1, ..., 20}");
    }

    getLevelPointsAndWeights(level, pcoordinates, pweights);

    // scale coordinates
    for (size_t i = 0; i < level; i++) {
        // [-1, 1] -> [0, 1]
        pcoordinates[i] = (pcoordinates[i] + 1.0f) / 2.0f;
    }

    // transform weights according to the volume of the linear transformation
    pweights.mult(0.5f);
}

} /* namespace base */
} /* namespace SGPP */
