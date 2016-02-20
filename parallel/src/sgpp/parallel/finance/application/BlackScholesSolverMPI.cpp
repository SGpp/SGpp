// Copyright (C) 2008-today The SG++ Project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/parallel/tools/MPI/SGppMPITools.hpp>
#include <sgpp/parallel/solver/sle/BiCGStabMPI.hpp>
#include <sgpp/parallel/solver/sle/ConjugateGradientsMPI.hpp>
#include <sgpp/parallel/finance/algorithm/BlackScholesParabolicPDESolverSystemEuroAmerParallelMPI.hpp>
#include <sgpp/parallel/finance/algorithm/BlackScholesPATParabolicPDESolverSystemEuroAmerParallelMPI.hpp>
#include <sgpp/parallel/finance/application/BlackScholesSolverMPI.hpp>
#include <sgpp/parallel/finance/algorithm/BlackScholesPATParabolicPDESolverSystemEuroAmerVectorizedMPI.hpp>

#include <sgpp/solver/ode/Euler.hpp>
#include <sgpp/solver/ode/CrankNicolson.hpp>
#include <sgpp/solver/sle/ConjugateGradients.hpp>

#include <sgpp/base/grid/Grid.hpp>
#include <sgpp/base/exception/application_exception.hpp>
#include <sgpp/base/operation/BaseOpFactory.hpp>
#include <sgpp/base/tools/StdNormalDistribution.hpp>
#include <cstdlib>
#include <cstring>
#include <sstream>
#include <cmath>
#include <fstream>
#include <iomanip>

#include <sgpp/globaldef.hpp>


namespace SGPP {
namespace parallel {

BlackScholesSolverMPI::BlackScholesSolverMPI(bool useLogTransform,
    bool usePAT) : SGPP::pde::ParabolicPDESolver() {

  this->bStochasticDataAlloc = false;
  this->bGridConstructed = false;
  this->myScreen = NULL;
  this->useCoarsen = false;
  this->coarsenThreshold = 0.0;
  this->adaptSolveMode = "none";
  this->refineMode = "classic";
  this->numCoarsenPoints = -1;
  this->useLogTransform = useLogTransform;
  this->usePAT = usePAT;

  if (this->usePAT == true) {
    this->useLogTransform = true;
  }

  this->refineMaxLevel = 0;
  this->nNeededIterations = 0;
  this->dNeededTime = 0.0;
  this->staInnerGridSize = 0;
  this->finInnerGridSize = 0;
  this->avgInnerGridSize = 0;
  this->current_time = 0.0;
  this->tBoundaryType = "freeBoundaries";
}

BlackScholesSolverMPI::~BlackScholesSolverMPI() {
  if (this->bStochasticDataAlloc) {
    delete this->mus;
    delete this->sigmas;
    delete this->rhos;
  }

  if (this->myScreen != NULL) {
    delete this->myScreen;
  }
}

void BlackScholesSolverMPI::getGridNormalDistribution(SGPP::base::DataVector&
    alpha, std::vector<double>& norm_mu, std::vector<double>& norm_sigma) {
  if (this->bGridConstructed) {
    double value;
    base::StdNormalDistribution myNormDistr;
    double* s_coords = new double[this->dim];

    for (size_t i = 0; i < this->myGrid->getStorage()->size(); i++) {
      std::string coords = this->myGridStorage->get(i)->getCoordsStringBB(*
                           (this->myBoundingBox));
      std::stringstream coordsStream(coords);

      for (size_t j = 0; j < this->dim; j++) {
        double tmp_load;
        coordsStream >> tmp_load;
        s_coords[j] = tmp_load;
      }

      value = 1.0;

      for (size_t j = 0; j < this->dim; j++) {
        if (this->useLogTransform == false) {
          value *= myNormDistr.getDensity(s_coords[j], norm_mu[j], norm_sigma[j]);
        } else {
          if (this->usePAT == true) {
            double inner_tmp = 0.0;

            for (size_t l = 0; l < dim; l++) {
              inner_tmp += this->eigvec_covar->get(j,
                                                   l) * (s_coords[l] - (this->current_time * this->mu_hat->get(l)));
            }

            value *= myNormDistr.getDensity(exp(inner_tmp), norm_mu[j], norm_sigma[j]);
          } else {
            value *= myNormDistr.getDensity(exp(s_coords[j]), norm_mu[j], norm_sigma[j]);
          }
        }
      }

      alpha[i] = value;
    }

    delete[] s_coords;
  } else {
    throw new base::application_exception("BlackScholesSolverMPI::getGridNormalDistribution : The grid wasn't initialized before!");
  }
}

void BlackScholesSolverMPI::constructGrid(SGPP::base::BoundingBox& BoundingBox,
    int level) {
  this->dim = BoundingBox.getDimensions();
  this->levels = level;

  this->myGrid = new base::LinearBoundaryGrid(BoundingBox);

  std::unique_ptr<SGPP::base::GridGenerator> myGenerator = this->myGrid->createGridGenerator();
  myGenerator->regular(this->levels);
  delete myGenerator;

  this->myBoundingBox = this->myGrid->getBoundingBox();
  this->myGridStorage = this->myGrid->getStorage();

  //std::string serGrid;
  //myGrid->serialize(serGrid);
  //std::cout << serGrid << std::endl;

  this->bGridConstructed = true;
}

void BlackScholesSolverMPI::setStochasticData(SGPP::base::DataVector& mus,
    SGPP::base::DataVector& sigmas, SGPP::base::DataMatrix& rhos, double r) {
  this->mus = new SGPP::base::DataVector(mus);
  this->sigmas = new SGPP::base::DataVector(sigmas);
  this->rhos = new SGPP::base::DataMatrix(rhos);
  this->r = r;

  // calculate eigenvalues, eigenvectors and mu_hat from stochastic data for PAT
  size_t mydim = this->mus->getSize();
  this->eigval_covar = new SGPP::base::DataVector(mydim);
  this->eigvec_covar = new SGPP::base::DataMatrix(mydim, mydim);
  this->mu_hat = new SGPP::base::DataVector(mydim);

  // 1d test case
  if (mydim == 1) {
    this->eigval_covar->set(0, this->sigmas->get(0)*this->sigmas->get(0));
    this->eigvec_covar->set(0, 0, 1.0);
  }

  // 2d test case
  if (mydim == 2) {
    // correlation -0.5, sigma_1 0.3, sigma_2 0.4
    this->eigval_covar->set(0, 0.0555377800527509792);
    this->eigval_covar->set(1, 0.194462219947249021);

    this->eigvec_covar->set(0, 0, -0.867142152569025494);
    this->eigvec_covar->set(0, 1, 0.498060726456078796);
    this->eigvec_covar->set(1, 0, -0.498060726456078796);
    this->eigvec_covar->set(1, 1, -0.867142152569025494);

    // correlation 0.1, sigma 0.4
    //    this->eigval_covar->set(0, 0.176);
    //    this->eigval_covar->set(1, 0.144);
    //
    //    this->eigvec_covar->set(0, 0, 0.707106781186547351);
    //    this->eigvec_covar->set(0, 1, -0.707106781186547573);
    //    this->eigvec_covar->set(1, 0, 0.707106781186547573);
    //    this->eigvec_covar->set(1, 1, 0.707106781186547351);

    // correlation 0.25, sigma 0.4
    //    this->eigval_covar->set(0, 0.20);
    //    this->eigval_covar->set(1, 0.12);
    //
    //    this->eigvec_covar->set(0, 0, 0.707106781186547573);
    //    this->eigvec_covar->set(0, 1, -0.707106781186547462);
    //    this->eigvec_covar->set(1, 0, 0.707106781186547462);
    //    this->eigvec_covar->set(1, 1, 0.707106781186547573);

    // correlation 0.5, sigma 0.4
    //    this->eigval_covar->set(0, 0.24);
    //    this->eigval_covar->set(1, 0.08);
    //
    //    this->eigvec_covar->set(0, 0, 0.707106781186547462);
    //    this->eigvec_covar->set(0, 1, -0.707106781186547462);
    //    this->eigvec_covar->set(1, 0, 0.707106781186547462);
    //    this->eigvec_covar->set(1, 1, 0.707106781186547462);

    // correlation -0.5, sigma 0.4
    //    this->eigval_covar->set(0, 0.24);
    //    this->eigval_covar->set(1, 0.08);
    //
    //    this->eigvec_covar->set(0, 0, 0.707106781186547462);
    //    this->eigvec_covar->set(0, 1, 0.707106781186547462);
    //    this->eigvec_covar->set(1, 0, -0.707106781186547462);
    //    this->eigvec_covar->set(1, 1, 0.707106781186547462);

    // correlation 0.0, sigma 0.4
    //    this->eigval_covar->set(0, 0.16);
    //    this->eigval_covar->set(1, 0.16);
    //
    //    this->eigvec_covar->set(0, 0, 1);
    //    this->eigvec_covar->set(0, 1, 0);
    //    this->eigvec_covar->set(1, 0, 0);
    //    this->eigvec_covar->set(1, 1, 1);
  }

  // 3d test case
  if (mydim == 3) {
    this->eigval_covar->set(0, 0.0161152062670340546);
    this->eigval_covar->set(1, 0.109759129882028184);
    this->eigval_covar->set(2, 0.164125663850937770);


    this->eigvec_covar->set(0, 0, -0.869836297894464927);
    this->eigvec_covar->set(0, 1, -0.472520156399458380);
    this->eigvec_covar->set(0, 2, -0.141808027493095845);

    this->eigvec_covar->set(1, 0, -0.493287590771869233);
    this->eigvec_covar->set(1, 1, 0.837245861704536964);
    this->eigvec_covar->set(1, 2, 0.235980337844304722);

    this->eigvec_covar->set(2, 0, -0.00722271802968979613);
    this->eigvec_covar->set(2, 1, -0.275216403680555832);
    this->eigvec_covar->set(2, 2, 0.961355170313971441);

    // american put test case, steffi
    //    this->eigval_covar->set(0, 0.0941905506508006474);
    //    this->eigval_covar->set(1, 0.211230588426859212);
    //    this->eigval_covar->set(2, 0.174578860922339874);
    //
    //
    //    this->eigvec_covar->set(0, 0, 0.592851303426220277);
    //    this->eigvec_covar->set(0, 1, -0.610834162257256286);
    //    this->eigvec_covar->set(0, 2, 0.524794205613311915);
    //
    //    this->eigvec_covar->set(1, 0, -0.658103087561460964);
    //    this->eigvec_covar->set(1, 1, -0.743068675017919733);
    //    this->eigvec_covar->set(1, 2, -0.121446574052896814);
    //
    //    this->eigvec_covar->set(2, 0, -0.464141851342780687);
    //    this->eigvec_covar->set(2, 1, 0.273368927324575406);
    //    this->eigvec_covar->set(2, 2, 0.842521080689087154);
  }

  // 4d test case
  if (mydim == 4) {
    this->eigval_covar->set(0, 0.203896808612126890);
    this->eigval_covar->set(1, 0.143228838600683389);
    this->eigval_covar->set(2, 0.0369289052706551282);
    this->eigval_covar->set(3, 0.884454475165345477e-1);


    this->eigvec_covar->set(0, 0, -0.758784156527507303);
    this->eigvec_covar->set(0, 1, -0.441609644035335480);
    this->eigvec_covar->set(0, 2, -0.329723883556484409);
    this->eigvec_covar->set(0, 3, 0.347145051398191185);

    this->eigvec_covar->set(1, 0, 0.0381555338704252095);
    this->eigvec_covar->set(1, 1, -0.502978780370165746e-1);
    this->eigvec_covar->set(1, 2, 0.713853106787664671);
    this->eigvec_covar->set(1, 3, 0.697443919343796126);

    this->eigvec_covar->set(2, 0, 0.327315044904916141);
    this->eigvec_covar->set(2, 1, 0.376967846478724333);
    this->eigvec_covar->set(2, 2, -0.600862700798999505);
    this->eigvec_covar->set(2, 3, 0.624278879098610684);

    this->eigvec_covar->set(3, 0, -0.561832377508447611);
    this->eigvec_covar->set(3, 1, 0.812616938342507478);
    this->eigvec_covar->set(3, 2, 0.143735581296216858);
    this->eigvec_covar->set(3, 3, -0.577769311359955864e-1);
  }

  // 5d test case
  if (mydim == 5) {
    this->eigval_covar->set(0, 0.248090694157781677);
    this->eigval_covar->set(1, 0.181003240820949346);
    this->eigval_covar->set(2, 0.0132179416451147155);
    this->eigval_covar->set(3, 0.0939549786605669707);
    this->eigval_covar->set(4, 0.0587331447155870698);


    this->eigvec_covar->set(0, 0, 0.263790550378285305);
    this->eigvec_covar->set(0, 1, 0.859923642135395405);
    this->eigvec_covar->set(0, 2, 0.334323851170300113);
    this->eigvec_covar->set(0, 3, -0.280069559202008711);
    this->eigvec_covar->set(0, 4, 0.0271012873267933163);

    this->eigvec_covar->set(1, 0, -0.460592496016481306e-1);
    this->eigvec_covar->set(1, 1, 0.127520536615915014e-1);
    this->eigvec_covar->set(1, 2, -0.377657780213884409);
    this->eigvec_covar->set(1, 3, -0.528472957191017390);
    this->eigvec_covar->set(1, 4, -0.758819389061222926);

    this->eigvec_covar->set(2, 0, -0.164069200805689486);
    this->eigvec_covar->set(2, 1, -0.378822519988808670);
    this->eigvec_covar->set(2, 2, 0.470800766942319815);
    this->eigvec_covar->set(2, 3, -0.728869775137000020);
    this->eigvec_covar->set(2, 4, 0.276893994941337263);

    this->eigvec_covar->set(3, 0, 0.749451659161123107);
    this->eigvec_covar->set(3, 1, -0.143542996441549164);
    this->eigvec_covar->set(3, 2, -0.467549699477188219);
    this->eigvec_covar->set(3, 3, -0.257715994147244276);
    this->eigvec_covar->set(3, 4, 0.364276493384776856);

    this->eigvec_covar->set(4, 0, -0.582834967194725273);
    this->eigvec_covar->set(4, 1, 0.310254123817747751);
    this->eigvec_covar->set(4, 2, -0.552581288090625233);
    this->eigvec_covar->set(4, 3, -0.211207700566551748);
    this->eigvec_covar->set(4, 4, 0.462699694124277694);
  }

  // 6d test case
  if (mydim == 6) {
    this->eigval_covar->set(0, 0.295078790670591784);
    this->eigval_covar->set(1, 0.208971923682093358);
    this->eigval_covar->set(2, 0.0128769346283334751);
    this->eigval_covar->set(3, 0.0498094136802732304);
    this->eigval_covar->set(4, 0.0774284377364329035);
    this->eigval_covar->set(5, 0.110834499602274789);


    this->eigvec_covar->set(0, 0, 0.0901970509435878476);
    this->eigvec_covar->set(0, 1, 0.825424523396568799);
    this->eigvec_covar->set(0, 2, 0.355559506316694030);
    this->eigvec_covar->set(0, 3, 0.198721715730790532);
    this->eigvec_covar->set(0, 4, 0.363046763193636213);
    this->eigvec_covar->set(0, 5, 0.113238743705900313);

    this->eigvec_covar->set(1, 0, 0.0365722370191640458);
    this->eigvec_covar->set(1, 1, -0.0144109917849721478);
    this->eigvec_covar->set(1, 2, -0.385578116247248359);
    this->eigvec_covar->set(1, 3, -0.505921195549970393);
    this->eigvec_covar->set(1, 4, 0.742524001100560604);
    this->eigvec_covar->set(1, 5, -0.206121718286344846);

    this->eigvec_covar->set(2, 0, 0.0739556327669701891);
    this->eigvec_covar->set(2, 1, -0.312594988132788687);
    this->eigvec_covar->set(2, 2, 0.484247830906325782);
    this->eigvec_covar->set(2, 3, 0.317369431761546084);
    this->eigvec_covar->set(2, 4, 0.263216340455829090);
    this->eigvec_covar->set(2, 5, -0.701650039506314549);

    this->eigvec_covar->set(3, 0, -.593995441134266944);
    this->eigvec_covar->set(3, 1, 0.333166368988965789);
    this->eigvec_covar->set(3, 2, -0.468239680986198614);
    this->eigvec_covar->set(3, 3, 0.270759276883645317);
    this->eigvec_covar->set(3, 4, -0.153188802604424607);
    this->eigvec_covar->set(3, 5, -0.469194834130158833);

    this->eigvec_covar->set(4, 0, 0.575964981307517609);
    this->eigvec_covar->set(4, 1, -0.0554637711185953924);
    this->eigvec_covar->set(4, 2, -0.516994216391904682);
    this->eigvec_covar->set(4, 3, 0.613357677157853276);
    this->eigvec_covar->set(4, 4, 0.135825834003240276);
    this->eigvec_covar->set(4, 5, 0.0569981868738118971);

    this->eigvec_covar->set(5, 0, 0.548170853002664993);
    this->eigvec_covar->set(5, 1, 0.326611380012385655);
    this->eigvec_covar->set(5, 2, -0.0622862213401522563);
    this->eigvec_covar->set(5, 3, -0.392825458984802756);
    this->eigvec_covar->set(5, 4, -0.453494043197070207);
    this->eigvec_covar->set(5, 5, -0.478524252837880804);
  }

  // 7d test case
  if (mydim == 7) {
    this->eigval_covar->set(0, 0.297277846593271722);
    this->eigval_covar->set(1, 0.208971930701126907);
    this->eigval_covar->set(2, 0.127627200649076539);
    this->eigval_covar->set(3, 0.0128255347503872102);
    this->eigval_covar->set(4, 0.0306260409692614491);
    this->eigval_covar->set(5, 0.0811043629942213989);
    this->eigval_covar->set(6, 0.0590670833426543418);


    this->eigvec_covar->set(0, 0, 0.0843288248079574199);
    this->eigvec_covar->set(0, 1, 0.825397719934479079);
    this->eigvec_covar->set(0, 2, -0.168709288751078074);
    this->eigvec_covar->set(0, 3, -0.363307928519013612);
    this->eigvec_covar->set(0, 4, -0.202137777735263119);
    this->eigvec_covar->set(0, 5, -0.331775600546961924);
    this->eigvec_covar->set(0, 6, -0.0147325453249897778);

    this->eigvec_covar->set(1, 0, 0.0388837766210140823);
    this->eigvec_covar->set(1, 1, -0.0144026678855673478);
    this->eigvec_covar->set(1, 2, 0.138436093520141523);
    this->eigvec_covar->set(1, 3, 0.391084279235813126);
    this->eigvec_covar->set(1, 4, 0.196762223947347770);
    this->eigvec_covar->set(1, 5, -0.670298189519680343);
    this->eigvec_covar->set(1, 6, 0.581510201090990941);

    this->eigvec_covar->set(2, 0, 0.0816282868053433697);
    this->eigvec_covar->set(2, 1, -0.312561527401625017);
    this->eigvec_covar->set(2, 2, 0.544757290477682066);
    this->eigvec_covar->set(2, 3, -0.481767290108688295);
    this->eigvec_covar->set(2, 4, -0.0150731089385307045);
    this->eigvec_covar->set(2, 5, -0.480781834905232630);
    this->eigvec_covar->set(2, 6, -0.367972166945297441);

    this->eigvec_covar->set(3, 0, -0.583758726401433226);
    this->eigvec_covar->set(3, 1, 0.333233472084286164);
    this->eigvec_covar->set(3, 2, 0.441249303136342808);
    this->eigvec_covar->set(3, 3, 0.468515871811986007);
    this->eigvec_covar->set(3, 4, -0.0558306699252684793);
    this->eigvec_covar->set(3, 5, -0.00669728666705872531);
    this->eigvec_covar->set(3, 6, -0.361678331424167054);

    this->eigvec_covar->set(4, 0, 0.567646273735049278);
    this->eigvec_covar->set(4, 1, -0.0555178986758357862);
    this->eigvec_covar->set(4, 2, -0.156878448552707023);
    this->eigvec_covar->set(4, 3, 0.504829137586297549);
    this->eigvec_covar->set(4, 4, -0.401843814498330509);
    this->eigvec_covar->set(4, 5, -0.202551001291546379);
    this->eigvec_covar->set(4, 6, -0.439006612353535552);

    this->eigvec_covar->set(5, 0, 0.558697722750611581);
    this->eigvec_covar->set(5, 1, 0.326633776766320549);
    this->eigvec_covar->set(5, 2, 0.476573019522854413);
    this->eigvec_covar->set(5, 3, 0.0821587688846129316);
    this->eigvec_covar->set(5, 4, 0.494200581047168797);
    this->eigvec_covar->set(5, 5, 0.320995293884535349);
    this->eigvec_covar->set(5, 6, 0.00480949306111711893);

    this->eigvec_covar->set(6, 0, 0.0978459803650174598);
    this->eigvec_covar->set(6, 1, 0.000224850494160951678);
    this->eigvec_covar->set(6, 2, 0.457360339661911197);
    this->eigvec_covar->set(6, 3, -0.0430240652502010823);
    this->eigvec_covar->set(6, 4, -0.715100005708458997);
    this->eigvec_covar->set(6, 5, 0.255664574715722459);
    this->eigvec_covar->set(6, 6, 0.450182374325719670);
  }

  for (size_t i = 0; i < mydim; i++) {
    double tmp = 0.0;

    for (size_t j = 0; j < mydim; j++) {
      tmp += ((this->mus->get(j) - (0.5 * this->sigmas->get(j) * this->sigmas->get(
                                      j))) * this->eigvec_covar->get(j, i));
    }

    this->mu_hat->set(i, tmp);
  }

  bStochasticDataAlloc = true;
}

void BlackScholesSolverMPI::solveExplicitEuler(size_t numTimesteps,
    double timestepsize, size_t maxCGIterations, double epsilonCG,
    SGPP::base::DataVector& alpha, bool verbose, bool generateAnimation,
    size_t numEvalsAnimation) {
  throw new base::application_exception("BlackScholesSolverMPI::solveExplicitEuler : Explicit Euler is not supported in this setting!");
}

void BlackScholesSolverMPI::solveImplicitEuler(size_t numTimesteps,
    double timestepsize, size_t maxCGIterations, double epsilonCG,
    SGPP::base::DataVector& alpha, bool verbose, bool generateAnimation,
    size_t numEvalsAnimation) {
  if (this->bGridConstructed && this->bStochasticDataAlloc) {
    solver::Euler* myEuler = new solver::Euler("ImEul", numTimesteps, timestepsize,
        generateAnimation, numEvalsAnimation, myScreen);
    solver::SLESolver* myCG;
    solver::OperationParabolicPDESolverSystem* myBSSystem = NULL;

    if (this->usePAT == false) {
      myCG = new parallel::BiCGStabMPI(maxCGIterations, epsilonCG);
      myBSSystem = new
      parallel::BlackScholesParabolicPDESolverSystemEuroAmerParallelMPI(*this->myGrid,
          alpha, *this->mus, *this->sigmas, *this->rhos, this->r, timestepsize, "ImEul",
          this->dStrike, this->payoffType, this->useLogTransform, this->useCoarsen,
          this->coarsenThreshold, this->adaptSolveMode, this->numCoarsenPoints,
          this->refineThreshold, this->refineMode,
          static_cast<SGPP::base::GridIndex::level_type>(this->refineMaxLevel));
    } else {
      // read env variable, which solver type should be selected
      char* alg_selector = getenv("SGPP_PDE_SOLVER_ALG");

      if (alg_selector != NULL) {
        if (! strcmp(alg_selector, "X86SIMD")) {
          myCG = new solver::ConjugateGradients(maxCGIterations, epsilonCG);
          myBSSystem = new BlackScholesPATParabolicPDESolverSystemEuroAmerVectorizedMPI(
            *this->myGrid, alpha, *this->eigval_covar, *this->eigvec_covar, *this->mu_hat,
            timestepsize, "ImEul", this->dStrike, this->payoffType, this->r,
            this->useCoarsen, this->coarsenThreshold, this->adaptSolveMode,
            this->numCoarsenPoints, this->refineThreshold, this->refineMode,
            this->refineMaxLevel);
        } else if (! strcmp(alg_selector, "OCL")) {
          myCG = new solver::ConjugateGradients(maxCGIterations, epsilonCG);
          myBSSystem = new BlackScholesPATParabolicPDESolverSystemEuroAmerVectorizedMPI(
            *this->myGrid, alpha, *this->eigval_covar, *this->eigvec_covar, *this->mu_hat,
            timestepsize, "ImEul", this->dStrike, this->payoffType, this->r,
            this->useCoarsen, this->coarsenThreshold, this->adaptSolveMode,
            this->numCoarsenPoints, this->refineThreshold, this->refineMode,
            this->refineMaxLevel);
        } else {
          throw new base::application_exception("BlackScholesSolverMPI::solveImplicitEuler : You have selected an unsupport vectorization method!");
        }
      } else {
        myCG = new parallel::ConjugateGradientsMPI(maxCGIterations, epsilonCG);
        myBSSystem = new
        parallel::BlackScholesPATParabolicPDESolverSystemEuroAmerParallelMPI(
          *this->myGrid, alpha, *this->eigval_covar, *this->eigvec_covar, *this->mu_hat,
          timestepsize, "ImEul", this->dStrike, this->payoffType, this->r,
          this->useCoarsen, this->coarsenThreshold, this->adaptSolveMode,
          this->numCoarsenPoints, this->refineThreshold, this->refineMode,
          static_cast<SGPP::base::GridIndex::level_type>(this->refineMaxLevel));
      }
    }

    base::SGppStopwatch* myStopwatch = new base::SGppStopwatch();
    this->staInnerGridSize = getNumberInnerGridPoints();

    if (myGlobalMPIComm->getMyRank() == 0) {
      std::cout << "Using Implicit Euler to solve " << numTimesteps << " timesteps:"
                << std::endl;
    }

    MPI_Barrier(MPI_COMM_WORLD);
    myStopwatch->start();
    myEuler->solve(*myCG, *myBSSystem, true, verbose);
    MPI_Barrier(MPI_COMM_WORLD);
    this->dNeededTime = myStopwatch->stop();

    if (myGlobalMPIComm->getMyRank() == 0) {
      std::cout << "Using Implicit Euler to solve " << numTimesteps << " timesteps:"
                << std::endl;
      std::cout << std::endl << "Final Grid size: " << getNumberGridPoints() <<
                std::endl;
      std::cout << "Final Grid size (inner): " << getNumberInnerGridPoints() <<
                std::endl << std::endl << std::endl;
      std::cout << "Average Grid size: " << static_cast<double>
                (myBSSystem->getSumGridPointsComplete()) / static_cast<double>
                (numTimesteps) << std::endl;
      std::cout << "Average Grid size (Inner): " << static_cast<double>
                (myBSSystem->getSumGridPointsInner()) / static_cast<double>
                (numTimesteps) << std::endl << std::endl << std::endl;

      if (this->myScreen != NULL) {
        std::cout << "Time to solve: " << this->dNeededTime << " seconds" << std::endl;
        this->myScreen->writeEmptyLines(2);
      }
    }

    this->finInnerGridSize = getNumberInnerGridPoints();
    this->avgInnerGridSize = static_cast<size_t>((static_cast<double>
                             (myBSSystem->getSumGridPointsInner()) / static_cast<double>
                             (numTimesteps)) + 0.5);
    this->nNeededIterations = myEuler->getNumberIterations();

    delete myBSSystem;
    delete myCG;
    delete myEuler;
    delete myStopwatch;

    this->current_time += (static_cast<double>(numTimesteps) * timestepsize);
  } else {
    throw new base::application_exception("BlackScholesSolver::solveImplicitEuler : A grid wasn't constructed before or stochastic parameters weren't set!");
  }
}

void BlackScholesSolverMPI::solveCrankNicolson(size_t numTimesteps,
    double timestepsize, size_t maxCGIterations, double epsilonCG,
    SGPP::base::DataVector& alpha, size_t NumImEul) {
  if (this->bGridConstructed && this->bStochasticDataAlloc) {
    solver::SLESolver* myCG;
    solver::OperationParabolicPDESolverSystem* myBSSystem = NULL;

    if (this->usePAT == false) {
      myCG = new parallel::BiCGStabMPI(maxCGIterations, epsilonCG);
      myBSSystem = new
      parallel::BlackScholesParabolicPDESolverSystemEuroAmerParallelMPI(*this->myGrid,
          alpha, *this->mus, *this->sigmas, *this->rhos, this->r, timestepsize, "ImEul",
          this->dStrike, this->payoffType, this->useLogTransform, this->useCoarsen,
          this->coarsenThreshold, this->adaptSolveMode, this->numCoarsenPoints,
          this->refineThreshold, this->refineMode,
          static_cast<SGPP::base::GridIndex::level_type>(this->refineMaxLevel));
    } else {
      // read env variable, which solver type should be selected
      char* alg_selector = getenv("SGPP_PDE_SOLVER_ALG");

      if (alg_selector != NULL) {
        if (! strcmp(alg_selector, "X86SIMD")) {
          myCG = new solver::ConjugateGradients(maxCGIterations, epsilonCG);
          myBSSystem = new BlackScholesPATParabolicPDESolverSystemEuroAmerVectorizedMPI(
            *this->myGrid, alpha, *this->eigval_covar, *this->eigvec_covar, *this->mu_hat,
            timestepsize, "ImEul", this->dStrike, this->payoffType, this->r,
            this->useCoarsen, this->coarsenThreshold, this->adaptSolveMode,
            this->numCoarsenPoints, this->refineThreshold, this->refineMode,
            this->refineMaxLevel);
        } else if (! strcmp(alg_selector, "OCL")) {
          myCG = new solver::ConjugateGradients(maxCGIterations, epsilonCG);
          myBSSystem = new BlackScholesPATParabolicPDESolverSystemEuroAmerVectorizedMPI(
            *this->myGrid, alpha, *this->eigval_covar, *this->eigvec_covar, *this->mu_hat,
            timestepsize, "ImEul", this->dStrike, this->payoffType, this->r,
            this->useCoarsen, this->coarsenThreshold, this->adaptSolveMode,
            this->numCoarsenPoints, this->refineThreshold, this->refineMode,
            this->refineMaxLevel);
        } else {
          throw new base::application_exception("BlackScholesSolverMPI::solveCrankNicolson : You have selected an unsupport vectorization method!");
        }
      } else {
        myCG = new parallel::ConjugateGradientsMPI(maxCGIterations, epsilonCG);
        myBSSystem = new
        parallel::BlackScholesPATParabolicPDESolverSystemEuroAmerParallelMPI(
          *this->myGrid, alpha, *this->eigval_covar, *this->eigvec_covar, *this->mu_hat,
          timestepsize, "ImEul", this->dStrike, this->payoffType, this->r,
          this->useCoarsen, this->coarsenThreshold, this->adaptSolveMode,
          this->numCoarsenPoints, this->refineThreshold, this->refineMode,
          static_cast<SGPP::base::GridIndex::level_type>(this->refineMaxLevel));
      }
    }

    base::SGppStopwatch* myStopwatch = new base::SGppStopwatch();
    this->staInnerGridSize = getNumberInnerGridPoints();

    size_t numCNSteps;
    size_t numIESteps;

    numCNSteps = numTimesteps;

    if (numTimesteps > NumImEul) {
      numCNSteps = numTimesteps - NumImEul;
    }

    numIESteps = NumImEul;

    solver::Euler* myEuler = new solver::Euler("ImEul", numIESteps, timestepsize,
        false, 0, this->myScreen);
    solver::CrankNicolson* myCN = new solver::CrankNicolson(numCNSteps,
        timestepsize, this->myScreen);

    MPI_Barrier(MPI_COMM_WORLD);
    myStopwatch->start();

    if (numIESteps > 0) {
      if (myGlobalMPIComm->getMyRank() == 0) {
        std::cout << "Using Implicit Euler to solve " << numIESteps << " timesteps:" <<
                  std::endl;
      }

      myBSSystem->setODESolver("ImEul");
      myEuler->solve(*myCG, *myBSSystem, false, false);
    }

    if (myGlobalMPIComm->getMyRank() == 0) {
      std::cout << "Using Crank Nicolson to solve " << numCNSteps << " timesteps:" <<
                std::endl << std::endl << std::endl << std::endl;
    }

    myBSSystem->setODESolver("CrNic");
    myCN->solve(*myCG, *myBSSystem, true, false);
    MPI_Barrier(MPI_COMM_WORLD);
    this->dNeededTime = myStopwatch->stop();

    if (myGlobalMPIComm->getMyRank() == 0) {
      std::cout << std::endl << "Final Grid size: " << getNumberGridPoints() <<
                std::endl;
      std::cout << "Final Grid size (inner): " << getNumberInnerGridPoints() <<
                std::endl << std::endl << std::endl;
      std::cout << "Average Grid size: " << static_cast<double>
                (myBSSystem->getSumGridPointsComplete()) / static_cast<double>
                (numTimesteps) << std::endl;
      std::cout << "Average Grid size (Inner): " << static_cast<double>
                (myBSSystem->getSumGridPointsInner()) / static_cast<double>
                (numTimesteps) << std::endl << std::endl << std::endl;

      if (this->myScreen != NULL) {
        std::cout << "Time to solve: " << this->dNeededTime << " seconds" << std::endl;
        this->myScreen->writeEmptyLines(2);
      }
    }

    this->finInnerGridSize = getNumberInnerGridPoints();
    this->avgInnerGridSize = static_cast<size_t>((static_cast<double>
                             (myBSSystem->getSumGridPointsInner()) / static_cast<double>
                             (numTimesteps)) + 0.5);
    this->nNeededIterations = myEuler->getNumberIterations() +
                              myCN->getNumberIterations();

    delete myBSSystem;
    delete myCG;
    delete myCN;
    delete myEuler;
    delete myStopwatch;

    this->current_time += (static_cast<double>(numTimesteps) * timestepsize);
  } else {
    throw new base::application_exception("BlackScholesSolver::solveCrankNicolson : A grid wasn't constructed before or stochastic parameters weren't set!");
  }
}


void BlackScholesSolverMPI::solveAdamsBashforth(size_t numTimesteps,
    double timestepsize, size_t maxCGIterations, double epsilonCG,
    SGPP::base::DataVector& alpha, bool verbose) {
  throw new base::application_exception("BlackScholesSolverMPI::solveAdamsBashforth : An unsupported ODE Solver type has been chosen!");
}

void BlackScholesSolverMPI::solveSCAC(size_t numTimesteps, double timestepsize,
                                      double epsilon, size_t maxCGIterations, double epsilonCG,
                                      SGPP::base::DataVector& alpha, bool verbose) {
  throw new base::application_exception("BlackScholesSolverMPI::solveSCAC : An unsupported ODE Solver type has been chosen!");
}

void BlackScholesSolverMPI::solveSCH(size_t numTimesteps, double timestepsize,
                                     double epsilon, size_t maxCGIterations, double epsilonCG,
                                     SGPP::base::DataVector& alpha, bool verbose) {
  throw new base::application_exception("BlackScholesSolverMPI::solveSCH : An unsupported ODE Solver type has been chosen!");
}

void BlackScholesSolverMPI::solveSCBDF(size_t numTimesteps, double timestepsize,
                                       double epsilon, size_t maxCGIterations, double epsilonCG,
                                       SGPP::base::DataVector& alpha, bool verbose) {
  throw new base::application_exception("BlackScholesSolverMPI::solveSCBDF : An unsupported ODE Solver type has been chosen!");
}

void BlackScholesSolverMPI::solveSCEJ(size_t numTimesteps, double timestepsize,
                                      double epsilon, double myAlpha, size_t maxCGIterations, double epsilonCG,
                                      SGPP::base::DataVector& alpha, bool verbose) {
  throw new base::application_exception("BlackScholesSolverMPI::solveSCEJ : An unsupported ODE Solver type has been chosen!");
}

void BlackScholesSolverMPI::solveX(size_t numTimesteps, double timestepsize,
                                   size_t maxCGIterations, double epsilonCG, SGPP::base::DataVector& alpha,
                                   bool verbose, void* myODESolverV, std::string Solver) {
  throw new base::application_exception("BlackScholesSolverMPI::solveX : An unsupported ODE Solver type has been chosen!");
}

void BlackScholesSolverMPI::setPayoffType(std::string payoffType) {
  this->payoffType = payoffType;
}

void BlackScholesSolverMPI::initGridWithPayoff(SGPP::base::DataVector& alpha,
    double strike, std::string payoffType) {
  this->dStrike = strike;
  this->payoffType = payoffType;

  if (payoffType == "std_euro_call" || payoffType == "std_euro_put"
      || payoffType == "std_amer_put") {
    this->tBoundaryType = "Dirichlet";
  }

  if (this->useLogTransform) {
    if (this->usePAT == true) {
      initPATTransformedGridWithPayoff(alpha, strike, payoffType);
    } else {
      initLogTransformedGridWithPayoff(alpha, strike, payoffType);
    }
  } else {
    initCartesianGridWithPayoff(alpha, strike, payoffType);
  }
}

double BlackScholesSolverMPI::get1DEuroCallPayoffValue(double assetValue,
    double strike) {
  if (assetValue <= strike) {
    return 0.0;
  } else {
    return assetValue - strike;
  }
}

std::vector<size_t> BlackScholesSolverMPI::getAlgorithmicDimensions() {
  return this->myGrid->getAlgorithmicDimensions();
}

void BlackScholesSolverMPI::setAlgorithmicDimensions(std::vector<size_t>
    newAlgoDims) {
  if (this->tBoundaryType == "freeBoundaries") {
    this->myGrid->setAlgorithmicDimensions(newAlgoDims);
  } else {
    throw new base::application_exception("BlackScholesSolver::setAlgorithmicDimensions : Set algorithmic dimensions is only supported when choosing option type all!");
  }
}

void BlackScholesSolverMPI::initScreen() {
  this->myScreen = new base::ScreenOutput();
  this->myScreen->writeTitle("SGpp - Black Scholes Solver, 2.0.0",
                             "The SG++ Project (C) 2009-2010, by Alexander Heinecke");
  this->myScreen->writeStartSolve("Multidimensional Black Scholes Solver");
}

void BlackScholesSolverMPI::setEnableCoarseningData(std::string adaptSolveMode,
    std::string refineMode, int refineMaxLevel, int numCoarsenPoints,
    double coarsenThreshold, double refineThreshold) {
  this->useCoarsen = true;
  this->coarsenThreshold = coarsenThreshold;
  this->refineThreshold = refineThreshold;
  this->refineMaxLevel = refineMaxLevel;
  this->adaptSolveMode = adaptSolveMode;
  this->refineMode = refineMode;
  this->numCoarsenPoints = numCoarsenPoints;
}

void BlackScholesSolverMPI::initCartesianGridWithPayoff(
  SGPP::base::DataVector& alpha, double strike, std::string payoffType) {
  double tmp;

  if (this->bGridConstructed) {
    for (size_t i = 0; i < this->myGrid->getStorage()->size(); i++) {
      std::string coords = this->myGridStorage->get(i)->getCoordsStringBB(
                             *this->myBoundingBox);
      std::stringstream coordsStream(coords);
      double* dblFuncValues = new double[dim];

      for (size_t j = 0; j < this->dim; j++) {
        coordsStream >> tmp;

        dblFuncValues[j] = tmp;
      }

      if (payoffType == "std_euro_call") {
        tmp = 0.0;

        for (size_t j = 0; j < dim; j++) {
          tmp += dblFuncValues[j];
        }

        alpha[i] = std::max<double>(((tmp / static_cast<double>(dim)) - strike), 0.0);
      } else if (payoffType == "std_euro_put" || payoffType == "std_amer_put") {
        tmp = 0.0;

        for (size_t j = 0; j < dim; j++) {
          tmp += dblFuncValues[j];
        }

        alpha[i] = std::max<double>(strike - ((tmp / static_cast<double>(dim))), 0.0);
      } else {
        throw new base::application_exception("BlackScholesSolver::initCartesianGridWithPayoff : An unknown payoff-type was specified!");
      }

      delete[] dblFuncValues;
    }

    base::OperationHierarchisation* myHierarchisation =
      SGPP::op_factory::createOperationHierarchisation(*this->myGrid);
    myHierarchisation->doHierarchisation(alpha);
    delete myHierarchisation;
  } else {
    throw new base::application_exception("BlackScholesSolver::initCartesianGridWithPayoff : A grid wasn't constructed before!");
  }
}

void BlackScholesSolverMPI::initLogTransformedGridWithPayoff(
  base::DataVector& alpha, double strike, std::string payoffType) {
  double tmp;

  if (this->bGridConstructed) {
    for (size_t i = 0; i < this->myGrid->getStorage()->size(); i++) {
      std::string coords = this->myGridStorage->get(i)->getCoordsStringBB(
                             *this->myBoundingBox);
      std::stringstream coordsStream(coords);
      double* dblFuncValues = new double[dim];

      for (size_t j = 0; j < this->dim; j++) {
        coordsStream >> tmp;

        dblFuncValues[j] = tmp;
      }

      if (payoffType == "std_euro_call") {
        tmp = 0.0;

        for (size_t j = 0; j < dim; j++) {
          tmp += exp(dblFuncValues[j]);
        }

        alpha[i] = std::max<double>(((tmp / static_cast<double>(dim)) - strike), 0.0);
      } else if (payoffType == "std_euro_put" || payoffType == "std_amer_put") {
        tmp = 0.0;

        for (size_t j = 0; j < dim; j++) {
          tmp += exp(dblFuncValues[j]);
        }

        alpha[i] = std::max<double>(strike - ((tmp / static_cast<double>(dim))), 0.0);
      } else {
        throw new base::application_exception("BlackScholesSolver::initLogTransformedGridWithPayoff : An unknown payoff-type was specified!");
      }

      delete[] dblFuncValues;
    }

    base::OperationHierarchisation* myHierarchisation =
      SGPP::op_factory::createOperationHierarchisation(*this->myGrid);
    myHierarchisation->doHierarchisation(alpha);
    delete myHierarchisation;
  } else {
    throw new base::application_exception("BlackScholesSolver::initLogTransformedGridWithPayoff : A grid wasn't constructed before!");
  }
}

void BlackScholesSolverMPI::initPATTransformedGridWithPayoff(
  base::DataVector& alpha, double strike, std::string payoffType) {
  double tmp;

  if (this->bGridConstructed) {
    for (size_t i = 0; i < this->myGrid->getStorage()->size(); i++) {
      std::string coords = this->myGridStorage->get(i)->getCoordsStringBB(
                             *this->myBoundingBox);
      std::stringstream coordsStream(coords);
      double* dblFuncValues = new double[dim];

      for (size_t j = 0; j < this->dim; j++) {
        coordsStream >> tmp;

        dblFuncValues[j] = tmp;
      }

      if (payoffType == "std_euro_call") {
        tmp = 0.0;

        for (size_t j = 0; j < dim; j++) {
          double inner_tmp = 0.0;

          for (size_t l = 0; l < dim; l++) {
            inner_tmp += this->eigvec_covar->get(j, l) * dblFuncValues[l];
          }

          tmp += exp(inner_tmp);
        }

        alpha[i] = std::max<double>(((tmp / static_cast<double>(dim)) - strike), 0.0);
      } else if (payoffType == "std_euro_put" || payoffType == "std_amer_put") {
        tmp = 0.0;

        for (size_t j = 0; j < dim; j++) {
          double inner_tmp = 0.0;

          for (size_t l = 0; l < dim; l++) {
            inner_tmp += this->eigvec_covar->get(j, l) * dblFuncValues[l];
          }

          tmp += exp(inner_tmp);
        }

        alpha[i] = std::max<double>(strike - ((tmp / static_cast<double>(dim))), 0.0);
      } else {
        throw new base::application_exception("BlackScholesSolverMPI::initPATTransformedGridWithPayoff : An unknown payoff-type was specified!");
      }

      delete[] dblFuncValues;
    }

    base::OperationHierarchisation* myHierarchisation =
      SGPP::op_factory::createOperationHierarchisation(*this->myGrid);
    myHierarchisation->doHierarchisation(alpha);
    delete myHierarchisation;
  } else {
    throw new base::application_exception("BlackScholesSolverMPI::initPATTransformedGridWithPayoff : A grid wasn't constructed before!");
  }
}

double BlackScholesSolverMPI::evalOption(std::vector<double>& eval_point,
    SGPP::base::DataVector& alpha) {
  std::vector<double> trans_eval = eval_point;

  // apply needed coordinate transformations
  if (this->useLogTransform) {
    if (this->usePAT) {
      for (size_t i = 0; i < eval_point.size(); i++) {
        double trans_point = 0.0;

        for (size_t j = 0; j < this->dim; j++) {
          trans_point += (this->eigvec_covar->get(j, i) * (log(eval_point[j])));
        }

        trans_point += (this->current_time * this->mu_hat->get(i));

        trans_eval[i] = trans_point;
      }
    } else {
      for (size_t i = 0; i < eval_point.size(); i++) {
        trans_eval[i] = log(trans_eval[i]);
      }
    }
  }

  SGPP::base::OperationEval* myEval = SGPP::op_factory::createOperationEval(
                                        *this->myGrid);
  double result = myEval->eval(alpha, trans_eval);
  delete myEval;

  // discounting, if PAT is used
  // discounting, if PAT is used
  if (this->usePAT == true && this->payoffType != "std_amer_put") {
    result *= exp(((-1.0) * (this->r * this->current_time)));
  }

  return result;
}

void BlackScholesSolverMPI::transformPoint(SGPP::base::DataVector& point) {
  SGPP::base::DataVector tmp_point(point);

  // apply needed coordinate transformations
  if (this->useLogTransform) {
    if (this->usePAT) {
      for (size_t i = 0; i < point.getSize(); i++) {
        double trans_point = 0.0;

        for (size_t j = 0; j < point.getSize(); j++) {
          trans_point += (this->eigvec_covar->get(j, i) * (log(point[j])));
        }

        trans_point += (this->current_time * this->mu_hat->get(i));

        tmp_point[i] = trans_point;
      }
    } else {
      for (size_t i = 0; i < point.getSize(); i++) {
        tmp_point[i] = log(point[i]);
      }
    }
  }

  point = tmp_point;
}

void BlackScholesSolverMPI::printSparseGridPAT(SGPP::base::DataVector& alpha,
    std::string tfilename, bool bSurplus) const {
  base::DataVector temp(alpha);
  double tmp = 0.0;
  size_t dim = myGrid->getStorage()->dim();
  std::ofstream fileout;

  // Do Dehierarchisation, is specified
  if (bSurplus == false) {
    base::OperationHierarchisation* myHier =
      SGPP::op_factory::createOperationHierarchisation(*myGrid);
    myHier->doDehierarchisation(temp);
    delete myHier;
  }

  // Open filehandle
  fileout.open(tfilename.c_str());

  for (size_t i = 0; i < myGrid->getStorage()->size(); i++) {
    std::string coords =  myGrid->getStorage()->get(i)->getCoordsStringBB(
                            *myGrid->getBoundingBox());
    std::stringstream coordsStream(coords);

    double* dblFuncValues = new double[dim];

    for (size_t j = 0; j < dim; j++) {
      coordsStream >> tmp;
      dblFuncValues[j] = tmp;
    }

    for (size_t l = 0; l < dim; l++) {
      double trans_point = 0.0;

      for (size_t j = 0; j < dim; j++) {
        trans_point += this->eigvec_covar->get(l,
                                               j) * (dblFuncValues[j] - (this->current_time * this->mu_hat->get(j)));
      }

      fileout << exp(trans_point) << " ";
    }

    fileout << temp[i] << std::endl;

    delete[] dblFuncValues;
  }

  fileout.close();
}

void BlackScholesSolverMPI::resetSolveTime() {
  this->current_time = 0.0;
}


size_t BlackScholesSolverMPI::getNeededIterationsToSolve() {
  return this->nNeededIterations;
}

double BlackScholesSolverMPI::getNeededTimeToSolve() {
  return this->dNeededTime;
}

size_t BlackScholesSolverMPI::getStartInnerGridSize() {
  return this->staInnerGridSize;
}

size_t BlackScholesSolverMPI::getFinalInnerGridSize() {
  return this->finInnerGridSize;
}

size_t BlackScholesSolverMPI::getAverageInnerGridSize() {
  return this->avgInnerGridSize;
}

}

}
