// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/datadriven/algorithm/AlgorithmAdaBoostBase.hpp>
#include <sgpp/base/exception/operation_exception.hpp>
#include <sgpp/base/operation/BaseOpFactory.hpp>

#include <sgpp/globaldef.hpp>


namespace SGPP {
  namespace datadriven {
    AlgorithmAdaBoostBase::AlgorithmAdaBoostBase(SGPP::base::Grid& SparseGrid, size_t gridType, SGPP::base::HashGenerator::level_t gridLevel, SGPP::base::DataMatrix& trainData, SGPP::base::DataVector& trainDataClass, size_t NUM, float_t lambda, size_t IMAX, float_t eps, size_t IMAX_final, float_t eps_final, float_t firstLabel, float_t secondLabel, float_t threshold, float_t maxLambda, float_t minLambda, size_t searchNum, bool refine, size_t refineMode, size_t refineNum, size_t numberOfAda, float_t percentOfAda, size_t mode) {
      if (refine && (gridType != 1 && gridType != 2 && gridType != 3)) {
        throw new SGPP::base::operation_exception("AlgorithmAdaBoostBase : Only 1 or 2 or 3 are supported gridType(1 = Linear Grid, 2 = LinearL0Boundary Grid, 3 = ModLinear Grid)!");
      }

      if (refine && (percentOfAda >= 1.0 || percentOfAda <= 0.0)) {
        throw new SGPP::base::operation_exception("AlgorithmAdaBoostBase : Only number between 0 and 1 is the supported percent to Adaptive!");
      }

      if (refineMode != 1 && refineMode != 2) {
        throw new SGPP::base::operation_exception("AlgorithmAdaBoostBase : Only 1 or 2 are supported refine mode(1 : use grid point number, 2: use grid point percentage)!");
      }

      SGPP::base::GridStorage* gridStorage = SparseGrid.getStorage();
      this->grid = &SparseGrid;
      this->type = gridType;
      this->gridPoint = gridStorage->size();
      this->level = static_cast<SGPP::base::HashGenerator::level_t>(gridLevel);
      this->lamb = lambda;
      this->data = &trainData;
      this->classes = &trainDataClass;
      this->numData = trainData.getNrows();
      this->dim = gridStorage->dim();
      this->numBaseLearners = NUM;
      this->imax = IMAX;
      this->epsilon = eps;
      this->imax_final = IMAX_final;
      this->epsilon_final = eps_final;
      this->labelOne = firstLabel;
      this->labelTwo = secondLabel;
      this->threshold = threshold;
      this->lambLogMax = log(maxLambda);
      this->lambSteps = searchNum;

      if (searchNum == 1)
        this->lambStepsize = (log((float_t)maxLambda) - log((float_t)minLambda)) / 2;
      else
        this->lambStepsize = (log((float_t)maxLambda) - log((float_t)minLambda)) / ((float_t)searchNum - 1);

      this->actualBaseLearners = 0;
      this->refinement = refine;
      this->refineMode = refineMode;
      this->refineTimes = refineNum;
      this->numOfAda = numberOfAda;
      this->perOfAda = percentOfAda;
      this->maxGridPoint = new SGPP::base::DataVector(NUM);
      this->sumGridPoint = new SGPP::base::DataVector(NUM);
      this->boostMode = mode;
    }

    AlgorithmAdaBoostBase::~AlgorithmAdaBoostBase() {

    }

    void AlgorithmAdaBoostBase::doDiscreteAdaBoost(SGPP::base::DataVector& hypoWeight, SGPP::base::DataVector& weightError, SGPP::base::DataMatrix& weights, SGPP::base::DataMatrix& decision, SGPP::base::DataMatrix& testData, SGPP::base::DataMatrix& algorithmValueTrain, SGPP::base::DataMatrix& algorithmValueTest) {
      SGPP::base::DataVector weight(this->numData);
      weight.setAll(1.0 / float_t(this->numData));
      SGPP::base::OperationEval* opEval = SGPP::op_factory::createOperationEval(*this->grid);
      // to store certain train data point
      SGPP::base::DataVector p_train(this->dim);
      // to store certain train data point
      SGPP::base::DataVector p_test(this->dim);
      // create vector to store the hypothesis of the training data according to certain alpha vector(base learner)
      SGPP::base::DataVector newclasses(this->numData);
      // create vector to store the identity of comparing hypothesis of training data with class of
      SGPP::base::DataVector identity(this->numData);

      SGPP::base::DataVector tmpweight(this->numData);

      for (size_t count = 0; count < this->numBaseLearners; count++) {
        (this->actualBaseLearners)++;
        std::cout << std::endl;
        std::cout << "This is the " << this->actualBaseLearners << "th weak learner." << std::endl;
        std::cout << std::endl;

        // create coefficient vector
        SGPP::base::DataVector alpha_train(this->gridPoint);
        SGPP::base::DataVector alpha_learn(this->gridPoint);
        std::cout << "gridPoint: " << this->gridPoint << std::endl;

        if (this->maxGridPoint->get(count) < this->gridPoint)
          this->maxGridPoint->set(count, (float_t)this->gridPoint);

        if (!this->refinement) {
          if (count == 0)
            this->sumGridPoint->set(count, (float_t)gridPoint);
          else
            this->sumGridPoint->set(count, this->sumGridPoint->get(count - 1) + (float_t)gridPoint);
        }

        alpha_train.setAll(0.0);
        weights.setColumn(count, weight);

        bool final_step = false;

        if (this->refinement == 0)
          final_step = true;

        // calculate alpha
        alphaSolver(this->lamb, weight, alpha_train, final_step);

        if (this->refinement) {
          doRefinement(alpha_train, weight, count + 1);
          opEval = SGPP::op_factory::createOperationEval(*this->grid);
          alpha_learn.resizeZero(alpha_train.getSize());
        }

        //set the alpha for testing data(copy of alpha for training data)
        alpha_learn.copyFrom(alpha_train);

        #pragma omp parallel for schedule(static)

        for (size_t i = 0; i < this->numData; i++) {
          SGPP::base::DataVector p_train_private(this->dim);
          this->data->getRow(i, p_train_private);
          float_t value_train = opEval->eval(alpha_train, p_train_private);
          newclasses.set(i, hValue(value_train));
        }

        for (size_t i = 0; i < this->numData; i++) {
          if (newclasses.get(i) == this->classes->get(i)) {
            identity.set(i, 0.0);
            decision.set(i, count, 1.0);
          } else {
            identity.set(i, 1.0);
            decision.set(i, count, 0.0);
          }
        }

        // calculate the weight error
        weightError.set(count, weight.dotProduct(identity));

        // find the optimal lambda to minimize the weighted error
        if (this->lambSteps > 0 && count > 0) {
          float_t cur_lambda;
          float_t weighterror;
          float_t minWeightError = weightError.get(count);

          for (size_t it = 0; it < this->lambSteps; it++) {
            std::cout << std::endl;
            std::cout << "This is the " << it + 1 << "th search of " << this->actualBaseLearners << "th weak learner." << std::endl;
            std::cout << std::endl;
            alpha_train.setAll(0.0);
            cur_lambda = exp(this->lambLogMax - (float_t)it * this->lambStepsize);

            alphaSolver(cur_lambda, weight, alpha_train, true);

            for (size_t i = 0; i < this->numData; i++) {
              this->data->getRow(i, p_train);
              float_t value_seach = opEval->eval(alpha_train, p_train);
              newclasses.set(i, hValue(value_seach));
            }

            for (size_t i = 0; i < this->numData; i++) {
              if (newclasses.get(i) == this->classes->get(i)) {
                identity.set(i, 0.0);
              } else {
                identity.set(i, 1.0);
              }
            }

            // compare the weight error we need the minimum weight error
            weighterror = weight.dotProduct(identity);

            if (weighterror < minWeightError) {
              minWeightError = weighterror;
              weightError.set(count, weighterror);

              //reset the alpha for testing data(copy of alpha for training data)
              for (size_t index = 0; index < alpha_train.getSize(); index++)
                alpha_learn.set(index, alpha_train.get(index));

              for (size_t i = 0; i < this->numData; i++) {
                if (newclasses.get(i) == this->classes->get(i)) {
                  decision.set(i, count, 1.0);
                } else {
                  decision.set(i, count, 0.0);
                }
              }
            }
          }
        }

        // to judge the classification whether match the simple condition of Adaboost
        if (weightError.get(count) >= 0.50) {
          std::cout << std::endl << "The training error rate exceeds 0.5 after " << count + 1 << " iterations" << std::endl;
          (this->actualBaseLearners)--;
          break;
        }

        // calculate the weight of this weak classif
        float_t hypoweight;

        if (weightError.get(count) == 0) {
          hypoweight = log(1e+10);
          hypoWeight.set(count, hypoweight);
        } else {
          hypoweight = 0.5 * (log((1.0 - weightError.get(count)) / weightError.get(count)));
          hypoWeight.set(count, hypoweight);
        }

        // calculate the algorithm value of the testing data and training data
        // for training data
        #pragma omp parallel for schedule(static)

        for (size_t i = 0; i < numData; i++) {
          SGPP::base::DataVector p_train_private(this->dim);
          this->data->getRow(i, p_train_private);
          float_t value_train = opEval->eval(alpha_learn, p_train_private);

          // when there is only one baselearner actually, we do as following, just use normal classify to get the value
          if (this->numBaseLearners == 1)
            algorithmValueTrain.set(i, count, value_train);
          else if (count == 0)
            algorithmValueTrain.set(i, count, hypoweight * hValue(value_train));
          // each column is the sum value of baselearner respect to the column index
          else
            algorithmValueTrain.set(i, count, algorithmValueTrain.get(i, count - 1) + hypoweight * hValue(value_train));
        }

        // for testing data
        #pragma omp parallel for schedule(static)

        for (size_t i = 0; i < testData.getNrows(); i++) {
          SGPP::base::DataVector p_test_private(this->dim);
          testData.getRow(i, p_test_private);
          float_t value_test = opEval->eval(alpha_learn, p_test_private);

          // when there is only one baselearner actually, we do as following, just use normal classify to get the value
          if (this->numBaseLearners == 1)
            algorithmValueTest.set(i, count, value_test);
          else if (count == 0)
            algorithmValueTest.set(i, count, hypoweight * hValue(value_test));
          // each column is the sum value of baselearner respect to the column index
          else
            algorithmValueTest.set(i, count, algorithmValueTest.get(i, count - 1) + hypoweight * hValue(value_test));
        }

        float_t helper;

        for (size_t i = 0; i < this->numData; i++) {
          // helper = weight.get(i)*exp(-hypoWeight.get(count) * newclasses.get(i)*this->classes->get(i));
          if (newclasses.get(i) == this->classes->get(i))
            helper = weight.get(i) * exp(-hypoWeight.get(count));
          else
            helper = weight.get(i) * exp(hypoWeight.get(count));

          tmpweight.set(i, helper);
        }

        // normalization constant, this expression equals to normalizer = 2 * sqrt((weightError.get(count)) * (1.0 - weightError.get(count)));
        float_t normalizer = tmpweight.sum();

        // get new weights vector
        // SGPP::base::DataVector tmpweighthelp(this->numData);

        // tmpweighthelp = tmpweight;
        tmpweight.mult(1.0 / normalizer);
        weight = tmpweight;
        // for (size_t i = 0; i < this->numData; i++)
        // {
        //  weight.set(i, ((count+1)*tmpweighthelp.get(i) + tmpweight.get(i)) / (count + 2));
        // }

        if (count < this->numBaseLearners - 1 && this->refinement) {
          //reset the grid to the regular grid
          if (this->type == 1) {
            this->grid = SGPP::base::Grid::createLinearGrid(this->dim);
            std::cout << std::endl;
            std::cout << "Reset to the regular LinearGrid" << std::endl;
          } else if (this->type == 2) {
            this->grid = SGPP::base::Grid::createLinearBoundaryGrid(this->dim);
            std::cout << std::endl;
            std::cout << "Reset to the regular LinearBoundaryGrid" << std::endl;
          } else if (this->type == 3) {
            this->grid = SGPP::base::Grid::createModLinearGrid(this->dim);
            std::cout << std::endl;
            std::cout << "Reset to the regular ModLinearGrid" << std::endl;
          }
          // should not happen because this exception should have been thrown some lines upwards!
          else {
            throw new SGPP::base::operation_exception("AlgorithmAdaboost : Only 1 or 2 or 3 are supported gridType(1 = Linear Grid, 2 = LinearL0Boundary Grid, 3 = ModLinear Grid)!");
          }

          SGPP::base::GridGenerator* gridGen = this->grid->createGridGenerator();
          gridGen->regular(this->level);
          std::cout << std::endl;
          delete gridGen;
        }
      }

      delete opEval;
    }

    void AlgorithmAdaBoostBase::doRealAdaBoost(SGPP::base::DataMatrix& weights, SGPP::base::DataMatrix& testData, SGPP::base::DataMatrix& algorithmValueTrain, SGPP::base::DataMatrix& algorithmValueTest) {
      SGPP::base::DataVector weight(this->numData);
      weight.setAll(1.0 / float_t(this->numData));
      SGPP::base::OperationEval* opEval = SGPP::op_factory::createOperationEval(*this->grid);
      // to store certain train data point
      SGPP::base::DataVector p_train(this->dim);
      // to store certain train data point
      SGPP::base::DataVector p_test(this->dim);

      SGPP::base::DataVector tmpweight(this->numData);

      for (size_t count = 0; count < this->numBaseLearners; count++) {
        (this->actualBaseLearners)++;
        std::cout << std::endl;
        std::cout << "This is the " << this->actualBaseLearners << "th weak learner." << std::endl;
        std::cout << std::endl;

        // create coefficient vector
        SGPP::base::DataVector alpha_train(this->gridPoint);
        SGPP::base::DataVector alpha_learn(this->gridPoint);
        std::cout << "gridPoint: " << this->gridPoint << std::endl;

        if (this->maxGridPoint->get(count) < this->gridPoint)
          this->maxGridPoint->set(count, (float_t)this->gridPoint);

        if (!this->refinement) {
          if (count == 0)
            this->sumGridPoint->set(count, (float_t)gridPoint);
          else
            this->sumGridPoint->set(count, this->sumGridPoint->get(count - 1) + (float_t)gridPoint);
        }

        alpha_train.setAll(0.0);
        weights.setColumn(count, weight);

        bool final_step = false;

        if (this->refinement == 0)
          final_step = true;

        // calculate alpha
        alphaSolver(this->lamb, weight, alpha_train, final_step);

        if (this->refinement) {
          doRefinement(alpha_train, weight, count + 1);
          opEval = SGPP::op_factory::createOperationEval(*this->grid);
          alpha_learn.resizeZero(alpha_train.getSize());
        }

        //set the alpha for testing data(copy of alpha for training data)
        alpha_learn.copyFrom(alpha_train);

        // calculate the algorithm value of the testing data and training data
        // for training data
        #pragma omp parallel for schedule(static)

        for (size_t i = 0; i < numData; i++) {
          SGPP::base::DataVector p_train_private(this->dim);
          this->data->getRow(i, p_train_private);
          float_t value_train = opEval->eval(alpha_learn, p_train_private);
          float_t helper = weight.get(i) * exp(-this->classes->get(i) * value_train);
          tmpweight.set(i, helper);

          // when there is only one baselearner actually, we do as following, just use normal classify to get the value
          if (this->numBaseLearners == 1)
            algorithmValueTrain.set(i, count, 0.5 * value_train); // 0.5 as the coefficient (the original algorithm)
          else if (count == 0)
            algorithmValueTrain.set(i, count, 0.5 * value_train);
          // each column is the sum value of baselearner respect to the column index
          else
            algorithmValueTrain.set(i, count, algorithmValueTrain.get(i, count - 1) + 0.5 * value_train);
        }

        // normalize weight
        // update the weight
        float_t normalizer = tmpweight.sum();
        tmpweight.mult(1.0 / normalizer);
        weight = tmpweight;

        // for testing data
        #pragma omp parallel for schedule(static)

        for (size_t i = 0; i < testData.getNrows(); i++) {
          SGPP::base::DataVector p_test_private(this->dim);
          testData.getRow(i, p_test_private);
          float_t value_test = opEval->eval(alpha_learn, p_test_private);

          // when there is only one baselearner actually, we do as following, just use normal classify to get the value
          if (this->numBaseLearners == 1)
            algorithmValueTest.set(i, count, 0.5 * value_test);
          else if (count == 0)
            algorithmValueTest.set(i, count, 0.5 * value_test);
          // each column is the sum value of baselearner respect to the column index
          else
            algorithmValueTest.set(i, count, algorithmValueTest.get(i, count - 1) + 0.5 * value_test);
        }

        if (count < this->numBaseLearners - 1 && this->refinement) {
          //reset the grid to the regular grid
          if (this->type == 1) {
            this->grid = SGPP::base::Grid::createLinearGrid(this->dim);
            std::cout << std::endl;
            std::cout << "Reset to the regular LinearGrid" << std::endl;
          } else if (this->type == 2) {
            this->grid = SGPP::base::Grid::createLinearBoundaryGrid(this->dim);
            std::cout << std::endl;
            std::cout << "Reset to the regular LinearBoundaryGrid" << std::endl;
          } else if (this->type == 3) {
            this->grid = SGPP::base::Grid::createModLinearGrid(this->dim);
            std::cout << std::endl;
            std::cout << "Reset to the regular ModLinearGrid" << std::endl;
          }
          // should not happen because this exception should have been thrown some lines upwards!
          else {
            throw new SGPP::base::operation_exception("AlgorithmAdaboost : Only 1 or 2 or 3 are supported gridType(1 = Linear Grid, 2 = LinearL0Boundary Grid, 3 = ModLinear Grid)!");
          }

          SGPP::base::GridGenerator* gridGen = this->grid->createGridGenerator();
          gridGen->regular(this->level);
          std::cout << std::endl;
          delete gridGen;
        }
      }

      delete opEval;
    }

    void AlgorithmAdaBoostBase::doAdaBoostR2(SGPP::base::DataMatrix& weights, SGPP::base::DataMatrix& testData, SGPP::base::DataMatrix& algorithmValueTrain, SGPP::base::DataMatrix& algorithmValueTest, std::string lossFucType) {
      if (lossFucType != "linear" && lossFucType != "square" && lossFucType != "exponential") {
        throw new SGPP::base::operation_exception("AlgorithmAdaBoostBase::doAdaBoostR2 : An unknown loss function type was specified!");
      }

      SGPP::base::DataVector weight(this->numData);
      weight.setAll(1.0 / float_t(this->numData));
      SGPP::base::OperationEval* opEval = SGPP::op_factory::createOperationEval(*this->grid);
      // to store certain train data point
      SGPP::base::DataVector p_train(this->dim);
      // to store certain train data point
      SGPP::base::DataVector p_test(this->dim);
      SGPP::base::DataVector tmpweight(this->numData);

      // to store the loss
      SGPP::base::DataVector loss(this->numData);
      // to store the loss function value
      SGPP::base::DataVector lossFuc(this->numData);
      // to store the prediction training values
      SGPP::base::DataVector value_train(this->numData);
      // to store the prediction testing values
      SGPP::base::DataVector value_test(testData.getNrows());
      float_t maxloss;
      float_t meanloss;
      SGPP::base::DataVector beta(this->numBaseLearners); //[0,1]
      SGPP::base::DataVector logBetaSumR(this->numBaseLearners); //[0,1]

      for (size_t count = 0; count < this->numBaseLearners; count++) {
        (this->actualBaseLearners)++;
        std::cout << std::endl;
        std::cout << "This is the " << this->actualBaseLearners << "th weak learner." << std::endl;
        std::cout << std::endl;

        // create coefficient vector
        SGPP::base::DataVector alpha_train(this->gridPoint);
        SGPP::base::DataVector alpha_learn(this->gridPoint);
        std::cout << "gridPoint: " << this->gridPoint << std::endl;

        if (this->maxGridPoint->get(count) < this->gridPoint)
          this->maxGridPoint->set(count, (float_t)this->gridPoint);

        if (!this->refinement) {
          if (count == 0)
            this->sumGridPoint->set(count, (float_t)gridPoint);
          else
            this->sumGridPoint->set(count, this->sumGridPoint->get(count - 1) + (float_t)gridPoint);
        }

        alpha_train.setAll(0.0);
        weights.setColumn(count, weight);

        bool final_step = false;

        if (this->refinement == 0)
          final_step = true;

        // calculate alpha
        alphaSolver(this->lamb, weight, alpha_train, final_step);

        if (this->refinement) {
          doRefinement(alpha_train, weight, count + 1);
          opEval = SGPP::op_factory::createOperationEval(*this->grid);
          alpha_learn.resizeZero(alpha_train.getSize());
        }

        //set the alpha for testing data(copy of alpha for training data)
        alpha_learn.copyFrom(alpha_train);

        // calculate the algorithm value of the testing data and training data
        // for training data
        #pragma omp parallel for schedule(static)

        for (size_t i = 0; i < numData; i++) {
          SGPP::base::DataVector p_train_private(this->dim);
          this->data->getRow(i, p_train_private);
          value_train.set(i, opEval->eval(alpha_learn, p_train_private));
          loss.set(i, this->classes->get(i) - value_train.get(i));
        }

        loss.abs();
        maxloss = loss.max();

        if (lossFucType == "linear") {
          lossFuc = loss;
          lossFuc.mult(1 / maxloss);
        } else if (lossFucType == "square") {
          lossFuc = loss;
          lossFuc.mult(1 / maxloss);
          lossFuc.sqr();
        } else if (lossFucType == "exponential") {
          loss.mult(-1 / maxloss);
          #pragma omp parallel for schedule(static)

          for (size_t i = 0; i < numData; i++)
            lossFuc.set(i, 1 - exp(loss.get(i)));
        } else
          throw new SGPP::base::operation_exception("AlgorithmAdaBoostBase::doAdaBoostR2 : An unknown loss function type was specified!");

        meanloss = lossFuc.dotProduct(weight);

        // to judge the regressor whether match the simple condition of AdaboostR2
        if (meanloss >= 0.50) {
          std::cout << std::endl << "The average loss exceeds 0.5 after " << count + 1 << " iterations" << std::endl;
          (this->actualBaseLearners)--;
          break;
        }

        beta.set(count, meanloss / (1 - meanloss));

        SGPP::base::DataVector TrValueHelper(this->numData);
        float_t loghelp = log(1 / beta.get(count));

        if (count == 0) {
          logBetaSumR.set(count, loghelp);
          algorithmValueTrain.setColumn(count, value_train);
        }
        // each column is the sum value of baselearner respect to the column index
        else {
          logBetaSumR.set(count, logBetaSumR.get(count - 1) + loghelp);
          algorithmValueTrain.getColumn(count - 1, TrValueHelper);
          TrValueHelper.mult(logBetaSumR.get(count - 1));
          TrValueHelper.axpy(loghelp, value_train);
          TrValueHelper.mult(1 / logBetaSumR.get(count));
          algorithmValueTrain.setColumn(count, TrValueHelper);
        }

        #pragma omp parallel for schedule(static)

        for (size_t i = 0; i < numData; i++)
          tmpweight.set(i, std::pow(beta.get(count), 1 - lossFuc.get(i)));

        tmpweight.componentwise_mult(weight);
        // normalize weight
        // update the weight
        float_t normalizer = tmpweight.sum();
        tmpweight.mult(1.0 / normalizer);
        weight = tmpweight;

        // for testing data
        #pragma omp parallel for schedule(static)

        for (size_t i = 0; i < testData.getNrows(); i++) {
          SGPP::base::DataVector p_test_private(this->dim);
          testData.getRow(i, p_test_private);
          value_test.set(i, opEval->eval(alpha_learn, p_test_private));
        }

        SGPP::base::DataVector TeValueHelper(testData.getNrows());

        if (count == 0) {
          algorithmValueTest.setColumn(count, value_test);
        }
        // each column is the sum value of baselearner respect to the column index
        else {
          algorithmValueTest.getColumn(count - 1, TeValueHelper);
          TeValueHelper.mult(logBetaSumR.get(count - 1));
          TeValueHelper.axpy(loghelp, value_test);
          TeValueHelper.mult(1 / logBetaSumR.get(count));
          algorithmValueTest.setColumn(count, TeValueHelper);
        }

        if (count < this->numBaseLearners - 1 && this->refinement) {
          //reset the grid to the regular grid
          if (this->type == 1) {
            this->grid = SGPP::base::Grid::createLinearGrid(this->dim);
            std::cout << std::endl;
            std::cout << "Reset to the regular LinearGrid" << std::endl;
          } else if (this->type == 2) {
            this->grid = SGPP::base::Grid::createLinearBoundaryGrid(this->dim);
            std::cout << std::endl;
            std::cout << "Reset to the regular LinearBoundaryGrid" << std::endl;
          } else if (this->type == 3) {
            this->grid = SGPP::base::Grid::createModLinearGrid(this->dim);
            std::cout << std::endl;
            std::cout << "Reset to the regular ModLinearGrid" << std::endl;
          }
          // should not happen because this exception should have been thrown some lines upwards!
          else {
            throw new SGPP::base::operation_exception("AlgorithmAdaboost : Only 1 or 2 or 3 are supported gridType(1 = Linear Grid, 2 = LinearL0Boundary Grid, 3 = ModLinear Grid)!");
          }

          SGPP::base::GridGenerator* gridGen = this->grid->createGridGenerator();
          gridGen->regular(this->level);
          std::cout << std::endl;
          delete gridGen;
        }
      }

      delete opEval;
    }

    void AlgorithmAdaBoostBase::doAdaBoostRT(SGPP::base::DataMatrix& weights, SGPP::base::DataMatrix& testData, SGPP::base::DataMatrix& algorithmValueTrain, SGPP::base::DataMatrix& algorithmValueTest, float_t Tvalue, std::string powerType) {
      if (Tvalue >= 1 || Tvalue <= 0) {
        throw new SGPP::base::operation_exception("AlgorithmAdaBoostBase::doAdaBoostRT : the Tvalue must lie between 0 and 1!");
      }

      if (powerType != "linear" && powerType != "square" && powerType != "cubic") {
        throw new SGPP::base::operation_exception("AlgorithmAdaBoostBase::doAdaBoostRT : An unknown power type was specified!");
      }

      SGPP::base::DataVector weight(this->numData);
      weight.setAll(1.0 / float_t(this->numData));
      SGPP::base::OperationEval* opEval = SGPP::op_factory::createOperationEval(*this->grid);
      // to store certain train data point
      SGPP::base::DataVector p_train(this->dim);
      // to store certain train data point
      SGPP::base::DataVector p_test(this->dim);
      SGPP::base::DataVector tmpweight(this->numData);

      // to store the absolute relative error
      SGPP::base::DataVector ARE(this->numData);
      // to store the prediction training values
      SGPP::base::DataVector value_train(this->numData);
      // to store the prediction testing values
      SGPP::base::DataVector value_test(testData.getNrows());

      SGPP::base::DataVector beta(this->numBaseLearners);
      SGPP::base::DataVector logBetaSumR(this->numBaseLearners);
      float_t errorRate;

      for (size_t count = 0; count < this->numBaseLearners; count++) {
        (this->actualBaseLearners)++;
        std::cout << std::endl;
        std::cout << "This is the " << this->actualBaseLearners << "th weak learner." << std::endl;
        std::cout << std::endl;

        // create coefficient vector
        SGPP::base::DataVector alpha_train(this->gridPoint);
        SGPP::base::DataVector alpha_learn(this->gridPoint);
        std::cout << "gridPoint: " << this->gridPoint << std::endl;

        if (this->maxGridPoint->get(count) < this->gridPoint)
          this->maxGridPoint->set(count, (float_t)this->gridPoint);

        if (!this->refinement) {
          if (count == 0)
            this->sumGridPoint->set(count, (float_t)gridPoint);
          else
            this->sumGridPoint->set(count, this->sumGridPoint->get(count - 1) + (float_t)gridPoint);
        }

        alpha_train.setAll(0.0);
        weights.setColumn(count, weight);

        bool final_step = false;

        if (this->refinement == 0)
          final_step = true;

        // calculate alpha
        alphaSolver(this->lamb, weight, alpha_train, final_step);

        if (this->refinement) {
          doRefinement(alpha_train, weight, count + 1);
          opEval = SGPP::op_factory::createOperationEval(*this->grid);
          alpha_learn.resizeZero(alpha_train.getSize());
        }

        //set the alpha for testing data(copy of alpha for training data)
        alpha_learn.copyFrom(alpha_train);

        // calculate the algorithm value of the testing data and training data
        // for training data
        errorRate = 0;
        #pragma omp parallel for schedule(static)

        for (size_t i = 0; i < numData; i++) {
          SGPP::base::DataVector p_train_private(this->dim);
          this->data->getRow(i, p_train_private);
          value_train.set(i, opEval->eval(alpha_learn, p_train_private));
          ARE.set(i, std::abs((this->classes->get(i) - value_train.get(i)) / this->classes->get(i)));

          if (ARE.get(i) > Tvalue)
            errorRate += weight.get(i);
        }

        if (powerType == "linear")
          beta.set(count, errorRate);
        else if (powerType == "square")
          beta.set(count, errorRate * errorRate);
        else if (powerType == "cubic")
          beta.set(count, errorRate * errorRate * errorRate);
        else
          throw new SGPP::base::operation_exception("AlgorithmAdaBoostBase::doAdaBoostRT : An unknown power type was specified!");

        SGPP::base::DataVector TrValueHelper(this->numData);
        float_t loghelp = log(1 / beta.get(count));

        if (count == 0) {
          logBetaSumR.set(count, loghelp);
          algorithmValueTrain.setColumn(count, value_train);
        }
        // each column is the sum value of baselearner respect to the column index
        else {
          logBetaSumR.set(count, logBetaSumR.get(count - 1) + loghelp);
          algorithmValueTrain.getColumn(count - 1, TrValueHelper);
          TrValueHelper.mult(logBetaSumR.get(count - 1));
          TrValueHelper.axpy(loghelp, value_train);
          TrValueHelper.mult(1 / logBetaSumR.get(count));
          algorithmValueTrain.setColumn(count, TrValueHelper);
        }

        #pragma omp parallel for schedule(static)

        for (size_t i = 0; i < numData; i++) {
          if (ARE.get(i) <= Tvalue)
            tmpweight.set(i, weight.get(i)*beta.get(count));
          else
            tmpweight.set(i, weight.get(i));
        }

        // normalize weight
        // update the weight
        float_t normalizer = tmpweight.sum();
        tmpweight.mult(1.0 / normalizer);
        weight = tmpweight;

        // for testing data
        #pragma omp parallel for schedule(static)

        for (size_t i = 0; i < testData.getNrows(); i++) {
          SGPP::base::DataVector p_test_private(this->dim);
          testData.getRow(i, p_test_private);
          value_test.set(i, opEval->eval(alpha_learn, p_test_private));
        }

        SGPP::base::DataVector TeValueHelper(testData.getNrows());

        if (count == 0) {
          algorithmValueTest.setColumn(count, value_test);
        }
        // each column is the sum value of baselearner respect to the column index
        else {
          algorithmValueTest.getColumn(count - 1, TeValueHelper);
          TeValueHelper.mult(logBetaSumR.get(count - 1));
          TeValueHelper.axpy(loghelp, value_test);
          TeValueHelper.mult(1 / logBetaSumR.get(count));
          algorithmValueTest.setColumn(count, TeValueHelper);
        }

        if (count < this->numBaseLearners - 1 && this->refinement) {
          //reset the grid to the regular grid
          if (this->type == 1) {
            this->grid = SGPP::base::Grid::createLinearGrid(this->dim);
            std::cout << std::endl;
            std::cout << "Reset to the regular LinearGrid" << std::endl;
          } else if (this->type == 2) {
            this->grid = SGPP::base::Grid::createLinearBoundaryGrid(this->dim);
            std::cout << std::endl;
            std::cout << "Reset to the regular LinearBoundaryGrid" << std::endl;
          } else if (this->type == 3) {
            this->grid = SGPP::base::Grid::createModLinearGrid(this->dim);
            std::cout << std::endl;
            std::cout << "Reset to the regular ModLinearGrid" << std::endl;
          }
          // should not happen because this exception should have been thrown some lines upwards!
          else {
            throw new SGPP::base::operation_exception("AlgorithmAdaboost : Only 1 or 2 or 3 are supported gridType(1 = Linear Grid, 2 = LinearL0Boundary Grid, 3 = ModLinear Grid)!");
          }

          SGPP::base::GridGenerator* gridGen = this->grid->createGridGenerator();
          gridGen->regular(this->level);
          std::cout << std::endl;
          delete gridGen;
        }
      }

      delete opEval;
    }

    void AlgorithmAdaBoostBase::eval(SGPP::base::DataMatrix& testData, SGPP::base::DataMatrix& algorithmValueTrain, SGPP::base::DataMatrix& algorithmValueTest) {
      SGPP::base::DataMatrix weightsMatrix(this->numData, this->numBaseLearners);
      weightsMatrix.setAll(0.0);

      if (this->boostMode == 1) {
        SGPP::base::DataVector theHypoWeight(this->numBaseLearners);
        SGPP::base::DataVector theWeightError(this->numBaseLearners);
        SGPP::base::DataMatrix decisionMatrix(this->numData, this->numBaseLearners);
        theHypoWeight.setAll(0.0);
        theWeightError.setAll(0.0);
        decisionMatrix.setAll(0.0);
        doDiscreteAdaBoost(theHypoWeight, theWeightError, weightsMatrix, decisionMatrix, testData, algorithmValueTrain, algorithmValueTest);
      } else if (this->boostMode == 2)
        doRealAdaBoost(weightsMatrix, testData, algorithmValueTrain, algorithmValueTest);
      else {
        throw new SGPP::base::operation_exception("AlgorithmAdaboost : Only 1 or 2 for the boost mode(1 = Discrete Adaboost, 2 = Real Adaboost)!");
      }
    }

    void AlgorithmAdaBoostBase::classif(SGPP::base::DataMatrix& testData, SGPP::base::DataVector& algorithmClassTrain, SGPP::base::DataVector& algorithmClassTest, SGPP::base::DataMatrix& algorithmValueTrain, SGPP::base::DataMatrix& algorithmValueTest) {
      eval(testData, algorithmValueTrain, algorithmValueTest);

      for (size_t i = 0; i < this->numData; i++) {
        algorithmClassTrain.set(i, hValue(algorithmValueTrain.get(i, this->actualBaseLearners - 1)));
      }

      for (size_t i = 0; i < testData.getNrows(); i++) {
        algorithmClassTest.set(i, hValue(algorithmValueTest.get(i, this->actualBaseLearners - 1)));
      }
    }

    void AlgorithmAdaBoostBase::getAccuracy(SGPP::base::DataMatrix& testData, SGPP::base::DataVector& testDataClass, float_t* accuracy_train, float_t* accuracy_test) {
      /* get the accuracy */
      size_t right_test = 0;
      size_t right_train = 0;
      SGPP::base::DataVector classTrain(this->numData);
      SGPP::base::DataVector classTest(testDataClass.getSize());
      SGPP::base::DataMatrix valueTrain(this->numData, this->numBaseLearners);
      SGPP::base::DataMatrix valueTest(testDataClass.getSize(), this->numBaseLearners);
      classif(testData, classTrain, classTest, valueTrain, valueTest);

      // for training data
      for (size_t i = 0; i < this->numData; i++) {
        if (classTrain.get(i) == this->classes->get(i))
          right_train = right_train + 1;
      }

      *accuracy_train = float_t(right_train) / float_t(this->numData);

      // for testing data
      for (size_t i = 0; i < testData.getNrows(); i++) {
        if (classTest.get(i) == testDataClass.get(i))
          right_test = right_test + 1;
      }

      *accuracy_test = float_t(right_test) / float_t(classTest.getSize());
    }

    void AlgorithmAdaBoostBase::getROC(SGPP::base::DataMatrix& validationData, SGPP::base::DataVector& validationDataClass, float_t* acc, float_t* sensitivity, float_t* specificity, float_t* precision, float_t* recall, float_t* fOneScore) {
      size_t truePos = 0;
      size_t predictPos = 0;
      size_t trueNeg = 0;
      size_t predictNeg = 0;
      size_t actualPos = 0;
      size_t actualNeg = 0;

      SGPP::base::DataVector classTrain(this->numData);
      SGPP::base::DataVector classValidation(validationDataClass.getSize());
      SGPP::base::DataMatrix valueTrain(this->numData, this->numBaseLearners);
      SGPP::base::DataMatrix valueValidation(validationDataClass.getSize(), this->numBaseLearners);
      classif(validationData, classTrain, classValidation, valueTrain, valueValidation);

      // for validation data
      for (size_t i = 0; i < validationData.getNrows(); i++) {
        if (validationDataClass.get(i) == this->labelOne)
          actualPos += 1;
        else
          actualNeg += 1;

        if (classValidation.get(i) == this->labelOne)
          predictPos += 1;
        else
          predictNeg += 1;

        if (classValidation.get(i) == validationDataClass.get(i) && validationDataClass.get(i) == this->labelOne)
          truePos += 1;

        if (classValidation.get(i) == validationDataClass.get(i) && validationDataClass.get(i) == this->labelTwo)
          trueNeg += 1;
      }

      *acc = float_t(truePos + trueNeg) / float_t(classValidation.getSize());
      *sensitivity = float_t(truePos) / float_t(actualPos);
      *specificity = float_t(trueNeg) / float_t(actualNeg);
      *precision = float_t(truePos) / float_t(predictPos);
      *recall = *sensitivity;
      *fOneScore = 2 * (*precision) * (*recall) / ((*precision) + (*recall));
    }

    void AlgorithmAdaBoostBase::getAccuracyBL(SGPP::base::DataMatrix& testData, SGPP::base::DataVector& testDataClass, SGPP::base::DataMatrix& algorithmValueTrain, SGPP::base::DataMatrix& algorithmValueTest, float_t* accuracy_train, float_t* accuracy_test, size_t yourBaseLearner) {
      size_t right_test = 0;
      size_t right_train = 0;

      // for training data
      for (size_t i = 0; i < this->numData; i++) {
        if (hValue(algorithmValueTrain.get(i, yourBaseLearner - 1)) == this->classes->get(i))
          right_train = right_train + 1;
      }

      *accuracy_train = float_t(right_train) / float_t(this->numData);

      // for testing data
      for (size_t i = 0; i < testData.getNrows(); i++) {
        if (hValue(algorithmValueTest.get(i, yourBaseLearner - 1)) == testDataClass.get(i))
          right_test = right_test + 1;
      }

      *accuracy_test = float_t(right_test) / float_t(testDataClass.getSize());
    }

    void AlgorithmAdaBoostBase::doRefinement(SGPP::base::DataVector& alpha_ada, SGPP::base::DataVector& weight_ada, size_t curBaseLearner) {
      bool final_ada = false;

      for (size_t adaptiveStep = 1; adaptiveStep <= this->refineTimes; adaptiveStep++) {

        SGPP::base::GridGenerator* myGenerator = this->grid->createGridGenerator();
        size_t refineNumber;

        if (this->refineMode == 1) {
          if (this->numOfAda > this->grid->getSize())
            refineNumber = this->grid->getSize();
          else
            refineNumber = this->numOfAda;
        } else if (this->refineMode == 2) {
          refineNumber = (size_t)(this->perOfAda * (float_t)(this->grid->getSize()));

          //force to refine at least one point
          if (refineNumber == 0)
            refineNumber = 1;
        }
        // should not happen because this exception should have been thrown some lines upwards!
        else {
          throw new SGPP::base::operation_exception("AlgorithmAdaBoost : Only 1 or 2 are supported refine mode(1 : use grid point number, 2: use grid point percentage)!");
        }

        SGPP::base::SurplusRefinementFunctor* myRefineFunc = new SGPP::base::SurplusRefinementFunctor(&alpha_ada, refineNumber, 0.0);
        myGenerator->refine(myRefineFunc);
        delete myRefineFunc;
        delete myGenerator;

        SGPP::base::GridStorage* gridStorage_ada = this->grid->getStorage();
        size_t gridPts = gridStorage_ada->size();

        std::cout << std::endl;
        std::cout << "Refinement time step: " << adaptiveStep << ", new grid size: " << gridPts << ", refined number of grid points: " << refineNumber << std::endl;

        if (adaptiveStep == this->refineTimes) {
          final_ada = true;

          if (curBaseLearner == 1) {
            this->maxGridPoint->set(curBaseLearner - 1, (float_t)gridPts);
            this->sumGridPoint->set(curBaseLearner - 1, (float_t)gridPts);
          } else {
            if (gridPts > this->maxGridPoint->get(curBaseLearner - 2))
              this->maxGridPoint->set(curBaseLearner - 1, (float_t)gridPts);
            else
              this->maxGridPoint->set(curBaseLearner - 1, this->maxGridPoint->get(curBaseLearner - 2));

            this->sumGridPoint->set(curBaseLearner - 1, this->sumGridPoint->get(curBaseLearner - 2) + (float_t)gridPts);
          }
        }

        // extend alpha vector (new entries uninitialized)
        alpha_ada.resizeZero(gridPts);

        // calculate new alpha
        alphaSolver(this->lamb, weight_ada, alpha_ada, final_ada);
      }
    }

    float_t AlgorithmAdaBoostBase::hValue(float_t realValue) {
      if (realValue >= this->threshold) {
        if (labelOne > labelTwo)
          return labelOne;
        else
          return labelTwo;
      } else {
        if (labelOne > labelTwo)
          return labelTwo;
        else
          return labelOne;
      }
    }

    size_t AlgorithmAdaBoostBase::getActualBL() {
      return this->actualBaseLearners;
    }

    size_t AlgorithmAdaBoostBase::getMeanGridPoint(size_t baselearner) {
      size_t mean = (size_t)(this->sumGridPoint->get(baselearner - 1) / (float_t)baselearner);
      return mean;
    }

    size_t AlgorithmAdaBoostBase::getMaxGridPoint(size_t baselearner) {
      size_t max = (size_t)(this->maxGridPoint->get(baselearner - 1));
      return max;
    }

    size_t AlgorithmAdaBoostBase::getSumGridPoint(size_t baselearner) {
      size_t sum = (size_t)(this->sumGridPoint->get(baselearner - 1));
      return sum;
    }
  }
}