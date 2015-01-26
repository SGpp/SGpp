/******************************************************************************
* Copyright (C) 2011 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
******************************************************************************/
// @author Zhongwen Song (songz@in.tum.de)
// @author Benjamin (pehersto@in.tum.de)
// @auther Alexander Heinecke (alexander.heinecke@mytum.de)

#include <sgpp/datadriven/algorithm/AlgorithmAdaBoostBase.hpp>
#include <sgpp/base/exception/operation_exception.hpp>
#include <sgpp/base/operation/BaseOpFactory.hpp>

namespace sg {
  namespace datadriven {
    AlgorithmAdaBoostBase::AlgorithmAdaBoostBase(sg::base::Grid& SparseGrid, size_t gridType, sg::base::HashGenerator::level_t gridLevel, sg::base::DataMatrix& trainData, sg::base::DataVector& trainDataClass, size_t NUM, double lambda, size_t IMAX, double eps, size_t IMAX_final, double eps_final, double firstLabel, double secondLabel, double threshold, double maxLambda, double minLambda, size_t searchNum, bool refine, size_t refineMode, size_t refineNum, size_t numberOfAda, double percentOfAda, size_t mode) {
      if (refine && (gridType != 1 && gridType != 2 && gridType != 3)) {
        throw new sg::base::operation_exception("AlgorithmAdaBoostBase : Only 1 or 2 or 3 are supported gridType(1 = Linear Grid, 2 = LinearBoundary Grid, 3 = ModLinear Grid)!");
      }

      if (refine && (percentOfAda >= 1.0 || percentOfAda <= 0.0)) {
        throw new sg::base::operation_exception("AlgorithmAdaBoostBase : Only number between 0 and 1 is the supported percent to Adaptive!");
      }

      if (refineMode != 1 && refineMode != 2) {
        throw new sg::base::operation_exception("AlgorithmAdaBoostBase : Only 1 or 2 are supported refine mode(1 : use grid point number, 2: use grid point percentage)!");
      }

      sg::base::GridStorage* gridStorage = SparseGrid.getStorage();
      this->grid = &SparseGrid;
      this->type = gridType;
      this->gridPoint = gridStorage->size();
      this->level = static_cast<sg::base::HashGenerator::level_t>(gridLevel);
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
        this->lambStepsize = (log((double)maxLambda) - log((double)minLambda)) / 2;
      else
        this->lambStepsize = (log((double)maxLambda) - log((double)minLambda)) / ((double)searchNum - 1);

      this->actualBaseLearners = 0;
      this->refinement = refine;
      this->refineMode = refineMode;
      this->refineTimes = refineNum;
      this->numOfAda = numberOfAda;
      this->perOfAda = percentOfAda;
      this->maxGridPoint = new sg::base::DataVector(NUM);
      this->sumGridPoint = new sg::base::DataVector(NUM);
      this->boostMode = mode;
    }

    AlgorithmAdaBoostBase::~AlgorithmAdaBoostBase() {

    }

    void AlgorithmAdaBoostBase::doDiscreteAdaBoost(sg::base::DataVector& hypoWeight, sg::base::DataVector& weightError, sg::base::DataMatrix& weights, sg::base::DataMatrix& decision, sg::base::DataMatrix& testData, sg::base::DataMatrix& algorithmValueTrain, sg::base::DataMatrix& algorithmValueTest) {
      sg::base::DataVector weight(this->numData);
      weight.setAll(1.0 / double(this->numData));
      sg::base::OperationEval* opEval = sg::op_factory::createOperationEval(*this->grid);
      // to store certain train data point
      sg::base::DataVector p_train(this->dim);
      // to store certain train data point
      sg::base::DataVector p_test(this->dim);
      // create vector to store the hypothesis of the training data according to certain alpha vector(base learner)
      sg::base::DataVector newclasses(this->numData);
      // create vector to store the identity of comparing hypothesis of training data with class of
      sg::base::DataVector identity(this->numData);

      sg::base::DataVector tmpweight(this->numData);

      for (size_t count = 0; count < this->numBaseLearners; count++) {
        (this->actualBaseLearners)++;
        std::cout << std::endl;
        std::cout << "This is the " << this->actualBaseLearners << "th weak learner." << std::endl;
        std::cout << std::endl;

        // create coefficient vector
        sg::base::DataVector alpha_train(this->gridPoint);
        sg::base::DataVector alpha_learn(this->gridPoint);
        std::cout << "gridPoint: " << this->gridPoint << std::endl;

        if (this->maxGridPoint->get(count) < this->gridPoint)
          this->maxGridPoint->set(count, (double)this->gridPoint);

        if (!this->refinement) {
          if (count == 0)
            this->sumGridPoint->set(count, (double)gridPoint);
          else
            this->sumGridPoint->set(count, this->sumGridPoint->get(count - 1) + (double)gridPoint);
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
          opEval = sg::op_factory::createOperationEval(*this->grid);
          alpha_learn.resizeZero(alpha_train.getSize());
        }

        //set the alpha for testing data(copy of alpha for training data)
        alpha_learn.copyFrom(alpha_train);

        #pragma omp parallel for schedule(static)

        for (size_t i = 0; i < this->numData; i++) {
          sg::base::DataVector p_train_private(this->dim);
          this->data->getRow(i, p_train_private);
          double value_train = opEval->eval(alpha_train, p_train_private);
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
          double cur_lambda;
          double weighterror;
          double minWeightError = weightError.get(count);

          for (size_t it = 0; it < this->lambSteps; it++) {
            std::cout << std::endl;
            std::cout << "This is the " << it + 1 << "th search of " << this->actualBaseLearners << "th weak learner." << std::endl;
            std::cout << std::endl;
            alpha_train.setAll(0.0);
            cur_lambda = exp(this->lambLogMax - (double)it * this->lambStepsize);

            alphaSolver(cur_lambda, weight, alpha_train, true);

            for (size_t i = 0; i < this->numData; i++) {
              this->data->getRow(i, p_train);
              double value_seach = opEval->eval(alpha_train, p_train);
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
        double hypoweight;

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
          sg::base::DataVector p_train_private(this->dim);
          this->data->getRow(i, p_train_private);
          double value_train = opEval->eval(alpha_learn, p_train_private);

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
          sg::base::DataVector p_test_private(this->dim);
          testData.getRow(i, p_test_private);
          double value_test = opEval->eval(alpha_learn, p_test_private);

          // when there is only one baselearner actually, we do as following, just use normal classify to get the value
          if (this->numBaseLearners == 1)
            algorithmValueTest.set(i, count, value_test);
          else if (count == 0)
            algorithmValueTest.set(i, count, hypoweight * hValue(value_test));
          // each column is the sum value of baselearner respect to the column index
          else
            algorithmValueTest.set(i, count, algorithmValueTest.get(i, count - 1) + hypoweight * hValue(value_test));
        }

        double helper;

        for (size_t i = 0; i < this->numData; i++) {
          // helper = weight.get(i)*exp(-hypoWeight.get(count) * newclasses.get(i)*this->classes->get(i));
          if (newclasses.get(i) == this->classes->get(i))
            helper = weight.get(i) * exp(-hypoWeight.get(count));
          else
            helper = weight.get(i) * exp(hypoWeight.get(count));

          tmpweight.set(i, helper);
        }

        // normalization constant, this expression equals to normalizer = 2 * sqrt((weightError.get(count)) * (1.0 - weightError.get(count)));
        double normalizer = tmpweight.sum();

        // get new weights vector
        // sg::base::DataVector tmpweighthelp(this->numData);

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
            this->grid = sg::base::Grid::createLinearGrid(this->dim);
            std::cout << std::endl;
            std::cout << "Reset to the regular LinearGrid" << std::endl;
          } else if (this->type == 2) {
            this->grid = sg::base::Grid::createLinearTrapezoidBoundaryGrid(this->dim);
            std::cout << std::endl;
            std::cout << "Reset to the regular LinearTrapezoidBoundaryGrid" << std::endl;
          } else if (this->type == 3) {
            this->grid = sg::base::Grid::createModLinearGrid(this->dim);
            std::cout << std::endl;
            std::cout << "Reset to the regular ModLinearGrid" << std::endl;
          }
          // should not happen because this exception should have been thrown some lines upwards!
          else {
            throw new sg::base::operation_exception("AlgorithmAdaboost : Only 1 or 2 or 3 are supported gridType(1 = Linear Grid, 2 = LinearBoundary Grid, 3 = ModLinear Grid)!");
          }

          sg::base::GridGenerator* gridGen = this->grid->createGridGenerator();
          gridGen->regular(this->level);
          std::cout << std::endl;
          delete gridGen;
        }
      }

      delete opEval;
    }

    void AlgorithmAdaBoostBase::doRealAdaBoost(sg::base::DataMatrix& weights, sg::base::DataMatrix& testData, sg::base::DataMatrix& algorithmValueTrain, sg::base::DataMatrix& algorithmValueTest) {
      sg::base::DataVector weight(this->numData);
      weight.setAll(1.0 / double(this->numData));
      sg::base::OperationEval* opEval = sg::op_factory::createOperationEval(*this->grid);
      // to store certain train data point
      sg::base::DataVector p_train(this->dim);
      // to store certain train data point
      sg::base::DataVector p_test(this->dim);

      sg::base::DataVector tmpweight(this->numData);

      for (size_t count = 0; count < this->numBaseLearners; count++) {
        (this->actualBaseLearners)++;
        std::cout << std::endl;
        std::cout << "This is the " << this->actualBaseLearners << "th weak learner." << std::endl;
        std::cout << std::endl;

        // create coefficient vector
        sg::base::DataVector alpha_train(this->gridPoint);
        sg::base::DataVector alpha_learn(this->gridPoint);
        std::cout << "gridPoint: " << this->gridPoint << std::endl;

        if (this->maxGridPoint->get(count) < this->gridPoint)
          this->maxGridPoint->set(count, (double)this->gridPoint);

        if (!this->refinement) {
          if (count == 0)
            this->sumGridPoint->set(count, (double)gridPoint);
          else
            this->sumGridPoint->set(count, this->sumGridPoint->get(count - 1) + (double)gridPoint);
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
          opEval = sg::op_factory::createOperationEval(*this->grid);
          alpha_learn.resizeZero(alpha_train.getSize());
        }

        //set the alpha for testing data(copy of alpha for training data)
        alpha_learn.copyFrom(alpha_train);

        // calculate the algorithm value of the testing data and training data
        // for training data
        #pragma omp parallel for schedule(static)

        for (size_t i = 0; i < numData; i++) {
          sg::base::DataVector p_train_private(this->dim);
          this->data->getRow(i, p_train_private);
          double value_train = opEval->eval(alpha_learn, p_train_private);
          double helper = weight.get(i) * exp(-this->classes->get(i) * value_train);
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
        double normalizer = tmpweight.sum();
        tmpweight.mult(1.0 / normalizer);
        weight = tmpweight;

        // for testing data
        #pragma omp parallel for schedule(static)

        for (size_t i = 0; i < testData.getNrows(); i++) {
          sg::base::DataVector p_test_private(this->dim);
          testData.getRow(i, p_test_private);
          double value_test = opEval->eval(alpha_learn, p_test_private);

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
            this->grid = sg::base::Grid::createLinearGrid(this->dim);
            std::cout << std::endl;
            std::cout << "Reset to the regular LinearGrid" << std::endl;
          } else if (this->type == 2) {
            this->grid = sg::base::Grid::createLinearTrapezoidBoundaryGrid(this->dim);
            std::cout << std::endl;
            std::cout << "Reset to the regular LinearTrapezoidBoundaryGrid" << std::endl;
          } else if (this->type == 3) {
            this->grid = sg::base::Grid::createModLinearGrid(this->dim);
            std::cout << std::endl;
            std::cout << "Reset to the regular ModLinearGrid" << std::endl;
          }
          // should not happen because this exception should have been thrown some lines upwards!
          else {
            throw new sg::base::operation_exception("AlgorithmAdaboost : Only 1 or 2 or 3 are supported gridType(1 = Linear Grid, 2 = LinearBoundary Grid, 3 = ModLinear Grid)!");
          }

          sg::base::GridGenerator* gridGen = this->grid->createGridGenerator();
          gridGen->regular(this->level);
          std::cout << std::endl;
          delete gridGen;
        }
      }

      delete opEval;
    }

    void AlgorithmAdaBoostBase::doAdaBoostR2(sg::base::DataMatrix& weights, sg::base::DataMatrix& testData, sg::base::DataMatrix& algorithmValueTrain, sg::base::DataMatrix& algorithmValueTest, std::string lossFucType) {
      if (lossFucType != "linear" && lossFucType != "square" && lossFucType != "exponential") {
        throw new sg::base::operation_exception("AlgorithmAdaBoostBase::doAdaBoostR2 : An unknown loss function type was specified!");
      }

      sg::base::DataVector weight(this->numData);
      weight.setAll(1.0 / double(this->numData));
      sg::base::OperationEval* opEval = sg::op_factory::createOperationEval(*this->grid);
      // to store certain train data point
      sg::base::DataVector p_train(this->dim);
      // to store certain train data point
      sg::base::DataVector p_test(this->dim);
      sg::base::DataVector tmpweight(this->numData);

      // to store the loss
      sg::base::DataVector loss(this->numData);
      // to store the loss function value
      sg::base::DataVector lossFuc(this->numData);
      // to store the prediction training values
      sg::base::DataVector value_train(this->numData);
      // to store the prediction testing values
      sg::base::DataVector value_test(testData.getNrows());
      double maxloss;
      double meanloss;
      sg::base::DataVector beta(this->numBaseLearners); //[0,1]
      sg::base::DataVector logBetaSumR(this->numBaseLearners); //[0,1]

      for (size_t count = 0; count < this->numBaseLearners; count++) {
        (this->actualBaseLearners)++;
        std::cout << std::endl;
        std::cout << "This is the " << this->actualBaseLearners << "th weak learner." << std::endl;
        std::cout << std::endl;

        // create coefficient vector
        sg::base::DataVector alpha_train(this->gridPoint);
        sg::base::DataVector alpha_learn(this->gridPoint);
        std::cout << "gridPoint: " << this->gridPoint << std::endl;

        if (this->maxGridPoint->get(count) < this->gridPoint)
          this->maxGridPoint->set(count, (double)this->gridPoint);

        if (!this->refinement) {
          if (count == 0)
            this->sumGridPoint->set(count, (double)gridPoint);
          else
            this->sumGridPoint->set(count, this->sumGridPoint->get(count - 1) + (double)gridPoint);
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
          opEval = sg::op_factory::createOperationEval(*this->grid);
          alpha_learn.resizeZero(alpha_train.getSize());
        }

        //set the alpha for testing data(copy of alpha for training data)
        alpha_learn.copyFrom(alpha_train);

        // calculate the algorithm value of the testing data and training data
        // for training data
        #pragma omp parallel for schedule(static)

        for (size_t i = 0; i < numData; i++) {
          sg::base::DataVector p_train_private(this->dim);
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
          throw new sg::base::operation_exception("AlgorithmAdaBoostBase::doAdaBoostR2 : An unknown loss function type was specified!");

        meanloss = lossFuc.dotProduct(weight);

        // to judge the regressor whether match the simple condition of AdaboostR2
        if (meanloss >= 0.50) {
          std::cout << std::endl << "The average loss exceeds 0.5 after " << count + 1 << " iterations" << std::endl;
          (this->actualBaseLearners)--;
          break;
        }

        beta.set(count, meanloss / (1 - meanloss));

        sg::base::DataVector TrValueHelper(this->numData);
        double loghelp = log(1 / beta.get(count));

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
        double normalizer = tmpweight.sum();
        tmpweight.mult(1.0 / normalizer);
        weight = tmpweight;

        // for testing data
        #pragma omp parallel for schedule(static)

        for (size_t i = 0; i < testData.getNrows(); i++) {
          sg::base::DataVector p_test_private(this->dim);
          testData.getRow(i, p_test_private);
          value_test.set(i, opEval->eval(alpha_learn, p_test_private));
        }

        sg::base::DataVector TeValueHelper(testData.getNrows());

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
            this->grid = sg::base::Grid::createLinearGrid(this->dim);
            std::cout << std::endl;
            std::cout << "Reset to the regular LinearGrid" << std::endl;
          } else if (this->type == 2) {
            this->grid = sg::base::Grid::createLinearTrapezoidBoundaryGrid(this->dim);
            std::cout << std::endl;
            std::cout << "Reset to the regular LinearTrapezoidBoundaryGrid" << std::endl;
          } else if (this->type == 3) {
            this->grid = sg::base::Grid::createModLinearGrid(this->dim);
            std::cout << std::endl;
            std::cout << "Reset to the regular ModLinearGrid" << std::endl;
          }
          // should not happen because this exception should have been thrown some lines upwards!
          else {
            throw new sg::base::operation_exception("AlgorithmAdaboost : Only 1 or 2 or 3 are supported gridType(1 = Linear Grid, 2 = LinearBoundary Grid, 3 = ModLinear Grid)!");
          }

          sg::base::GridGenerator* gridGen = this->grid->createGridGenerator();
          gridGen->regular(this->level);
          std::cout << std::endl;
          delete gridGen;
        }
      }

      delete opEval;
    }

    void AlgorithmAdaBoostBase::doAdaBoostRT(sg::base::DataMatrix& weights, sg::base::DataMatrix& testData, sg::base::DataMatrix& algorithmValueTrain, sg::base::DataMatrix& algorithmValueTest, double Tvalue, std::string powerType) {
      if (Tvalue >= 1 || Tvalue <= 0) {
        throw new sg::base::operation_exception("AlgorithmAdaBoostBase::doAdaBoostRT : the Tvalue must lie between 0 and 1!");
      }

      if (powerType != "linear" && powerType != "square" && powerType != "cubic") {
        throw new sg::base::operation_exception("AlgorithmAdaBoostBase::doAdaBoostRT : An unknown power type was specified!");
      }

      sg::base::DataVector weight(this->numData);
      weight.setAll(1.0 / double(this->numData));
      sg::base::OperationEval* opEval = sg::op_factory::createOperationEval(*this->grid);
      // to store certain train data point
      sg::base::DataVector p_train(this->dim);
      // to store certain train data point
      sg::base::DataVector p_test(this->dim);
      sg::base::DataVector tmpweight(this->numData);

      // to store the absolute relative error
      sg::base::DataVector ARE(this->numData);
      // to store the prediction training values
      sg::base::DataVector value_train(this->numData);
      // to store the prediction testing values
      sg::base::DataVector value_test(testData.getNrows());

      sg::base::DataVector beta(this->numBaseLearners);
      sg::base::DataVector logBetaSumR(this->numBaseLearners);
      double errorRate;

      for (size_t count = 0; count < this->numBaseLearners; count++) {
        (this->actualBaseLearners)++;
        std::cout << std::endl;
        std::cout << "This is the " << this->actualBaseLearners << "th weak learner." << std::endl;
        std::cout << std::endl;

        // create coefficient vector
        sg::base::DataVector alpha_train(this->gridPoint);
        sg::base::DataVector alpha_learn(this->gridPoint);
        std::cout << "gridPoint: " << this->gridPoint << std::endl;

        if (this->maxGridPoint->get(count) < this->gridPoint)
          this->maxGridPoint->set(count, (double)this->gridPoint);

        if (!this->refinement) {
          if (count == 0)
            this->sumGridPoint->set(count, (double)gridPoint);
          else
            this->sumGridPoint->set(count, this->sumGridPoint->get(count - 1) + (double)gridPoint);
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
          opEval = sg::op_factory::createOperationEval(*this->grid);
          alpha_learn.resizeZero(alpha_train.getSize());
        }

        //set the alpha for testing data(copy of alpha for training data)
        alpha_learn.copyFrom(alpha_train);

        // calculate the algorithm value of the testing data and training data
        // for training data
        errorRate = 0;
        #pragma omp parallel for schedule(static)

        for (size_t i = 0; i < numData; i++) {
          sg::base::DataVector p_train_private(this->dim);
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
          throw new sg::base::operation_exception("AlgorithmAdaBoostBase::doAdaBoostRT : An unknown power type was specified!");

        sg::base::DataVector TrValueHelper(this->numData);
        double loghelp = log(1 / beta.get(count));

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
        double normalizer = tmpweight.sum();
        tmpweight.mult(1.0 / normalizer);
        weight = tmpweight;

        // for testing data
        #pragma omp parallel for schedule(static)

        for (size_t i = 0; i < testData.getNrows(); i++) {
          sg::base::DataVector p_test_private(this->dim);
          testData.getRow(i, p_test_private);
          value_test.set(i, opEval->eval(alpha_learn, p_test_private));
        }

        sg::base::DataVector TeValueHelper(testData.getNrows());

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
            this->grid = sg::base::Grid::createLinearGrid(this->dim);
            std::cout << std::endl;
            std::cout << "Reset to the regular LinearGrid" << std::endl;
          } else if (this->type == 2) {
            this->grid = sg::base::Grid::createLinearTrapezoidBoundaryGrid(this->dim);
            std::cout << std::endl;
            std::cout << "Reset to the regular LinearTrapezoidBoundaryGrid" << std::endl;
          } else if (this->type == 3) {
            this->grid = sg::base::Grid::createModLinearGrid(this->dim);
            std::cout << std::endl;
            std::cout << "Reset to the regular ModLinearGrid" << std::endl;
          }
          // should not happen because this exception should have been thrown some lines upwards!
          else {
            throw new sg::base::operation_exception("AlgorithmAdaboost : Only 1 or 2 or 3 are supported gridType(1 = Linear Grid, 2 = LinearBoundary Grid, 3 = ModLinear Grid)!");
          }

          sg::base::GridGenerator* gridGen = this->grid->createGridGenerator();
          gridGen->regular(this->level);
          std::cout << std::endl;
          delete gridGen;
        }
      }

      delete opEval;
    }

    void AlgorithmAdaBoostBase::eval(sg::base::DataMatrix& testData, sg::base::DataMatrix& algorithmValueTrain, sg::base::DataMatrix& algorithmValueTest) {
      sg::base::DataMatrix weightsMatrix(this->numData, this->numBaseLearners);
      weightsMatrix.setAll(0.0);

      if (this->boostMode == 1) {
        sg::base::DataVector theHypoWeight(this->numBaseLearners);
        sg::base::DataVector theWeightError(this->numBaseLearners);
        sg::base::DataMatrix decisionMatrix(this->numData, this->numBaseLearners);
        theHypoWeight.setAll(0.0);
        theWeightError.setAll(0.0);
        decisionMatrix.setAll(0.0);
        doDiscreteAdaBoost(theHypoWeight, theWeightError, weightsMatrix, decisionMatrix, testData, algorithmValueTrain, algorithmValueTest);
      } else if (this->boostMode == 2)
        doRealAdaBoost(weightsMatrix, testData, algorithmValueTrain, algorithmValueTest);
      else {
        throw new sg::base::operation_exception("AlgorithmAdaboost : Only 1 or 2 for the boost mode(1 = Discrete Adaboost, 2 = Real Adaboost)!");
      }
    }

    void AlgorithmAdaBoostBase::classif(sg::base::DataMatrix& testData, sg::base::DataVector& algorithmClassTrain, sg::base::DataVector& algorithmClassTest, sg::base::DataMatrix& algorithmValueTrain, sg::base::DataMatrix& algorithmValueTest) {
      eval(testData, algorithmValueTrain, algorithmValueTest);

      for (size_t i = 0; i < this->numData; i++) {
        algorithmClassTrain.set(i, hValue(algorithmValueTrain.get(i, this->actualBaseLearners - 1)));
      }

      for (size_t i = 0; i < testData.getNrows(); i++) {
        algorithmClassTest.set(i, hValue(algorithmValueTest.get(i, this->actualBaseLearners - 1)));
      }
    }

    void AlgorithmAdaBoostBase::getAccuracy(sg::base::DataMatrix& testData, sg::base::DataVector& testDataClass, double* accuracy_train, double* accuracy_test) {
      /* get the accuracy */
      size_t right_test = 0;
      size_t right_train = 0;
      sg::base::DataVector classTrain(this->numData);
      sg::base::DataVector classTest(testDataClass.getSize());
      sg::base::DataMatrix valueTrain(this->numData, this->numBaseLearners);
      sg::base::DataMatrix valueTest(testDataClass.getSize(), this->numBaseLearners);
      classif(testData, classTrain, classTest, valueTrain, valueTest);

      // for training data
      for (size_t i = 0; i < this->numData; i++) {
        if (classTrain.get(i) == this->classes->get(i))
          right_train = right_train + 1;
      }

      *accuracy_train = double(right_train) / double(this->numData);

      // for testing data
      for (size_t i = 0; i < testData.getNrows(); i++) {
        if (classTest.get(i) == testDataClass.get(i))
          right_test = right_test + 1;
      }

      *accuracy_test = double(right_test) / double(classTest.getSize());
    }

    void AlgorithmAdaBoostBase::getROC(sg::base::DataMatrix& validationData, sg::base::DataVector& validationDataClass, double* acc, double* sensitivity, double* specificity, double* precision, double* recall, double* fOneScore) {
      size_t truePos = 0;
      size_t predictPos = 0;
      size_t trueNeg = 0;
      size_t predictNeg = 0;
      size_t actualPos = 0;
      size_t actualNeg = 0;

      sg::base::DataVector classTrain(this->numData);
      sg::base::DataVector classValidation(validationDataClass.getSize());
      sg::base::DataMatrix valueTrain(this->numData, this->numBaseLearners);
      sg::base::DataMatrix valueValidation(validationDataClass.getSize(), this->numBaseLearners);
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

      *acc = double(truePos + trueNeg) / double(classValidation.getSize());
      *sensitivity = double(truePos) / double(actualPos);
      *specificity = double(trueNeg) / double(actualNeg);
      *precision = double(truePos) / double(predictPos);
      *recall = *sensitivity;
      *fOneScore = 2 * (*precision) * (*recall) / ((*precision) + (*recall));
    }

    void AlgorithmAdaBoostBase::getAccuracyBL(sg::base::DataMatrix& testData, sg::base::DataVector& testDataClass, sg::base::DataMatrix& algorithmValueTrain, sg::base::DataMatrix& algorithmValueTest, double* accuracy_train, double* accuracy_test, size_t yourBaseLearner) {
      size_t right_test = 0;
      size_t right_train = 0;

      // for training data
      for (size_t i = 0; i < this->numData; i++) {
        if (hValue(algorithmValueTrain.get(i, yourBaseLearner - 1)) == this->classes->get(i))
          right_train = right_train + 1;
      }

      *accuracy_train = double(right_train) / double(this->numData);

      // for testing data
      for (size_t i = 0; i < testData.getNrows(); i++) {
        if (hValue(algorithmValueTest.get(i, yourBaseLearner - 1)) == testDataClass.get(i))
          right_test = right_test + 1;
      }

      *accuracy_test = double(right_test) / double(testDataClass.getSize());
    }

    void AlgorithmAdaBoostBase::doRefinement(sg::base::DataVector& alpha_ada, sg::base::DataVector& weight_ada, size_t curBaseLearner) {
      bool final_ada = false;

      for (size_t adaptiveStep = 1; adaptiveStep <= this->refineTimes; adaptiveStep++) {

        sg::base::GridGenerator* myGenerator = this->grid->createGridGenerator();
        size_t refineNumber;

        if (this->refineMode == 1) {
          if (this->numOfAda > this->grid->getSize())
            refineNumber = this->grid->getSize();
          else
            refineNumber = this->numOfAda;
        } else if (this->refineMode == 2) {
          refineNumber = (size_t)(this->perOfAda * (double)(this->grid->getSize()));

          //force to refine at least one point
          if (refineNumber == 0)
            refineNumber = 1;
        }
        // should not happen because this exception should have been thrown some lines upwards!
        else {
          throw new sg::base::operation_exception("AlgorithmAdaBoost : Only 1 or 2 are supported refine mode(1 : use grid point number, 2: use grid point percentage)!");
        }

        sg::base::SurplusRefinementFunctor* myRefineFunc = new sg::base::SurplusRefinementFunctor(&alpha_ada, refineNumber, 0.0);
        myGenerator->refine(myRefineFunc);
        delete myRefineFunc;
        delete myGenerator;

        sg::base::GridStorage* gridStorage_ada = this->grid->getStorage();
        size_t gridPts = gridStorage_ada->size();

        std::cout << std::endl;
        std::cout << "Refinement time step: " << adaptiveStep << ", new grid size: " << gridPts << ", refined number of grid points: " << refineNumber << std::endl;

        if (adaptiveStep == this->refineTimes) {
          final_ada = true;

          if (curBaseLearner == 1) {
            this->maxGridPoint->set(curBaseLearner - 1, (double)gridPts);
            this->sumGridPoint->set(curBaseLearner - 1, (double)gridPts);
          } else {
            if (gridPts > this->maxGridPoint->get(curBaseLearner - 2))
              this->maxGridPoint->set(curBaseLearner - 1, (double)gridPts);
            else
              this->maxGridPoint->set(curBaseLearner - 1, this->maxGridPoint->get(curBaseLearner - 2));

            this->sumGridPoint->set(curBaseLearner - 1, this->sumGridPoint->get(curBaseLearner - 2) + (double)gridPts);
          }
        }

        // extend alpha vector (new entries uninitialized)
        alpha_ada.resizeZero(gridPts);

        // calculate new alpha
        alphaSolver(this->lamb, weight_ada, alpha_ada, final_ada);
      }
    }

    double AlgorithmAdaBoostBase::hValue(double realValue) {
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
      size_t mean = (size_t)(this->sumGridPoint->get(baselearner - 1) / (double)baselearner);
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
