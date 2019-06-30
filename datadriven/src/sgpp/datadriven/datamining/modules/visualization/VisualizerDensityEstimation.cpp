/*
 * Copyright (C) 2008-today The SG++ project
 * This file is part of the SG++ project. For conditions of distribution and
 * use, please see the copyright notice provided with SG++ or at
 * sgpp.sparsegrids.org
 *
 * VisualizerDensityEstimation.hpp
 *
 *  Created on: 16th Jun 2019
 *      Author: Vincent Bautista
 */

#include <sgpp/datadriven/datamining/modules/visualization/VisualizerDensityEstimation.hpp>
#include <sgpp/datadriven/datamining/modules/visualization/algorithms/bhtsne/tsne.h>
#include <iostream>
#include <vector>


using sgpp::base::DataMatrix;
using sgpp::base::DataVector;
using sgpp::datadriven::TSNE;


namespace sgpp{
namespace datadriven{


 VisualizerDensityEstimation::VisualizerDensityEstimation(VisualizerConfiguration config){

  this->config = config;
 }

 void VisualizerDensityEstimation::visualize(ModelFittingBase &model)
 {
  createOutputDirectory();
  getLinearCuts(model);
  getHeatmap(model);
  runTsne(model);

  if( config.getGeneralConfig().targetFileType==VisualizationFileType::CSV)
  {
   storeGrid(model);
  }
 }

 void VisualizerDensityEstimation::storeGrid(ModelFittingBase &model)
 {

  ModelFittingBaseSingleGrid* gridModel = dynamic_cast<ModelFittingBaseSingleGrid*>(&model);

  auto grid = gridModel->getGrid().clone();

  DataMatrix gridMatrix;

  grid->getStorage().getCoordinateArrays(gridMatrix);


  CSVTools::writeMatrixToCSVFile(config.getGeneralConfig().targetFile+"/grid", gridMatrix);
 }
 void VisualizerDensityEstimation::runTsne(ModelFittingBase &model){

  DataMatrix data = model.getDataset()->getData();

  double* input = model.getDataset()->getData().data();

  DataVector evaluation(data.getNrows());

  model.evaluate(data, evaluation);


  int N = model.getDataset()->getNumberInstances();
  int D = model.getDataset()->getDimension();

  double* output;

  if(D>=config.getVisualizationParameters().targetDimension)
  {
   std::cout << "Compresiing with tsne to "<<std::to_string(config.getVisualizationParameters().targetDimension)
    <<" dimensions"<<std::endl;
  TSNE::run(input, N, D , output, config.getVisualizationParameters().targetDimension,
    config.getVisualizationParameters().perplexity, config.getVisualizationParameters().theta,
    config.getVisualizationParameters().seed, false, config.getVisualizationParameters().maxNunmberIterations);
   D = config.getVisualizationParameters().targetDimension;
  }
  else
  {
   output = input;
  }

  DataMatrix compressedModel(output, N, D);

  compressedModel.appendCol(evaluation);

  CSVTools::writeMatrixToCSVFile(config.getGeneralConfig().targetFile+"/tsne_compression", compressedModel);


 }

 void VisualizerDensityEstimation::getLinearCuts(ModelFittingBase &model)
 {
  std::cout <<"Generating the linear cuts"<<std::endl;
  std::string command("mkdir ");
  command.append(config.getGeneralConfig().targetFile);
  command.append("/LinearCuts");
  system(command.data());

  auto nDimensions  = model.getFitterConfiguration().getGridConfig().dim_;

  DataMatrix cutMatrix(0,nDimensions);

  if (nDimensions >=3)
  {
   getLinearCutsMore3D(cutMatrix, model);
  }
  else if(nDimensions==2)
  {
   getLinearCuts2D(cutMatrix, model);
  }
  else{
   getLinearCuts1D(cutMatrix, model);
  }

 }

 void VisualizerDensityEstimation::getHeatmap(ModelFittingBase &model)
{

  std::cout <<"Generating the heatmaps "<<std::endl;
  std::string command("mkdir ");
  command.append(config.getGeneralConfig().targetFile);
  command.append("/Heatmaps");
  system(command.data());
  auto nDimensions  = model.getFitterConfiguration().getGridConfig().dim_;

  if(nDimensions==1)
  {
   std::cout << "Heatmap generation is not available for models of 1 dimension" <<std::endl;
   return;
  }

  DataMatrix heatMapMatrix(0,nDimensions);

  if (nDimensions >=4)
  {
   getHeatmapMore4D(heatMapMatrix, model);
  }
  else if(nDimensions==3){
   getHeatmap3D(heatMapMatrix, model);
  }
  else
  {
   getHeatmap2D(heatMapMatrix, model);
  }



}
 void VisualizerDensityEstimation::translateColumns(DataMatrix &matrix, size_t maxColumns)
  {

   DataMatrix temp(matrix);

   DataVector column(matrix.getNrows());

   for (size_t dimension=0; dimension<maxColumns-1;dimension++)
   {

     matrix.getColumn(dimension,column);

     temp.setColumn(dimension+1,column);
   }

   matrix.getColumn(maxColumns-1,column);

   temp.setColumn(0,column);

   matrix.copyFrom(temp);

  }

 void VisualizerDensityEstimation::translateColumnsRight(DataMatrix &matrix, std::vector <size_t> indexes)
  {

   DataMatrix temp(matrix);

   DataVector column(matrix.getNrows());

   for (size_t dimension=0; dimension<indexes.size()-1;dimension++)
   {

     matrix.getColumn(indexes.at(dimension),column);

     temp.setColumn(indexes.at(dimension+1),column);
   }

   matrix.getColumn(indexes.back(),column);

   temp.setColumn(indexes.front(),column);

   matrix.copyFrom(temp);

  }

 void VisualizerDensityEstimation::translateColumnsLeft(DataMatrix &matrix, std::vector <size_t> indexes)
  {

   DataMatrix temp(matrix);

   DataVector column(matrix.getNrows());

   for (size_t dimension=indexes.size()-1; dimension>0;dimension--)
   {

     matrix.getColumn(indexes.at(dimension),column);

     temp.setColumn(indexes.at(dimension-1),column);
   }

   matrix.getColumn(indexes.front(),column);

   temp.setColumn(indexes.back(),column);

   matrix.copyFrom(temp);

  }

 void VisualizerDensityEstimation::updateIndexesCuts(std::vector <size_t> &indexes, DataMatrix &matrix)
 {
   if(indexes.at(2)<matrix.getNcols()-1)
   {
     indexes.at(2)++;
     swapColumns(matrix, indexes.at(2)-1, indexes.at(2));
   }
   else
   {
    if(indexes.at(1)<matrix.getNcols()-2)
    {
     indexes.at(1)++;
     indexes.at(2)=indexes.at(1)+1;

     swapColumns(matrix, indexes.at(2)-1, indexes.at(2));
     swapColumns(matrix, indexes.at(1)-1, indexes.at(1));

    }
    else
    {
     indexes.at(0)++;
     indexes.at(1)=indexes.at(0)+1;
     indexes.at(2) = indexes.at(1)+1;

     swapColumns(matrix, matrix.getNcols()-1, indexes.at(2));
     swapColumns(matrix, matrix.getNcols()-2, indexes.at(1));
     swapColumns(matrix, indexes.at(0)-1, indexes.at(0));
    }
   }
 }

 void VisualizerDensityEstimation::updateIndexesHeatmap(std::vector <size_t> &indexes, DataMatrix &matrix)
  {
    if(indexes.at(3)<matrix.getNcols()-1)
    {
      indexes.at(3)++;
      swapColumns(matrix, indexes.at(3)-1, indexes.at(3));
    }
    else
    {
     if(indexes.at(2)<matrix.getNcols()-2)
     {
      indexes.at(2)++;
      indexes.at(3)=indexes.at(2)+1;

      swapColumns(matrix, indexes.at(2)-1, indexes.at(2));
      swapColumns(matrix, matrix.getNcols()-1, indexes.at(3));


     }
     else
     {
      if(indexes.at(1)<matrix.getNcols()-3)
       {
        indexes.at(1)++;
        indexes.at(2)=indexes.at(1)+1;
        indexes.at(3)=indexes.at(2)+1;

        swapColumns(matrix, matrix.getNcols()-1, indexes.at(3));
        swapColumns(matrix, matrix.getNcols()-2, indexes.at(2));
        swapColumns(matrix, indexes.at(1)-1, indexes.at(1));

       }
      else
       {
       indexes.at(0)++;
       indexes.at(1)=indexes.at(0)+1;
       indexes.at(2) = indexes.at(1)+1;
       indexes.at(2)=indexes.at(1)+1;
       indexes.at(3)=indexes.at(2)+1;

       swapColumns(matrix, matrix.getNcols()-1, indexes.at(3));
       swapColumns(matrix, matrix.getNcols()-2, indexes.at(2));
       swapColumns(matrix, matrix.getNcols()-3, indexes.at(1));
       swapColumns(matrix, indexes.at(0)-1, indexes.at(0));
       }
     }
    }
  }

 void VisualizerDensityEstimation::swapColumns(DataMatrix &matrix, size_t col1, size_t col2){

     DataVector temp1(matrix.getNrows());

     DataVector temp2(matrix.getNrows());

     matrix.getColumn(col1,temp1);

     matrix.getColumn(col2, temp2);

     matrix.setColumn(col2, temp1);

     matrix.setColumn(col1, temp2);

 }
 void VisualizerDensityEstimation::getLinearCutsMore3D(DataMatrix &cutMatrix, ModelFittingBase &model)
 {

  std::string outputDir(config.getGeneralConfig().targetFile+"/LinearCuts/");

  for(double dim1=0; dim1<=1; dim1+=0.25)
   {
     for(double dim2=0; dim2<=1; dim2+=0.25)
     {
      for(double dim3=0; dim3<=1 ;dim3+=0.01)
          {
               DataVector row(cutMatrix.getNcols(),0.5);
               row.set(2,dim2);
               row.set(0,dim3);
               row.set(1, dim1);

               cutMatrix.appendRow(row);
          }

     }
   }
   std::vector <size_t> variableColumnIndexes = {0,1,2};

   while(variableColumnIndexes.at(2)<cutMatrix.getNcols())
   {

    std::string subfolder(outputDir+"dimensions_"+
      std::to_string(variableColumnIndexes.at(0)+1)+"_"+
      std::to_string(variableColumnIndexes.at(1)+1)+"_"+
      std::to_string(variableColumnIndexes.at(2)+1));

    std::string command("mkdir "+subfolder);

    system(command.data());


    for(size_t combination=0;combination<3;combination++)
      {
          DataMatrix cutResults(cutMatrix);
          DataVector evaluation(cutMatrix.getNrows());

          model.evaluate(cutMatrix,evaluation);

          cutResults.appendCol(evaluation);

          translateColumnsRight(cutMatrix, variableColumnIndexes);

          CSVTools::writeMatrixToCSVFile(subfolder+"/Cut_var_dimension_"+
            std::to_string(variableColumnIndexes.at(combination)+1),
                cutResults);

      }


    updateIndexesCuts(variableColumnIndexes, cutMatrix);

   }

 }

 void VisualizerDensityEstimation::getLinearCuts2D(DataMatrix &cutMatrix, ModelFittingBase &model)
 {

  std::string outputDir(config.getGeneralConfig().targetFile+"/LinearCuts/");

  for(double dim1=0; dim1<=1; dim1+=0.25)
  {
    for(double dim2=0; dim2<=1; dim2+=0.01)
    {
      DataVector row(2,0.5);
      row.set(0,dim2);
      row.set(1, dim1);
      cutMatrix.appendRow(row);
     }
  }

  for(size_t combination=0;combination<2;combination++)
  {
      DataMatrix cutResults(cutMatrix);
      DataVector evaluation(cutMatrix.getNrows());

      model.evaluate(cutMatrix,evaluation);

      cutResults.appendCol(evaluation);

      translateColumns(cutMatrix, cutMatrix.getNcols());

      CSVTools::writeMatrixToCSVFile(outputDir+"Cut_dimensions_1_2_variable_dimension"+
        std::to_string(combination+1),
            cutResults);
  }

 }

 void VisualizerDensityEstimation::getLinearCuts1D(DataMatrix &cutMatrix, ModelFittingBase &model)
 {
  std::string outputDir(config.getGeneralConfig().targetFile+"/LinearCuts/");
  for(double dim1=0; dim1<=1; dim1+=0.01)
   {
    DataVector row(1,0.5);
    row.set(0, dim1);
    cutMatrix.appendRow(row);
   }

    DataMatrix cutResults(cutMatrix);
    DataVector evaluation(cutMatrix.getNrows());

    model.evaluate(cutMatrix,evaluation);

    cutResults.appendCol(evaluation);

    CSVTools::writeMatrixToCSVFile(outputDir+"Cut_dimension 1",cutResults);
 }

 void VisualizerDensityEstimation::getHeatmapMore4D(DataMatrix &heatMapMatrix, ModelFittingBase &model)
 {
  std::string outputDir(config.getGeneralConfig().targetFile+"/Heatmaps/");
  for(double dim1=0; dim1<=1; dim1+=0.25)
   {
     for(double dim2=0; dim2<=1; dim2+=0.25)
     {
      for(double dim3=0; dim3<=1 ;dim3+=0.05)
          {
       for(double dim4=0; dim4<=1 ;dim4+=0.05)
           {
               DataVector row(heatMapMatrix.getNcols(),0.5);
               row.set(0,dim3);
               row.set(1,dim4);
               row.set(2, dim1);
               row.set(3, dim2);

               heatMapMatrix.appendRow(row);
          }
       }

     }
   }

   std::vector <size_t> variableColumnIndexes = {0,1,2,3};

   std::vector<size_t> workingIndexes={0,0,0};

   while(variableColumnIndexes.at(3)<heatMapMatrix.getNcols())
   {
    std::string subfolder(outputDir+"dimensions_"+
    std::to_string(variableColumnIndexes.at(0)+1)+"_"+
    std::to_string(variableColumnIndexes.at(1)+1)+"_"+
    std::to_string(variableColumnIndexes.at(2)+1)+"_"+
    std::to_string(variableColumnIndexes.at(3)+1));

    std::string command("mkdir "+subfolder);

    system(command.data());


    std::copy(variableColumnIndexes.begin()+1, variableColumnIndexes.end(),
      workingIndexes.begin());

    for(size_t iteration=0; iteration<2;iteration++)
    {

     for(size_t combination=0;combination<3;combination++)
       {
           DataMatrix heatMapResults(heatMapMatrix);
           DataVector evaluation(heatMapMatrix.getNrows());

           model.evaluate(heatMapMatrix,evaluation);

           heatMapResults.appendCol(evaluation);

           if(iteration==0)
           {
           CSVTools::writeMatrixToCSVFile(subfolder+"/Heatmap_var_dimensions_"
                      +std::to_string(variableColumnIndexes.at(0)+1)+"_"+
                      std::to_string(variableColumnIndexes.at(combination+1)+1), heatMapResults);
           translateColumnsRight(heatMapMatrix, workingIndexes);

           }
           else
           {

            CSVTools::writeMatrixToCSVFile(subfolder+"/Heatmap_var_dimensions_"
                              +std::to_string(variableColumnIndexes.at((combination<2)?1:2)+1)+"_"+
                              std::to_string(variableColumnIndexes.at((combination<1)?2:3)+1), heatMapResults);
            translateColumnsLeft(heatMapMatrix, workingIndexes);
           }

       }

    if(iteration==0)
    {
     translateColumnsRight(heatMapMatrix, variableColumnIndexes);
    }

    }
    translateColumnsLeft(heatMapMatrix, variableColumnIndexes);

    updateIndexesHeatmap(variableColumnIndexes, heatMapMatrix);
   }
 }


 void VisualizerDensityEstimation::getHeatmap3D(DataMatrix &heatMapMatrix, ModelFittingBase &model)
  {
  std::string outputDir(config.getGeneralConfig().targetFile+"/Heatmaps/");
   for(double dim1=0; dim1<=1; dim1+=0.25)
    {
     for(double dim2=0; dim2<=1 ;dim2+=0.05)
         {
      for(double dim3=0; dim3<=1 ;dim3+=0.05)
          {
              DataVector row(3,0.5);
              row.set(0,dim3);
              row.set(1,dim2);
              row.set(2, dim1);
              heatMapMatrix.appendRow(row);
         }
      }

    }

    for(size_t combination=1;combination<=3;combination++)
    {
        DataMatrix heatMapResults(heatMapMatrix);
        DataVector evaluation(heatMapMatrix.getNrows());

        model.evaluate(heatMapMatrix,evaluation);

        heatMapResults.appendCol(evaluation);

        translateColumns(heatMapMatrix, heatMapMatrix.getNcols());

        CSVTools::writeMatrixToCSVFile(outputDir+"dimensions_1_2_3_plotNumber_"+
          std::to_string(combination),heatMapResults);
    }
  }

 void VisualizerDensityEstimation::getHeatmap2D(DataMatrix &heatMapMatrix, ModelFittingBase &model)
 {
  std::string outputDir(config.getGeneralConfig().targetFile+"/Heatmaps/");
    for(double dim1=0; dim1<=1 ;dim1+=0.05)
    {
     for(double dim2=0; dim2<=1 ;dim2+=0.05)
         {
             DataVector row(2,0.5);
             row.set(0,dim1);
             row.set(1,dim2);
             heatMapMatrix.appendRow(row);
        }
     }

    DataMatrix heatMapResults(heatMapMatrix);
    DataVector evaluation(heatMapMatrix.getNrows());

    model.evaluate(heatMapMatrix,evaluation);

    heatMapResults.appendCol(evaluation);

    CSVTools::writeMatrixToCSVFile(outputDir+"dimensions_1_2",heatMapResults);

 }


} //namespace datadriven
}//namespace sgpp
