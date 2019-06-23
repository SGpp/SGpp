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
#include <iostream>
#include <vector>

using sgpp::base::DataMatrix;
using sgpp::base::DataVector;


namespace sgpp{
namespace datadriven{


 VisualizerDensityEstimation::VisualizerDensityEstimation(VisualizerConfiguration config){

  this->config = config;
 }

 void VisualizerDensityEstimation::visualize(ModelFittingBase &model)
 {
  getLinearCuts(model);
  getHeatmap(model);
  run_tsne();
 }

 void VisualizerDensityEstimation::run_tsne(){
  std::cout <<"Code to run tsne"<<std::endl;
 }

 void VisualizerDensityEstimation::getLinearCuts(ModelFittingBase &model)
 {
  std::cout <<"Generating the linear cuts"<<std::endl;
  auto nDimensions  = model.getFitterConfiguration().getGridConfig().dim_;

  DataMatrix cutMatrix(0,nDimensions);

  if (nDimensions >=3)
  {
   for(double dim1=0; dim1<=1; dim1+=0.25)
    {
      for(double dim2=0; dim2<=1; dim2+=0.25)
      {
       for(double dim3=0; dim3<=1 ;dim3+=0.01)
           {
                DataVector row(nDimensions,0.5);
                row.set(2,dim2);
                row.set(0,dim3);
                row.set(1, dim1);

                cutMatrix.appendRow(row);
           }

      }
    }
    std::vector <size_t> variableColumnIndexes = {0,1,2};
    while(variableColumnIndexes.at(2)<nDimensions)
    {
     for(size_t combination=0;combination<3;combination++)
       {
           DataMatrix cutResults(cutMatrix);
           DataVector evaluation(cutMatrix.getNrows());

           model.evaluate(cutMatrix,evaluation);

           cutResults.appendCol(evaluation);

           translateColumns(cutMatrix, variableColumnIndexes);

           CSVTools::writeMatrixToCSVFile("Cut_variable dimensions"+
              std::to_string(variableColumnIndexes.at(0))+"_"+
              std::to_string(variableColumnIndexes.at(1))+"_"+
              std::to_string(variableColumnIndexes.at(2))+"_plotNumber"
              +std::to_string(variableColumnIndexes.at(combination)),
                 cutResults);
       }

     if(variableColumnIndexes.at(1)==1)
     {
      translateColumns(cutMatrix, cutMatrix.getNcols());
     }
     else
     {
      translateColumns(cutMatrix, cutMatrix.getNcols()-1);
     }
     updateIndexes(variableColumnIndexes);


    }

  }


 }

 void VisualizerDensityEstimation::getHeatmap(ModelFittingBase &model)
{

  std::cout <<"Generating the heatmaps "<<std::endl;
  auto nDimensions  = model.getFitterConfiguration().getGridConfig().dim_;

  DataMatrix heatMapMatrix(0,nDimensions);

  if (nDimensions >=4)
  {
   for(double dim1=0; dim1<=1; dim1+=0.25)
    {
      for(double dim2=0; dim2<=1; dim2+=0.25)
      {
       for(double dim3=0; dim3<=1 ;dim3+=0.05)
           {
        for(double dim4=0; dim4<=1 ;dim4+=0.05)
            {
                DataVector row(nDimensions,0.5);
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

    while(variableColumnIndexes.at(3)<nDimensions)
    {
     std::vector<size_t> workingIndexes={1,2,3};
     for(size_t iteration=1; iteration<=2;iteration++)
     {

      for(size_t combination=1;combination<=3;combination++)
        {
            DataMatrix heatMapResults(heatMapMatrix);
            DataVector evaluation(heatMapMatrix.getNrows());

            model.evaluate(heatMapMatrix,evaluation);

            heatMapResults.appendCol(evaluation);

            translateColumns(heatMapMatrix, workingIndexes);

            CSVTools::writeMatrixToCSVFile("Heatmap_dimensions"+
               std::to_string(variableColumnIndexes.at(0))+"_"+
               std::to_string(variableColumnIndexes.at(1))+"_"+
               std::to_string(variableColumnIndexes.at(2))+"_"+
               std::to_string(variableColumnIndexes.at(3))+"_plotNumber"+
               std::to_string(iteration)+"_"+std::to_string(combination), heatMapResults);
        }
      workingIndexes.at(0)=0;
     }
     if(variableColumnIndexes.at(1)==1)
     {
      translateColumns(heatMapMatrix, heatMapMatrix.getNcols());
     }
     else
     {
      translateColumns(heatMapMatrix, heatMapMatrix.getNcols()-1);
     }
     updateIndexes(variableColumnIndexes);


    }

  }



}
 void VisualizerDensityEstimation::translateColumns(DataMatrix &matrix, size_t maxColumns)
  {
   std::cout << "Translating columns" << std::endl;
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

 void VisualizerDensityEstimation::translateColumns(DataMatrix &matrix, std::vector <size_t> indexes)
  {
   std::cout << "Translating columns" << std::endl;
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

 void VisualizerDensityEstimation::updateIndexes(std::vector <size_t> &indexes)
 {
  std::cout <<indexes.size()-2<<"_"<<indexes.at(indexes.size()-2)<<std::endl;

   if(indexes.at(indexes.size()-2)==indexes.size()-2)
   {

    for(size_t index=0; index<=indexes.size()-1;index++)
    {
     indexes.at(index)++;
    }

   }
   else
   {
    for(size_t index=0; index<=indexes.size()-2;index++)
    {
     if(indexes.at(index)==index+1)
     {
      indexes.at(index)=index;
      break;
     }
    }
   }
 }
} //namespace datadriven
}//namespace sgpp
