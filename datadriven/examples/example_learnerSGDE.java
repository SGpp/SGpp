// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

/**
 * \page example_example_learnerSGDE_java Learner SGDE
 * This tutorial demostrates the sparse grid density estimation 
 */

// import all packages
import sgpp.*;

import java.io.File;
import java.io.FileNotFoundException;
import java.util.Scanner;

/**
 * Simple Java SG++ example for the SGDE learner
 */
public class example_learnerSGDE {

    public static DataMatrix readARFF(String fileNameDefined) {
        /**
		* -File class needed to turn stringName to actual file
		*/
        File file = new File(fileNameDefined);
        DataMatrix ans = new DataMatrix(0, 0);
        try{
            /**
			* -read from filePooped with Scanner class
			*/
            Scanner inputStream = new Scanner(file);
            
            for (int i = 0; i < 18; i++) {
                inputStream.next();
            }

            if (inputStream.hasNext()) {
                String[] data = inputStream.next().split(",");
                ans.resize(1, data.length);
                for (int i = 0; i < data.length; i++) {
                    ans.set(0, i, Double.parseDouble(data[i]));
                }
            }
			
			/**
			* hashNext() loops line-by-line
			*/
            while(inputStream.hasNext()){
                /**
				* read single line, put in string
				*/
                String[] data = inputStream.next().split(",");
                long row = ans.appendRow();
                for (int i = 0; i < data.length; i++) {
                    ans.set(row, i, Double.parseDouble(data[i]));
                }
            }
            /**
			* after loop, close scanner
			*/
            inputStream.close();
        }catch (FileNotFoundException e){
            e.printStackTrace();
        }
        return ans;
    }

   /**
   * Main method.
   *
   * @param args ignored
   */
  public static void main(String[] args) {
    // Two working possiblities for loading the shared object file:
    //java.lang.System.load("/path/to/SGpp/trunk/lib/jsgpp/libjsgpp.so");
    sgpp.LoadJSGPPLib.loadJSGPPLib();

    String filename = "../datasets/friedman_4d_2000.arff";

    System.out.println("# loading file: " + filename);
    DataMatrix samples = readARFF(filename);

    /**
	* Configure grid
	*/
    System.out.println("# create grid config");
    RegularGridConfiguration gridConfig = new RegularGridConfiguration();
    gridConfig.setDim_(samples.getNcols());
    gridConfig.setLevel_(4);
    gridConfig.setType_(GridType.Linear);

    /** 
	* Configure adaptive refinement
	*/
    System.out.println("# create adaptive refinement config");
    AdaptivityConfiguration adaptConfig = new AdaptivityConfiguration();
    adaptConfig.setNoPoints_(5);
    adaptConfig.setNumRefinements_(0);

    /** 
	* Configure solver
	*/
    System.out.println("# create solver config");
    SLESolverConfiguration solverConfig = new SLESolverConfiguration();
    solverConfig.setType_(SLESolverType.CG);
    solverConfig.setMaxIterations_(100);
    solverConfig.setEps_(1e-10);
    solverConfig.setThreshold_(1e-10);

    /** 
	* Configure regularization
	*/
    System.out.println("# create regularization config");
    RegularizationConfiguration regularizationConfig = new RegularizationConfiguration();
    regularizationConfig.setRegType_(RegularizationType.Laplace);

    /**
	* Configure learner
	*/
    System.out.println("# create learner config");
    LearnerSGDEConfiguration learnerConfig = new LearnerSGDEConfiguration();
    learnerConfig.setDoCrossValidation_(false);
    learnerConfig.setKfold_(5);
    learnerConfig.setLambdaStart_(1e-1);
    learnerConfig.setLambdaEnd_(1e-10);
    learnerConfig.setLambdaSteps_(5);
    learnerConfig.setLogScale_(true);
    learnerConfig.setShuffle_(true);
    learnerConfig.setSeed_(1234567);
    learnerConfig.setSilent_(false);

	/**
	* Initialize and run learner
	*/
    System.out.println("# creating the learner");
    LearnerSGDE learner = new LearnerSGDE(gridConfig, adaptConfig, solverConfig, regularizationConfig, learnerConfig);
    learner.initialize(samples);

	/**
	* For comparison, run the sparse grid kernel-based learner
	*/
    KernelDensityEstimator kde = new KernelDensityEstimator(samples);
    DataVector x = new DataVector(3);

	/**
	* View resulting pdf in the middle of the domain
	*/
    for (int i = 0; i < x.getSize(); i++) {
        x.set(i, 0.5);
    }

    System.out.println("--------------------------------------------------------");
    System.out.println(learner.getGrid().getSize() + " -> " + learner.getAlpha().sum());
    System.out.println("pdf_SGDE(x) = " + learner.pdf(x) + " ~ " + kde.pdf(x) + " = pdf_KDE(x)");
    System.out.println("mean_SGDE(x) = " + learner.mean() + " ~ " + kde.mean() + " = mean_KDE(x)");
    System.out.println("var_SGDE(x) = " + learner.variance() + " ~ " + kde.variance() + " = var_KDE(x)");
  }
}
