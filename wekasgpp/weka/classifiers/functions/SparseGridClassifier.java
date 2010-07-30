/*****************************************************************************/
/* This file is part of jsgpp, a program package making use of spatially     */
/* adaptive sparse grids to solve numerical problems                         */
/*                                                                           */
/* Copyright (C) 2008 Joerg Blank (blankj@in.tum.de)                         */
/* Copyright (C) 2009 Alexander Heinecke (Alexander.Heinecke@mytum.de)       */
/*                                                                           */
/* jsgpp is free software; you can redistribute it and/or modify             */
/* it under the terms of the GNU General Public License as published by      */
/* the Free Software Foundation; either version 3 of the License, or         */
/* (at your option) any later version.                                       */
/*                                                                           */
/* jsgpp is distributed in the hope that it will be useful,                  */
/* but WITHOUT ANY WARRANTY; without even the implied warranty of            */
/* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the             */
/* GNU Lesser General Public License for more details.                       */
/*                                                                           */
/* You should have received a copy of the GNU General Public License         */
/* along with jsgpp; if not, write to the Free Software                      */
/* Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA */
/* or see <http://www.gnu.org/licenses/>.                                    */
/*****************************************************************************/

package weka.classifiers.functions;

import java.util.Enumeration;
import java.util.Vector;

import weka.de.tum.in.www5.sgpp.*;
import weka.classifiers.Classifier;
import weka.core.Capabilities;
import weka.core.Instance;
import weka.core.Instances;
import weka.core.Option;
import weka.core.OptionHandler;
import weka.core.Utils;
import weka.core.Capabilities.Capability;
import weka.filters.Filter;
//import weka.filters.Normalize;
//import weka.filters.NominalToNumeric;
import weka.filters.unsupervised.attribute.NominalToBinary;
import weka.filters.unsupervised.attribute.ReplaceMissingValues;
import java.util.Random;
import java.util.ArrayList;

/** 
<!-- globalinfo-start -->
 * <p>A Classifier that uses sparse grids to classify instances.<br>
 * </p>
<!-- globalinfo-end -->
 *
   <!-- options-start -->
 * <p>Valid options are: </p>
 * 
 * <pre> -L &lt;level&gt;
 *  grid level.
 *  (Default = 2).</pre>
 * 
 * <pre> -N &lt;numAdaptive&gt;
 *  Using an adaptive Grid with numAdaptive refinements 
 * (Default = 0).</pre>
 * 
 * <pre> -C &lt;cmode&gt;
 *  Value could be 0:laplace 1:laplaceadaptive 2:identity 3:ratio 4:levelsum 5:energy 6:copy and 7:pseudounit.
 *  (Default = 2).</pre>
 * 
 * <pre> -Y &lt;lambda&gt;
 * (Default = 0.1).</pre>
 * 
 * <pre> -I &lt;imax&gt;
 *  max number of iterations
 * (Default = 200).</pre>
 * 
 * <pre> -A &lt;accuracy&gt;
 *  specifies the accuracy of the CG-Iteration
 * (Default = 0.001).</pre>
 *
 * <pre> -M &lt;multiClassMode&gt;
 *  value can be:
 *  0:  learn NumClasses models: model i learns Output==i vs Output!=i
 *  1:  learn just one model\n"  
 * (Default = 0).</pre>
 * 
 *  <pre> -P &lt;polynom&gt;
 *  maximum degree for high order basis functions
 *  Set to 2 or larger to activate. Works only with 'identity' 
 *  (default: 0)</pre>
 *   
 *  <pre> -B
 *  Enables special border base functions
 *  </pre>
 *  
 *  <pre> -R
 *  Use regression approach
 *  </pre>  
 *  
   <!-- options-end -->
 *
 * @author Lennart Johansson
 * @version
 */

public class SparseGridClassifier
extends Classifier
implements OptionHandler
{

	/** for serialization */
	static final long serialVersionUID = 187395900385739571L;

	/**
	 * Call this function to build classifier for the training
	 * data provided.
	 * @param instances -< The training data.
	 * @throws Exception if can't build classification properly.
	 */
	public void buildClassifier(Instances instances) throws Exception
	{
	}	

	/**
	 * This will return a string describing the classifier.
	 * @return The string.
	 */
	public String globalInfo() {
		return "A Classifier that uses sparse grids to classify instances.";
	}

	/**
	 * Returns default capabilities of the classifier.
	 *
	 * @return the capabilities of this classifier
	 */
	public Capabilities getCapabilities() {
		Capabilities result = super.getCapabilities();

		// attributes
//		result.enable(Capability.NOMINAL_ATTRIBUTES);
		result.enable(Capability.NUMERIC_ATTRIBUTES);
//		result.enable(Capability.MISSING_VALUES);

		// class
		result.enable(Capability.NOMINAL_CLASS);
//		result.enable(Capability.MISSING_CLASS_VALUES);

		return result;
	}

	static {
		try {
			System.loadLibrary("jsgpp");
			System.out.println("SparseGridClassifier: jsgpp geladen");
		} catch (UnsatisfiedLinkError e) {
			System.err.println("Native code library failed to load. \n" + e);
			System.exit(1);
		}
	}
}