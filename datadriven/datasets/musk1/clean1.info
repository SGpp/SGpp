1. Title: MUSK "Clean1" database

2. Sources:
   (a) Creators:  AI Group at Arris Pharmaceutical Corporation
        contact:  David Chapman or Ajay Jain
                  Arris Pharmaceutical Corporation
                  385 Oyster Point Blvd.
                  South San Francisco, CA 94080
                  415-737-8600
                  zvona@arris.com, jain@arris.com
   (b) Donor:     Tom Dietterich
                  Department of Computer Science
                  Oregon State University
                  Corvallis, OR 97331
                  503-737-5559
                  tgd@cs.orst.edu
   (c) Date received: September 12, 1994

3. Past Usage:
   Dietterich, T. G., Lathrop, R. H., Lozano-Perez, T. (submitted)
   Solving the multiple-instance problem with axis-parallel rectangles.
   Submitted to Artificial Intelligence.

   This paper compares several axis-parallel rectangle algorithms and
   includes the following table:

        Algorithm                TP FN FP TN  errs  %correct [CI]
        iterated-discrim APR     42  5  2 43    7     92.4 [87.0--97.8]
        GFS elim-kde APR         46  1  7 38    8     91.3 [85.5--97.1]
        GFS elim-count APR       46  1  8 37    9     90.2 [84.2--96.3]
        GFS all-positive APR     47  0 15 30   15     83.7 [76.2--91.2]
        all-positive APR         36 11  7 38   18     80.4 [72.3--88.5]
        backpropagation          45  2 21 24   23     75.0 [66.2--83.8]
        C4.5 (pruned)            42  5 24 21   29     68.5 [40.9--61.3]

        key: TP = true positives
             FN = false negatives
             FP = false positives
             TN = true negatives
             errs = errors = FN+FP
             %correct = 10-fold cross-validation %correct.
             CI = 95% confidence interval on proportion of correct
             predictions.
             For explanations of the various algorithms, see the
             paper. 

        C4.5 and backprop were applied ignoring the multiple instance
        problem (see below) during training, but obeying it during
        testing.  

        This paper also gives more details on the construction of the
        data set. 

        This paper also describes an artificial generator that can
        generate data sets with statistics and properties similar to
        this one.

4. Relevant Information:
   This dataset describes a set of 92 molecules of which 47 are judged
   by human experts to be musks and the remaining 45 molecules are
   judged to be non-musks.  The goal is to learn to predict whether
   new molecules will be musks or non-musks.  However, the 166 features
   that describe these molecules depend upon the exact shape, or
   conformation, of the molecule.  Because bonds can rotate, a single
   molecule can adopt many different shapes.  To generate this data
   set, the low-energy conformations of the molecules were generated
   and then filtered to remove highly similar conformations. This left
   476 conformations.  Then, a feature vector was extracted that
   describes each conformation.

   This many-to-one relationship between feature vectors and molecules
   is called the "multiple instance problem".  When learning a
   classifier for this data, the classifier should classify a molecule
   as "musk" if ANY of its conformations is classified as a musk.  A
   molecule should be classified as "non-musk" if NONE of its
   conformations is classified as a musk.

5. Number of Instances  476

6. Number of Attributes 168 plus the class.

7. For Each Attribute:
   
   Attribute:           Description:
   molecule_name:       Symbolic name of each molecule.  Musks have names such
                        as MUSK-188.  Non-musks have names such as
                        NON-MUSK-jp13.
   conformation_name:   Symbolic name of each conformation.  These
                        have the format MOL_ISO+CONF, where MOL is the
                        molecule number, ISO is the stereoisomer
                        number (usually 1), and CONF is the
                        conformation number. 
   f1 through f162:     These are "distance features" along rays (see
                        paper cited above).  The distances are
                        measured in hundredths of Angstroms.  The
                        distances may be negative or positive, since
                        they are actually measured relative to an
                        origin placed along each ray.  The origin was
                        defined by a "consensus musk" surface that is
                        no longer used.  Hence, any experiments with
                        the data should treat these feature values as
                        lying on an arbitrary continuous scale.  In
                        particular, the algorithm should not make any
                        use of the zero point or the sign of each
                        feature value. 
   f163:                This is the distance of the oxygen atom in the
                        molecule to a designated point in 3-space.
                        This is also called OXY-DIS.
   f164:                OXY-X: X-displacement from the designated
                        point.
   f165:                OXY-Y: Y-displacement from the designated
                        point.
   f166:                OXY-Z: Z-displacement from the designated
                        point. 
   class:               0 => non-musk, 1 => musk

   Please note that the molecule_name and conformation_name attributes
   should not be used to predict the class.

8. Missing Attribute Values: none.

9. Class Distribution: 
   Musks:     47
   Non-musks: 45
