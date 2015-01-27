// Copyright (C) 2013 - Michael Baudin
//
// This file must be used under the terms of the 
// GNU Lesser General Public License license
// http://www.gnu.org/copyleft/lesser.html

#ifndef _LOWDISC_SSOBOL_H_
#define _LOWDISC_SSOBOL_H_


//! Class of Scrambled Sobol Sequence
class Ssobol {
public:
	/*
	Startup the Scrambled Sobol sequence.

	Parameters

	INPUTS : 
	dimen : the number of dimensions, DIMEN in {1,2,...,40}
	atmost : the maximum number of elements in the sequence, ATMOST in {1,2,...,2^30-1=1073741823}
	maxd : Maximum Digits of Scrambling Of Owen type Scrambling (suggestion : maxd=30)
	iflag: the scrambling type
	iflag = 0 : No Scrambling
	iflag = 1 : Owen type Scrambling
	iflag = 2 : Faure-Tezuka type Scrambling
	iflag = 3 : Owen + Faure-Tezuka type Scrambling

	OUTPUTS:
	isok = 1 if the parameters are OK.

	Description
	THIS IS MODIFIED ROUTINE OF "INSOBL". 
	FIRST CHECK WHETHER THE USER-SUPPLIED 
	DIMENSION "DIMEN" OF THE QUASI-RANDOM 
	VECTORS IS STRICTLY BETWEEN 1 AND 41. 

	CHECK "ATMOST", AN UPPER BOUND ON THE NUMBER 
	OF CALLS THE USER INTENDS TO MAKE ON "GOSOBL".  IF 
	THIS IS POSITIVE AND LESS THAN 2**30, THEN FLAG(2) = .TRUE. 
	(WE ASSUME WE ARE WORKING ON A COMPUTER WITH 
	WORD LENGTH AT LEAST 31 BITS EXCLUDING SIGN.) 
	THE NUMBER OF COLUMNS OF THE ARRAY V WHICH 
	ARE INITIALIZED IS 
	MAXCOL = NUMBER OF BITS IN ATMOST. 
	IN "GOSOBL" WE CHECK THAT THIS IS NOT EXCEEDED. 

	THE LEADING ELEMENTS OF EACH ROW OF V ARE 
	INITIALIZED USING "VINIT" FROM "BDSOBL". 
	EACH ROW CORRESPONDS TO A PRIMITIVE POLYNOMIAL 
	(AGAIN, SEE "BDSOBL").  IF THE POLYNOMIAL HAS 
	DEGREE M, ELEMENTS AFTER THE FIRST M ARE CALCULATED. 

	THE NUMBERS IN V ARE ACTUALLY BINARY FRACTIONS. 
	LSM ARE LOWER TRIAUGULAR SCRAMBLING MATRICES. 
	USM ARE UPPER TRIAUGULAR SCRMABLING MATRIX. 
	SV ARE SCAMBLING GENERATING MATRICES AND THE NUMBERS 
	ARE BINARY FRACTIONS. 
	"RECIPD" HOLDS 1/(THE COMMON DENOMINATOR OF ALL 
	OF THEM). 


	"INSSOBL" IMPLICITLY COMPUTES THE FIRST SHIFTED 
	VECTOR "LASTQ", AND RETURN IT TO THE CALLING 
	PROGRAM. SUBSEQUENT VECTORS COME FROM "GOSSOBL". 
	"LASTQ" HOLDS NUMERATORS OF THE LAST VECTOR GENERATED. 

	*/
	// This constructor calls seedreset.
	Ssobol(int dimen, int atmost, int iflag, int maxd, int *isok);

	//
	// Sets the seed of the random number generator. 
	// By default, the object always returns the same sequence of 
	// numbers, because the seed is reset at object creation.
	// This method allows to get different scramblings. 
	// newseed : an array of doubles (input), in the interval [0,1].
	// This constructor calls seedset.
	Ssobol(int dimen, int atmost, int iflag, int maxd, double seeds[24], int *isok);

	// Destructor (free the allocated memory)
	~Ssobol();

	// Next element in the Scrambled Sobol Sequence
	//
	// Parameters
	// quasi : an array of doubles (output), quasi[0,1,...,dimen-1]
	void next(double *quasi);

	// dim_num_get -- 
	// gets the spatial dimension for a leaped Halton subsequence.
	int dim_num_get ();

	// gettaus --
	// taus : to determine favorable number of calls
	/* "TAUS" IS FOR DETERMINING "FAVORABLE" VALUES. AS 
	DISCUSSED IN BRATLEY/FOX, THESE HAVE THE FORM 
	N = 2**K WHERE K .GE. (TAUS+S-1) FOR INTEGRATION 
	AND K .GT. TAUS FOR GLOBAL OPTIMIZATION. */
	int gettaus();
private:
	//
	// Fields
	//
	int ssobol_poly[39];
	int ssobol_vinit[40][8];
	double ssobol_recipd;
	int ssobol_lastq[40];
	int ssobol_maxcol;
	int ssobol_count;
	int ssobol_s;
	int ssobol_sv[40][31];
	int ssobol_tau[13];
	unsigned int ssobol_unifseed;

	// Variables for the random number generator.
	int ssobol_seedi;
	int ssobol_seedj;
	double ssobol_seedcarry;
	double ssobol_seedseeds[24];

	/*     THIS FUNCTION CALCULATES THE EXCLUSIVE-OR OF ITS */
	/*     TWO INPUT PARAMETERS */
	int exor(int *iin, int *jin);

	// genscrml --
	/* GENERATING LOWER TRIANULAR SCRAMBLING MATRICES AND SHIFT VECTORS. */
	/* INPUTS : ssobol_s, ssobol_maxcol, maxd */
	/* OUTPUTS : LSM, SHIFT */
	int genscrml(int maxd, int lsm[][31], int *shift);
	
	// genscrmu --
	/* GENERATING UPPER TRIANGULAR SCRAMBLING MATRICES AND */
	/* SHIFT VECTORS. */
	/* INPUTS : ssobol_s, ssobol_maxcol,  */
	/* OUTPUTS : USM, USHIFT */
	int genscrmu(int usm[][31], int *ushift);

	// unirnd --
	// Generates a uniform random number in [0,1]
	/*     Random number generator, adapted from F. James */
	/*     "A Review of Random Number Generators" */
	/*      Comp. Phys. Comm. 60 (1990), pp. 329-344. */
	double unirnd(void);

	// lbitbits --
	// From f2c
	int lbitbits(int a, int b, int len);

	// init --
	// Initialize the current object. 
	// This method is used in the constructor.
	void init(int dimen, int atmost, int iflag, int maxd, int *isok);
public:
	// seedset --
	// Set the seed of the random number generator.
	void seedset(double seeds[24]);

	// setreset --
	// Set the seed of the random number generator to the default seed.
	void seedreset();
};


#endif /* _LOWDISC_SSOBOL_H_ */

