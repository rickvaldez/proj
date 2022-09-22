
#include <Python.h>
#include <math.h>

//// to compile in ubuntu 12.04 use gcc -shared -fPIC -I/usr/include/python2.7/ -lpython2.7 -o myModule.so myModule.c





static PyObject* py_SimpsonRule(PyObject* self,PyObject* args)
{
	//added a capital S at the end to tell me its from Simpsons rule function
	// just so i dont accidentally edit another function with similar definitions
	double lenghtHS; //steplenght 
	double SumS = 0; //sumofterms
	double linecounterS; // counts the lines of function file to make array 
	int numlinesS = 0; //states the number of lines in function file
	int countAS = 0; //counter for loops with arrays 
//	double *weightAS = 0; // declaring pointers, will be used with numlines to make arrays
//	double *functS = 0;
//	double *xS = 0;


	//input file funtions
	//location varies
	
	FILE *fnctS;
	fnctS = fopen("functionfile.txt","r");
	
	//for array formating, counts lines

	while(fscanf(fnctS, "%lf %lf", &linecounterS,&linecounterS) != EOF)
	{
		numlinesS++;
	}
	
	fclose(fnctS);


	//opens again but to read the file this time and make arrays

	fnctS = fopen("functionfile.txt","r");
	
    /*
     * This creates arrays of size linecounterS++
     */	
	
	//array for the function file values
	
	double *weightAS = (double *) malloc (numlinesS * sizeof(double) ); 
	double *functS = (double *) malloc (numlinesS * sizeof(double) );
	double *xS = (double *) malloc (numlinesS * sizeof(double) );



//used to go through the file and retreive values and store then to arrays
	for(countAS = 0;countAS < numlinesS ; countAS++)
	{
		fscanf(fnctS, "%lf %lf", &xS[countAS],&functS[countAS]);
	}
		// used for stepsize, to calculate stepsize of file
	lenghtHS = xS[numlinesS -1 ] - xS[numlinesS - 2 ];

    /*
     * This loop creates the array for weighting the sum when doing the interval.
     * The equation is:
     * b
     * ? f(x) dx =  h * f1 + 4h * f2 + 2h * f3 + 4h *f3 + ... + h * fN
     * a            3        3         3         3              3
     * h is the interval length, N is the number of terms,
     * f1 = f(a) and fN = f(b)
     *
     * The weights are defined as the term that is not the function.
     * h    4h  2h  4h ... 2h  4h  h
     * 2    3   3   3      3   3   3
     *
     *
     */	
	
	
	//array for weighing sum
	if(numlinesS % 2==1) //modulo opperator 
	{
		for(countAS = 0; countAS < numlinesS; countAS++)
		{
			if(countAS == 0 || countAS == (numlinesS - 1))
			{
				weightAS[countAS] = lenghtHS/3;
			}
			else if ((countAS % 2) == 1)
			{
				weightAS[countAS] = 4*lenghtHS/3;
			}
			else
			{
				weightAS[countAS] = 2*lenghtHS/3;
			}
		}
	}
	//sums terms
	//needs to be even number of points as well
    /*
     * The equation is:
     * b
     * ? f(x) dx =  h * f1 + 4h * f2 + 2h * f3 + 4h *f3 + ... + h * fN
     * a            3         3         3         3             3
     * h is the interval length, N is the number of terms,
     * f1 = f(a) and fN = f(b)
     *
     * This loop sums the terms of the series, using the weightArray and
     * functionValue array appropriately.
     */

	if(numlinesS % 2 == 1 )
	{
		for(countAS = 0; countAS < numlinesS ; countAS++)
		{
			SumS += weightAS[countAS]*functS[countAS];
 	       /*
    	     * This adds the point where the 3.8ths rule and 1/3 rule combine twice.
   	      */
			if(countAS == 3)	
			{	
			
				SumS += lenghtHS/3 * functS[countAS];
			
			}

	}
	printf("Simpson's approx :  %f\n",SumS);
	

	PyArg_ParseTuple(args, "dd", &SumS);
	return Py_BuildValue("d", SumS);
	
	} else 
		{
		printf("error even # of points. Simpsons not applicable \n");
		char *s = "error even # of points. Simpsons not applicable \n";
	    return Py_BuildValue("s", s);
		}

free(functS);
free(xS);
free(weightAS);
Py_RETURN_NONE;

}





static PyObject* py_SimpsonError(PyObject* self, PyObject *args)
{
	double lenghtHSe;
	double SumSe = 0;
	
	double relerSE; //relative error
	double AnsSe = 1 - cos(1); //analytical answer.
	double functSE; 
	
	int nStepsSe; // number of steps
	int finStepsSe = 1000; // final value of steps
	//counterforloop with arrays, added E to indicate its from the error function
	//to prevent confusion while personally reading the code, so i wont accidentaly 
	//alter the wrong function 
	int countASe;
	
	
	double fromIntSe = 0.0; //left integral bound
	double toIntSe = 1.0;//right integral bound

//	double *weightASe ;//pointer for array added E as well
	int numlinesSe = 0;
//	double blah;//explained below, was used to troubleshoot error

    double integrand(double functE);
    
	


	FILE *Ersimp;
	Ersimp = fopen("relativEerrorSimp.txt","w");

	
	double *weightASe = (double *) malloc (finStepsSe * sizeof(double) );




    /*
     * The equation is:
     * b
     * ? f(x) dx =  h * f1 + 4h * f2 + 2h * f3 + 4h *f3 + ... + h * fN
     * a            3        3         3         3              3
     * h is the interval length, N is the number of terms,
     * f1 = f(a) and fN = f(b)
     *
     * The weights are defined as the term that is not the function.
     * h    4h  2h  4h ... 2h  4h  h
     * 2    3   3   3      3   3   3
     *
     * This loop sums the terms of the series (SumSe), using the weight Array (wieghtAS(impson)E(rror) and
     * passing the functionArgument to be evaluated appropriately. The number
     * of terms in the series is controlled, and incremented each time
     * until it reaches the maxNumberOfPoints.
     */


	for( nStepsSe = 2;	nStepsSe <  finStepsSe; nStepsSe++)
	{
	
			
		lenghtHSe = (toIntSe - fromIntSe)/(nStepsSe - 1);
		
	
		if(nStepsSe % 2==1) //modulo opperator 
		{

            /*
            * This loop creates the array for weighting the sum when doing the interval
            * it follows the weights as they are defined above, and assigns
            * weights based on the term number.
            */
		
			for(countASe = 0; countASe < nStepsSe ;countASe++)
			{
				if(countASe == 0 || countASe == (nStepsSe - 1))
				{
					weightASe[countASe] = lenghtHSe/3;
				}else if ((countASe % 2) == 1)
				{
					weightASe[countASe] = 4*lenghtHSe/3;
				}else
				{
					weightASe[countASe] = 2*lenghtHSe/3;
				}
			}
			//printf("counter: %d\n", countASe);

            /*
            * These two statements reset the parameters to initial conditions to
            * start the next loop.
            */
			
			SumSe = 0;
			functSE = fromIntSe;
			//print("# of points: %d\n", nStepsSe);
			
            /*
            * This loop actually does the summing from left most point (fromint(egral)) to right most point (toint).
            * It increments functionArgument by intervalLength each time
            */			

			for(countASe = 0; functSE < 1.0001*toIntSe; countASe++)
			{
				//printf("Funct Argument: %f\n", functSE);
				SumSe += weightASe[countASe] * integrand(functSE);
				functSE += lenghtHSe;
				//printf("counter in loop: "%d\n", countASe++);
			
			}
			
            /*
            * This calculates the relative error and prints it to file. 
            * Enter the analytic answer above.
            */			
			
			relerSE = fabs(( AnsSe - SumSe) / AnsSe );
			fprintf(Ersimp,"%d %e\n", nStepsSe ,relerSE);
			
			
			char *s = "Simpson Method Relative Error:";
	        //return Py_BuildValue("s", s);
			PyArg_ParseTuple(args, "dd", &relerSE);
            return Py_BuildValue("s, d",s, relerSE);
		
		}

	}

Py_RETURN_NONE;



}

double integrand(double functE)//integrand from example
{	
	
	
	double integrand ;

	integrand = exp(-functE) ;//6.25 can be changed
	
	return integrand ;
}

// Next Function to work on! 


static PyObject* py_functionFile(PyObject* self, PyObject* args) //makes the files read in for functions 
{
    double x;
    double nStep = .01; // Change the step size here
    double fromIntff = 0; //Change the left side of the integral here
    double toIntff = .99; // change the right side of the integral here
    double fVal; //functions value
   
    
	FILE *fp;
    fp = fopen("functionfile.txt", "w");
 
	/*
	*This loop writes the x tearms and functionValues to file.
	*/

       for(x = fromIntff; x < toIntff; x += nStep)
    {
        fVal = exp(-x); // functionfortextfile
        fprintf(fp, "%lf %lf\n", x, fVal);
    }
	 fclose(fp);
Py_RETURN_NONE;
}

/*
 * Gauss function. It maps the points, finds the weights,
 * allows the summation to be done in gaussQuadrature.
 */

void gauss(int npts, int job, double a, double b,double x[], double w[])
{
/*
*     npts     number of points
*     job = 0  rescaling uniformly between (a,b)
*     1  for integral (0,b) with 50% pts inside (0, ab/(a+b)
*     2  for integral (a,inf) with 50% inside (a,b+2a)
*     x, w     output grid points and weights.
*/

 int     m, i, j;
  double  t, t1, pp, p1, p2, p3;
  double  pi= 3.1415926535897932385E0;
  double  eps= 3.e-15; //adjust accuracy here
	  
	m= (npts+1)/2;
      for (i = 1; i<= m; i++) {
        t = cos(pi*(i-0.25)/(npts+0.5));
        t1= 1;
        while((fabs(t-t1))>= eps){
          p1= 1.0;
          p2= 0.0;
          for (j = 1; j<= npts; j++) {
            p3= p2;
            p2= p1;
            p1= ((2*j-1)*t*p2-(j-1)*p3)/j;
		   }
         pp= npts*(t*p1-p2)/(t*t-1);
         t1= t;
         t = t1 - p1/pp;
       }
       x[i-1]= -t;
       x[npts-i]= t;
       w[i-1]   = 2.0/((1-t*t)*pp*pp);
       w[npts-i]= w[i-1];
    }
    if (job == 0) {
      for (i = 0; i<npts ; i++) {
            x[i]= x[i]*(b-a)/2.0+(b+a)/2.0;
            w[i]= w[i]*(b-a)/2.0;
      }
    }
      if (job == 1)  {
        for (i = 0; i<npts; i++) {
          x[i]= a*b*(1+x[i]) / (b+a-(b-a)*x[i]);
          w[i]= w[i]*2*a*b*b /((b+a-(b-a)*x[i])*(b+a-(b-a)*x[i]));
       }
    }
      if (job == 2)   {
        for (i = 0; i<npts; i++)  {
          x[i]= (b*x[i]+b+a+a) / (1-x[i]);
          w[i]=  w[i]*2*(a+b)  /((1-x[i])*(1-x[i]));
         }
      }
}


/*
 * function estimates an integral using Gauss quadrature. It specifically
 * calls a function written by Landau for the weighting and mapping of points,
 * but it does the summation here.

 */
static PyObject* py_gaussQ(PyObject* self, PyObject* args)
{
    double sumOfTerms = 0;
    double leftMostPoint = 0; //Change the left of the integral here
    double rightMostPoint = 1; // Change the right side of the integral here
    int numberOfPoints = 100;
    int counterForArrays; // Used in loops for arrays
    void gauss(int npts, int job, double a, double b,double x[], double w[]);

    /*
     * These arrays will be passed to the gauss function and filled with
     * appropriate values.
     */
    double *yValue = (double *) malloc ( numberOfPoints * sizeof(double) );
    double *weightArray = (double *) malloc ( numberOfPoints * sizeof(double) );


    // This calls the gauss function to get the correct weights and points.

    gauss(numberOfPoints, 0, leftMostPoint, rightMostPoint, yValue, weightArray);


   /*
    * The following does job 0, in gauss, and then finds the sum.
    * Job 0 is called as it fits our arbitrary a to b integrals
    */

    for (counterForArrays = 0; counterForArrays < numberOfPoints; counterForArrays++)
    {
        sumOfTerms += weightArray[counterForArrays] * integrand(yValue[counterForArrays]);
    }

    printf("The sum for Gauss Job 0: %f\n", sumOfTerms);
    char *s = "Sum for Gauss Job 0:";
	PyArg_ParseTuple(args, "dd", &sumOfTerms);
    return Py_BuildValue("s, d",s, sumOfTerms);

    /*
     * The following would call jobs 1 and 2 and do the summation, but those
     * jobs are not applicable to the current problem. Job 1 is from 0 to b and
     * job 2 is from a to infinity.
     */

    // This calls the gauss function to get the correct weights and points.
    /*gauss(numberOfPoints, 1, leftMostPoint, rightMostPoint, yValue, weightArray);

    for (counterForArrays = 0; counterForArrays < numberOfPoints; counterForArrays++)
    {
        sumOfTerms+= weightArray[counterForArrays] * integrand(yValue[counterForArrays]);
    }

    printf("The sum for Gauss Job 1: %f\n", sumOfTerms);

    char *s = "Sum for Gauss Job 1:";
	PyArg_ParseTuple(args, "dd", &sumOfTerms);
    return Py_BuildValue("s, d",s, sumOfTerms);


    /*
     * The following does job 1, in gauss, and then finds the sum.
     *
    sumOfTerms = 0;
    // This calls the gauss function to get the correct weights and points.
    gauss(numberOfPoints, 2, leftMostPoint, rightMostPoint, yValue, weightArray);

    for (counterForArrays = 0; counterForArrays < numberOfPoints; counterForArrays++)
    {
        sumOfTerms+= weightArray[counterForArrays] * integrand(yValue[counterForArrays]);

    }

    printf("The sum for Gauss Job 2: %f\n", sumOfTerms);
    char *s = "Sum for Gauss Job 2:";
	PyArg_ParseTuple(args, "dd", &sumOfTerms);
    return Py_BuildValue("s, d",s, sumOfTerms);

    */



    free(yValue);
    free(weightArray);

}


/*
 * This function uses an increasing number of points to find the Gauss quadrature
 * approximations of the integral. It uses the function written by Landau and Paez. It
 * prints the number of points and relative error to file.
 */
static PyObject* py_gaussError(PyObject* self, PyObject* args)
{
    double sumOfTerms = 0;
    double leftMostPoint = 0; //Change the left of the integral here
    double rightMostPoint = 1; // Change the right side of the integral here
    int numberOfPoints = 10;
    int maxNumberOfPoints = 1000; // Change the max number of points here
    int counterForArrays; // Used in loops for arrays
    double relativeError;
    double analyticAnswer = 1 - exp(-1);
    void gauss(int npts, int job, double a, double b,double x[], double w[]);
    double integrand(double functE);
    
    
     //Change the file here
    FILE *fp;
    fp = fopen("relativeErrorGaussian.txt", "w");



    /*
     * These arrays will be passed to the gauss function and filled with
     * appropriate values.
     */
    double *yValue = (double *) malloc ( maxNumberOfPoints * sizeof(double) );
    double *weightArray = (double *) malloc ( maxNumberOfPoints * sizeof(double) );


    /*
     * This loop does succesive estimations of the integral, starting with a small
     * amount of points and incrementing each time. The results are printed to file
     * to be graphed.
     */

    for(numberOfPoints = 10; numberOfPoints < maxNumberOfPoints; numberOfPoints++)
    {
        // This calls the gauss function to get the correct weights and points.

        gauss(numberOfPoints, 0, leftMostPoint, rightMostPoint, yValue, weightArray);


       /*
        * The following does job 0, in gauss, and then finds the sum.
        * Job 0 is called as it fits our arbitrary a to b integrals
        */

        for (counterForArrays = 0; counterForArrays < numberOfPoints; counterForArrays++)
        {
            sumOfTerms += weightArray[counterForArrays] * integrand(yValue[counterForArrays]);
        }


        //This calculates the relative error
        relativeError = fabs((analyticAnswer - sumOfTerms)/analyticAnswer);

        fprintf(fp, "%d %e\n", numberOfPoints, relativeError); //Prints to file
        char *s = "Relative error";
	    PyArg_ParseTuple(args, "d, d", &relativeError );
        return Py_BuildValue("s, d, d", s ,relativeError );
       

        sumOfTerms = 0; // Resets sumOfTerms to 0 for the next loop.
    }

    /*
     * The following would call jobs 1 and 2 and do the summation, but those
     * jobs are not applicable to the current problem. Job 1 is from 0 to b
     * and job 2 is from a to infinity
     */

    /*// This calls the gauss function to get the correct weights and points.
    gauss(numberOfPoints, 1, leftMostPoint, rightMostPoint, yValue, weightArray);

    for (counterForArrays = 0; counterForArrays < numberOfPoints; counterForArrays++)
    {
        sumOfTerms+= weightArray[counterForArrays] * integrand(yValue[counterForArrays]);
    }

    printf("The sum for Gauss Job 1: %f\n", sumOfTerms);
    char *s = "Sum for Gauss Job 1:";
	PyArg_ParseTuple(args, "dd", &sumOfTerms);
    return Py_BuildValue("s, d",s, sumOfTerms);    

    /*
     * The following does job 1, in gauss, and then finds the sum.
     *
    sumOfTerms = 0;
    // This calls the gauss function to get the correct weights and points.
    gauss(numberOfPoints, 2, leftMostPoint, rightMostPoint, yValue, weightArray);

    for (counterForArrays = 0; counterForArrays < numberOfPoints; counterForArrays++)
    {
        sumOfTerms+= weightArray[counterForArrays] * integrand(yValue[counterForArrays]);

    }

    printf("The sum for Gauss Job 2: %f\n", sumOfTerms);
    char *s = "Sum for Gauss Job 2:";
	PyArg_ParseTuple(args, "dd", &sumOfTerms);
    return Py_BuildValue("s, d",s, sumOfTerms);
    */

    free(yValue);
    free(weightArray);
    fclose(fp);

}

/*
 * Integral of a function with monte carlo methods
 * Will integrate given amount of dimensions only from a to b on each dimension.
 */
static PyObject* py_MonteCarloMultiDIntegration (PyObject* self,PyObject* args)
{
    double meanValueSum = 0;
    double integrationValue = 0;
    double xValue;

    // Change the limits of integration here
    double leftMostPoint = 0;
    double rightMostPoint = 1;

    int counterForOuterLoop;
    int counterForInnerLoop;

    int maxNumberOfTrials = 80000;
    int numberOfXTerms = 10;
    int numberOfDimensions = 10;
    int seed = 48345;

    double errorOfIntegration;
    /*
     * analytic answer entered here
     */
    double analyticAnswer = 155/6;
    srand48(seed);
    
    double functionToBeIntegrated(double functionArgument);


    FILE *fp;
    fp = fopen("monteCarloIntegration.txt", "w");


     /* The equation:
      * b    d
      * ? dy ? f(x,y) dx = (b-a) * (d-c) * <f(x, y)>
      * a    c
      */
    for (counterForOuterLoop = 1; counterForOuterLoop <= maxNumberOfTrials; counterForOuterLoop++)
    {

        xValue = 0;

        for (counterForInnerLoop = 1; counterForInnerLoop <= numberOfXTerms; counterForInnerLoop++) 
//      xValue += drand48();
//      integrationValue += xValue*xValue;
       
       
        {
            xValue += (rightMostPoint - leftMostPoint) * drand48();
        }
        //meanValueSum divided by num of terms to find mean val of function
  	     meanValueSum += functionToBeIntegrated(xValue);

        /*
         * pow  = number of dimensions.
         */
  //     integrationValue = pow((rightMostPoint - leftMostPoint), numberOfDimensions) * meanValueSum / counterForOuterLoop;
     integrationValue = pow((rightMostPoint - leftMostPoint), numberOfDimensions) * meanValueSum / counterForOuterLoop;   
     
     errorOfIntegration = fabs(analyticAnswer - integrationValue)/analyticAnswer;

       // Print to file
     //  fprintf(fp, "%d %f\n", counterForOuterLoop, errorOfIntegration);
    }
//    printf("The monte carlo integration value is: %f\n", integrationValue/counterForOuterLoop);
	  printf("The monte carlo integration value: %f\n", integrationValue);
	    char *s = "Monte Carlo Integration Value";
	    PyArg_ParseTuple(args, "d", &integrationValue );
        return Py_BuildValue("s, d", s ,integrationValue );
	printf("Error of Integration is : %f\n ", errorOfIntegration); 


    fclose(fp);

}


/*
 * Error analysis program.
 * Accepts argument , returns the function Value of that argument. used for multi D monte carlo integration
 */
double functionToBeIntegrated(double functionArgument)
{
    double functionToBeIntegrated;

    functionToBeIntegrated = exp(functionArgument); //function varies here , used for 10 d monte carlo, function arg = sum of Xs

    return functionToBeIntegrated;
}

/*
 * This simple function accepts an argument and returns the function Value of
 * that argument. Used for error analysis and in the quadrature program.
 */
double functDiff (double value)
	{
		
		double functDiff;
		functDiff = sqrt(30 - value) * tan(sqrt(30 - value)) - sqrt(value);// User change function here**************************************************
		return functDiff; 
	}
	
	
/*
 * This function approximates the derivative using the extrapolate differentiation
 * method.
 */

static PyObject* py_extrapDiff(PyObject* self, PyObject* args)
{
	double stepsizee;
	double diffe;
	double difflocatione;
	double precisione ;
	double anse; //analytic answer 
	double relEre;//relative error
	char g;
	double functDiff(double value);
	
	
	FILE *extr;
	extr = fopen("extrSinEr01.txt","w");

	difflocatione = 1.0 ;
	anse = - sin(0.1);
	precisione = pow(10.0 , -15.0);

	fprintf(extr, "StepSize      Relative Error\n");

     /*
     * The extrap differentiation formula is:
     *
     * dy = ((8 * (y(x+h/4)-y(x-h/4)) - (y(x + h/2) -  y(x-h/2)))
     * dt                         3 * h
     *
     * where h is the interval length and x is the point to be differentiated at
     */



	for (stepsizee = .001; stepsizee > precisione; stepsizee /= 1.01)
	{


		diffe = (8 * (functDiff (difflocatione + stepsizee/4) - functDiff(difflocatione - stepsizee/4)) -  (functDiff (difflocatione + stepsizee/2) - functDiff(difflocatione - stepsizee/2))) / (stepsizee * 3);

        /*
         * Enable the commented printf statements to see the values on the console,
         * otherwise they are printed to file.
         */

		printf("extrapolate diff: step %e val %f\n", stepsizee, diffe); 
		relEre = fabs((anse - diffe)/anse);
	    char *s = "Extrapolate Differentiation answer :";
	    PyArg_ParseTuple(args, "d", &diffe );
        return Py_BuildValue("s, d", s ,diffe );
		
		
		printf("relError %e\n",relEre);
		//prints to file
		fprintf(extr, "%e %e\n", stepsizee , relEre);

	}

fclose(extr);

}
	
/*
 * This function approximates the second derivative using the central differentiation
 * more step method. Less prone to subtractive cancellation.
 */

static PyObject* py_secderivone(PyObject* self, PyObject* args) //moresteps
{

     /*
     * Please input the point to be differentiated at, and, for error analysis,
     * enter the analytic answer.
     */

	double stepsizedd;
	double diffdd;
	double diffddleast;
	double difflocationdd;
	double precisiondd ;
	double ansdd; //analytic answer 
	double relErdd;//relative error
	char g;
    double relErddpast = pow(10.0,12.0);
    
	 
	double functDiff(double value);
	FILE *ddone;
	ddone = fopen("sescDerivErrorMore01.txt","w");

	//uses forward diff formula 
	// 
	difflocationdd = 0.1 ;
	ansdd = - cos(0.1); //analytic Answer
	precisiondd = pow(10.0 , -14.0);

	fprintf(ddone, "StepSize      Relative Error\n");
     /*
     * The second derivative differentiation formula is:
     *
     * dy = y(x + h) - y(h) -  (y(h) - y(x-h))
     * dt                  h
     *
     * where h is the interval length and x is the point to be differentiated at
     */


	for (stepsizedd = .1; stepsizedd > precisiondd; stepsizedd /= 1.01)
	{

		diffdd = ((functDiff (difflocationdd + stepsizedd) - functDiff(difflocationdd)) - (functDiff(difflocationdd) - functDiff(difflocationdd - stepsizedd))) / (stepsizedd*stepsizedd);

        /*
         * Enable the commented printf statements to see the values on the console,
         * otherwise they are printed to file.
         */

         printf("Central differentiation: %f at step %e\n", diffdd , stepsizedd);
		//calculates relative error using the analytic answer and the differential value calculated
		relErdd = fabs((ansdd - diffdd)/ansdd);
		
		//printf("relError %e\n",relErdd);
		//prints to file
		fprintf(ddone, "%e %e\n", stepsizedd , relErdd);
		//scanf("%c",&g);
	 
	    if( relErdd <= relErddpast){
            relErddpast = relErdd;
        
            if(diffdd != 0.0){
        
                diffddleast = diffdd;
                }
        }
	 
	        
	
	}
	printf(" Least relError %e diffddt Value %e\n",relErddpast , diffddleast);
    char *s = "Second Derivative answer :";
	PyArg_ParseTuple(args, "d", &diffdd );
    return Py_BuildValue("s, d", s ,diffdd ); 
     
     
	fclose(ddone);
	//scanf("%c",&g);//so wont close
}

/*
 * This function approximates the second derivative using the central differentiation
 * less step method. More prone to subtractive cancellation
 */

static PyObject* py_secderivtwo(PyObject* self, PyObject* args) //less steps
{
     /*
     * Please input the point to be differentiated at, and, for error analysis,
     * enter the analytic answer.
     */
	double stepsizeddt;
	double diffddt;
	double diffddtleast;
	double difflocationddt = 0.1 ;
	double precisionddt = pow(10.0 , -14.0);
	double ansddt = - cos(0.1); //analytic answer 
	double relErddt;//relative error
	char g;
	double relErddtpast = pow(10.0,12.0);
	const char* test;
	//int ert;
	//char foobar[];

	
	FILE *ddtwo;
	ddtwo = fopen("secondDerivLess01.txt","w");
	double functDiff(double value);
	//uses forward diff formula 
	// 
	
	
	

	fprintf(ddtwo, "StepSizeRelative Error\n");


     /*
     * The second derivative differentiation formula is:
     *
     * dy = y(x + h) + y(x-h) - 2y(h)
     * dt               h
     *
     * where h is the interval length and x is the point to be differentiated at
     */



	for (stepsizeddt = .1; stepsizeddt > precisionddt; stepsizeddt /= 1.01)
	{

		diffddt = ((functDiff (difflocationddt + stepsizeddt) + functDiff(difflocationddt - stepsizeddt) - 2 * functDiff(difflocationddt))) / (stepsizeddt*stepsizeddt);

		//calculates relative error using the analytic answer and the differential value calculated
		
        /*
         * Enable the commented printf statements to see the values on the console,
         * otherwise they are printed to file.
         */	
      //   printf("Central differentiation: %f at step %e\n", diffddt , stepsizeddt);	
		
		relErddt = fabs((ansddt - diffddt)/ansddt);
		
	//	printf("relError %e diffddt vlue %e\n",relErddt , diffddt);
		//prints to file
		fprintf(ddtwo, "%e %e %e\n", stepsizeddt , relErddt , diffddtleast);
		//scanf("%c",&g);
        
        
        if( relErddt <= relErddtpast){
        relErddtpast = relErddt;
        
            if(diffddt != 0.0){
        
            diffddtleast = diffddt;
            }
        }


	}
    printf(" Least relError %e diffddt Value %e\n",relErddtpast , diffddtleast);
  //  PyArg_ParseTuple(args, "s", &test);
   // ert = system(test);

    char *s = "Second Derivative version two answer :";
	//PyArg_ParseTuple(args, "d", &diffddt );
   // return Py_BuildValue("s, d", s ,diffddt ); 
   
   	PyObject *m = PyImport_AddModule("__main__");
    PyObject *v = PyObject_GetAttrString(m,"foobar");

    int foobar = PyInt_AsLong(v); //reads value from python object foobar , use to get user information
    
    
    
    printf("input was %d \n", foobar ); 
    Py_DECREF(v);

   
   
   
   
   
   
    return Py_BuildValue("s", s ); 
	fclose(ddtwo);
	//scanf("%c",&g);//so wont close
}


static PyObject* py_newtonRaph(PyObject* self, PyObject* args)
{
    double fDiff ( double fdifflocation, double fxstepsize);
    /*
     * Please enter an initial guess near the zero.
     */
	
	double inVal = 28.25; // enter initial guess near zero
	double xstepsize = .01;
	double precision = pow(10,-13.0);

	/*
	* The iterative equation for the Newton Raphson method is:
	*
	* p_[n+1] = p_[n]  -  f(p_[n])
	*                     f'(p_[n])
	*
	* where p_[n+1] is the next guess, and p[n] is the current guess
	*/


	while (fabs(functDiff(inVal)) > precision)
	{
		inVal = - functDiff(inVal) / fDiff(inVal, xstepsize); //fDiff is forward differentiation 


               /*
                * This is the backtracking loop. It checks to make sure the guess is not
                * leading into an infinite loop. If it does, it makes a different
                * guess.
                */




            while ( pow(fabs(functDiff(inVal + xstepsize)), 2) > pow(fabs(functDiff(inVal)),2))
            {
                xstepsize /= 2 ;
            }	
	
	    inVal += xstepsize;
	}

	printf("zero point: %1.12f \n", inVal);
    char *s = "ZERO Point:";
	PyArg_ParseTuple(args, "d", &inVal );
    return Py_BuildValue("s, d", s , inVal ); 

}

/*
 * This simple function uses the forward differentiation method to find the
 * approximate derivative. It just returns the value.
 */


double fDiff ( double fdifflocation, double fxstepsize)
{

	double diffVal;
	double yVal;

     /*
     * The forward differentiation formula is:
     *
     * dy = y(x + h) -  y(x)
     * dt          h
     *
     * where h is the interval length and x is the point to be differentiated at
     */



	yVal = fdifflocation + fxstepsize;
	
	diffVal = (functDiff(yVal) - functDiff(fdifflocation)) / fxstepsize;

	return diffVal;

}









/* ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////






 








/*
 * Bind Python function names to our C functions
 */
static PyMethodDef myModule_methods[] = {
    {"newtonRaph", py_newtonRaph, METH_VARARGS},
    {"secderivtwo", py_secderivtwo,METH_VARARGS},
    {"secderivone", py_secderivone, METH_VARARGS},
    {"extrapDiff", py_extrapDiff, METH_VARARGS},
    {"MonteCarloMultiDIntegration", py_MonteCarloMultiDIntegration, METH_VARARGS},
    {"gaussError", py_gaussError, METH_VARARGS},
    {"gaussQ", py_gaussQ, METH_VARARGS},
    {"functionFile",py_functionFile, METH_VARARGS},
    {"SimpsonError",py_SimpsonError, METH_VARARGS},
	{"SimpsonRule",py_SimpsonRule, METH_VARARGS},
	{NULL, NULL}
};
/*
 * Python calls this to let us initialize our module
 */
void initmyModule()
{
	(void) Py_InitModule("myModule", myModule_methods);
}