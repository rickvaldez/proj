//Rick V. 
//440

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

double functionToBeIntegrated(double functionArgument);
void linearCongruentialRandomNumbers ();
void standardRandom();
void TwoDRandomWalk ();
void radioactiveDecay ();
void MonteCarloMultiDIntegration ();
void pond ();


/*
 * The main function calls other functions. 
 * Remove comments to run function and compile
 */
int main() {

   // linearCongruentialRandomNumbers();

   // standardRandom();

   // TwoDRandomWalk();

   //radioactiveDecay();

   //pond();

   MonteCarloMultiDIntegration();

}


/*
 *	Montecarlo methods function to determine the area of a pond by randomly tossing pebbles. 
 */


void pond()
{
    int max = 1000; // Change the max stones thrown here


    int seed = 68111; //Change the random seed here
    srand48(seed);
    int i, pi = 0;
    double x, y, area;

    FILE *fp;
    fp = fopen("pond.txt","w"); // writes to txt file

    for(i=1;i<=max;i++)
    {
    	x = drand48()*2-1;   //creates floats between 1 and -1
    	y = drand48()*2-1;
    	if((x*x + y*y) < 1) pi++;   //stone hit the pond condition
        area = 4 * (double)pi/i;     // calculates area

	fprintf(fp,"%i %f\n",i,area); //prints to txt file
    }
    fclose(fp);


}


/* 
 * Linear Congruent method for the generation of random number along with a moment generating randomness check. 
 */
void linearCongruentialRandomNumbers ()
{
    int seedNumber = 10; // Change the seed here
    int lastRandomNumber; // previous random number for check
    int constantMultiplyByA = 57; //parameter for lin cong meth
    int constantAddByC = 1; //parameter for lin cong meth
    int constantModByM = 256; //paremeter for lin cong meth
    int numberOfPointsToGenerate = 100001;
    int counter = 0;

    /*
     * Moment generating check Variables
     */
    double momentGeneratingCheck = 0;
    double momentGeneratingCheckExponent = 1;
    double momentGeneratingSum = 0;
    double scaledRandomNumber;


    FILE *fp;
    fp = fopen("randomNumberInfo.txt", "w"); // Open file to write into

    /*
     * Moment Generating equation
     *
     * r[next] = (a*r[last] + c) % M
     *
     * as previously said in the int declaration the parameters belonging to  lin cong method
     * a is constantMuliplyByA, c is constantAddByC, and M is constantModByM.
     * r[next] is nextRandomNumber, r[last] is lastRandomNumber
     *
     * The following code applied the equation using the seed as the initial value for the next number as the equation depends on the previous value.
     * code in the loop generates the succesive points.
     */

    lastRandomNumber = (constantMultiplyByA * seedNumber + constantAddByC) % constantModByM;
    //printf("The random number %d is %d. \n", counter, lastRandomNumber);
    fprintf(fp, "%d %d\n", counter, lastRandomNumber);


    for (counter = 1; counter < numberOfPointsToGenerate; counter++)
    {

        lastRandomNumber = (constantMultiplyByA * lastRandomNumber + constantAddByC) % constantModByM;

 
       /*
         * For a random number between 0-1,scaledRandomNumber is used.
         */
        scaledRandomNumber = (double) lastRandomNumber/constantModByM;


        //printf("The random number %d is %d. \n", counter, lastRandomNumber);
        fprintf(fp, "%d %d\n", counter, lastRandomNumber);



        /*
         * The equation does a check for randomness at 3 different N values.
         *                                N  k
         *         sqrt(N) * abs{ (1/N) * Σ x_i  -   [1 / (K+1)]   }
         *                               i=1        
         */

        momentGeneratingSum += pow(scaledRandomNumber, momentGeneratingCheckExponent);


        /*
         * The moment gen checks at done at 100, 10000, 100000
         */
        if (counter == 100 || counter == 10000 || counter == 100000)
        {
            momentGeneratingCheck = (double) sqrt(counter) * fabs(momentGeneratingSum/counter - 1 / (momentGeneratingCheckExponent + 1));
            printf("Linear Congruential: Moment Generating Check at N = %d: %e\n", counter, momentGeneratingCheck);

        }

    }

    fclose(fp);


}

/*
 * Uses rand48() function to generate random numbers. 
 * + moment gen check for randomness	
 */
void standardRandom()
{
    int counter;
    int numberOfPointsToGenerate = 100001;
    // Used in moment generating check
    double momentGeneratingCheck = 0;
    double momentGeneratingCheckExponent = 1;
    double momentGeneratingSum = 0;

    int seedForSrand48 = 123845; // change the seed here

    srand48(seedForSrand48); //Seeds for drand48

    FILE *fp;
    fp = fopen("standardRandomNumberInfo.txt", "w");

    /*
     *  Loop prints to file counter and rand() for graphing
     *  for randmomness check on drand48
     */

    for(counter = 0; counter < numberOfPointsToGenerate; counter++)
    {
        fprintf(fp, "%d %d\n", counter, rand());

        /*
         * The following statement does a check for randomness at three different N values.
         * The equation is:
         *                    N  k
         * sqrt(N) * abs( 1 * Σ x_i  -   1
         *                N   1         k+1
         *
         */
        momentGeneratingSum += pow(drand48(), momentGeneratingCheckExponent) ;

        /*
         * The moment generating check is only done at 100, 10000, 100000
         */
        if (counter == 100 || counter == 10000 || counter == 100000)
        {
            momentGeneratingCheck = (double) sqrt(counter) * fabs(momentGeneratingSum/counter  - 1/(momentGeneratingCheckExponent + 1));
            printf("Standard Random: Moment Generating Check at N = %d: %e\n", counter, momentGeneratingCheck);

        }


    }
}

/*
 * 2d random walk simulation for any given amount of steps. 
 *using drand48()
 *using random number between -1 and 1 adds that to the prior number, representing a random step in any direction for x and y axis.
 * determines average distance traveled.
 */
void TwoDRandomWalk ()
{
    /*
     * Set the initial 2D location here
     */
    double xLocation = 0;
    double yLocation = 0;

    double scalarDistanceTravelled;
    double sumOfScalarDistancesTravelled = 0;

    double averageOfScalarDistancesTravelled;

    int counterForInnerLoop;
    int numberOfWalks = 50; // num of walks
    int counterNumberOfWalks;
    int seedValue = 1234; // seed to change walk
    int numberOfSteps = 100;
    int squarerootNumsteps;
    int counterNumberOfSteps = 1;
    srand48(seedValue);

    FILE *fp;
    fp = fopen("randomWalk1.txt", "w");

    FILE *gp;
    gp = fopen("rootNumStepsvsAvgDist.txt" , "w");

/* 
 * Loop is used to find the average scalar distnce traveled for a given number of steps.
 */

for ( counterNumberOfSteps = 1; counterNumberOfSteps < numberOfSteps; counterNumberOfSteps++)
{
    /*
     * random walk loop using the random number gen. drand48
     * between 0 and 1. 
     * minus 0.5 gives a range between -0.5 and 0.5,
     * and multiplying by 2 gives range from  -1 and 1. 
     * The outer loop used to manage multiple walks.
     */
    for (counterNumberOfWalks = 0; counterNumberOfWalks < numberOfWalks; counterNumberOfWalks++)
    {
         /*
          * Resets location for neext walk
          */
        xLocation = 0;
        yLocation = 0;
        fprintf(fp, "\n");

        /*
         * loop executes the walking algorithm.
         */
       for (counterForInnerLoop = 0; counterForInnerLoop < numberOfSteps; counterForInnerLoop++)
        {
           /*
            * -.5 and *2 allow negative numbers to be possibilities
            */
            xLocation += (drand48() - 0.5) * 2;
            yLocation += (drand48() - 0.5) * 2;
            fprintf(fp, "%f %f\n", xLocation, yLocation);
        }

        // Print distance travelled per walk to screen
        scalarDistanceTravelled = xLocation * xLocation + yLocation * yLocation;

 //       fprintf(gp," %d %f \n", counterForInnerLoop, scalarDistanceTravelled);

        //Used for average distance travelled
        sumOfScalarDistancesTravelled += scalarDistanceTravelled;

        //Reseed
        seedValue += 100;
        srand48(seedValue);
    }

    squarerootNumsteps = sqrt(counterNumberOfSteps);
    averageOfScalarDistancesTravelled = sumOfScalarDistancesTravelled / numberOfWalks;
    fprintf(gp, "%d %f\n", squarerootNumsteps, averageOfScalarDistancesTravelled);

}

    fclose(gp);
    fclose(fp);



}

/*
 * radioactive decay simulation fucntion. 
 * algorithm loops through every nuclei in the mass.
 * evaluate random number vs the decay rate
 * If random number is smaller than decay rate, that nuclei decays. 
 * After every nuclei sorted by the loop the unit time measurement is increased by one. 
 */
void radioactiveDecay ()
{

    /*
     * Initial conditions totalNumberOfNucleiLeft should be equal numberOfNucleiRemaining 
     */
    int totalNumberOfNuclei = 10000;
    int totalNumberOfNucleiLeft = 10000;
    int counterNumberOfNuclei;
    int numberOfNucleiRemaining = 10000;
    int numberOfNucleiDecayedThisTime;
    double logNumNucRem = 0.0;
    double logNumNucDec = 0.0 ;
    double logFractOfNuc = 0.0 ;   
    double fractionOfNuclei = 0.0;
    double decayRate = 0.1;
    double randomDecay;
    int time;
    int maxTime = 100; // max time , use to vary time
    int seed = 1234; // random seed
    srand48(seed);


    FILE *fp;
    fp = fopen("decayDataNumRem.txt", "w");
   
    FILE *hp;
    hp = fopen("decayDataNumNucleiDec.txt","w");

    FILE *jp;
    jp = fopen("logFractRemVsTime.txt", "w");   

    FILE *kp;
    kp = fopen("logNumNucleiDecay.txt", "w");

    FILE *lp;
    lp = fopen("logNucleiLeftVsLogNucleiDec.txt","w");

    /*
     * Outer loop sorts through steps in time.
     * The inner loop executes a decay. 
     * Decays occur if the random number generated is less than decay rate.
     */
    for (time = 0; time < maxTime; time++)
    {
        numberOfNucleiDecayedThisTime = 0;

        /*
         * Sorts each atom individually, Decays occur if the random number generated is less than decay rate.
         */
        for (counterNumberOfNuclei = 1; counterNumberOfNuclei < totalNumberOfNucleiLeft; counterNumberOfNuclei++)
        {
            randomDecay = drand48();

            /*
             * This performs the decay if it is appropriate.
             */
            if(randomDecay < decayRate)
            {
                numberOfNucleiRemaining--;
                numberOfNucleiDecayedThisTime++;
            }
            

            //printf("%d\n", counterNumberOfNuclei);
        }

        totalNumberOfNucleiLeft = numberOfNucleiRemaining;

        // used for graph
        fractionOfNuclei = (double)numberOfNucleiRemaining/totalNumberOfNuclei;
		


        /*
         * Print to file.
         */
        fprintf(fp, "%d %d\n", time, numberOfNucleiRemaining);
	fprintf(hp, "%d %d\n", time , numberOfNucleiDecayedThisTime);

	logFractOfNuc = log( fractionOfNuclei);
	logNumNucRem = log( numberOfNucleiRemaining);
	logNumNucDec = log(numberOfNucleiDecayedThisTime);

	fprintf(jp, "%d %g\n", time , logFractOfNuc);
	fprintf(kp, "%d %g\n", time , logNumNucDec);
	fprintf(lp, "%g %g\n", logNumNucRem , logNumNucDec );

    }
    fclose(fp);
    fclose(hp);
    fclose(jp);
    fclose(kp);
    fclose(lp);
}


/*
 * Integral of a function with monte carlo methods
 * Will integrate given amount of dimensions only from a to b on each dimension.
 */
void MonteCarloMultiDIntegration ()
{
    double meanValueSum = 0;
    double integrationValue = 0;
    double xValue;

    // Change the limits of integration here
    double leftMostPoint = 0;
    double rightMostPoint = 1.0;

    int counterForOuterLoop;
    int counterForInnerLoop;

    int maxNumberOfTrials = 8192;
    int numberOfXTerms = 10;
    int numberOfDimensions = 10;
    int seed = 41345;

    double errorOfIntegration;
    /*
     * analytic answer entered here
     */
    double analyticAnswer = 155/6;
    srand48(seed);

    FILE *fp;
    fp = fopen("monteCarloIntegration.txt", "w");


     /* The equation:
      * b    d
      * ∫ dy ∫ f(x,y) dx = (b-a) * (d-c) * <f(x, y)>
      * a    c
      */
    for (counterForOuterLoop = 1; counterForOuterLoop < maxNumberOfTrials + 1; counterForOuterLoop++)
    {

        xValue = 0;

        for (counterForInnerLoop = 0; counterForInnerLoop < numberOfXTerms; counterForInnerLoop++)
        {
            xValue += (rightMostPoint - leftMostPoint) * drand48();
        }

        //meanValueSum divided by num of terms to find mean val of function
        meanValueSum += functionToBeIntegrated(xValue);

        /*
         * pow  = number of dimensions.
         */
       integrationValue = pow((rightMostPoint - leftMostPoint), numberOfDimensions) * meanValueSum / counterForOuterLoop;
       errorOfIntegration = fabs(analyticAnswer - integrationValue)/analyticAnswer;

       // Print to file
       fprintf(fp, "%d %f\n", counterForOuterLoop, errorOfIntegration);
    }
    printf("The monte carlo integration value is: %f\n", integrationValue);

    fclose(fp);

}


/*
 * Error analysis program.
 * Accepts argument , returns the function Value of that argument. 
 */
double functionToBeIntegrated(double functionArgument)
{
    double functionToBeIntegrated;

    functionToBeIntegrated = exp(functionArgument); //function varies here

    return functionToBeIntegrated;
}
