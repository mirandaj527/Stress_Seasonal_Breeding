// **********************************************************************************
// Dynamic programming model of stress response including seasonal breeding,
// predator autocorrelation and damage.
//
// 
// **********************************************************************************


//HEADER FILES

#include <cstdlib>
#include <stdio.h>
#include <ctime>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <fstream>
#include <cmath>
#include <random>
#include <vector>
#include <algorithm>
#include <chrono>
#include <string>



// constants, type definitions, etc.

using namespace std;

const int seed        = time(0); // pseudo-random seed
const double lambdaA  = 0.035;   // probability that predator arrives
const double lambdaL  = 0.065;   // probability that predator leaves
const double pAtt     = 0.5;     // probability that predator attacks if present
const double alpha    = 0.1;     // parameter controlling effect of hormone level on pKill
const double beta_b   = 0.0;     // parameter controlling effect of hormone level on reproductive rate
const double kappa    = 0.0;     // Parameter controlling affect of damage on mortality
const double mu       = 0.02;   // background mortality (independent of hormone level and predation risk)
const double rho      = 1.0;    // Fixed rate of repair
const double h0       = 20.0;   // Reference hormone level
const double omega    = 0.1;   // Effect of deviations from h0 on damage build-up
const double gamma_g  = 0.1;     // Effect of damage on reproductive output
const int maxI        = 100000; // maximum number of iterations
const int maxT        = 25;     // maximum number of time steps since last saw predator
const int maxD        = 100;    // Number of discrete damage levels?
const int maxH        = 100;    // maximum hormone level
const int maxS        = 5;       // Length of the breeding cycle
const int skip        = 10;       // interval between print-outs

// Create a random engine with your chosen seed
static std::mt19937 mt(seed);

ofstream outputfile; // output file
stringstream outfile; // for naming output file

//int hormone[maxT][maxD][maxS];          // hormone level (strategy)                                   
vector<                                   // vector 1: maxT
  vector<                                 // vector 2: maxD                                 
    vector<int>>>                         // vector 3: maxS
      hormone(maxT,                                        
        vector<                                
          vector<int>> 
            (maxD, 
              vector<int> (maxS, 0.0)     // initialising internal vector of length maxS to contain 0s
            )
);
// storing optimal hormone level for next time step
vector<                                   // vector 1: maxT
  vector<                                 // vector 2: maxD                                 
    vector<int>>>                         // vector 3: maxS
      hormone_next(maxT,                                        
        vector<                                
          vector<int>> 
            (maxD, 
              vector<int> (maxS, 0.0)     // initialising internal vector of length maxS to contain 0s
            )
);
//double pKill[maxH];                     // probability of being killed by an attacking predator
vector<double> pKill(maxH, 0.0);
//double repro[maxD][maxH][maxS];         // reproductive output
vector<                                   // vector 1: maxD
  vector<                                 // vector 2: maxH
    vector<double>>>                      // vector 3: maxS
      repro(maxD,
        vector<
          vector<double>>
            (maxH,
              vector<double> (maxS, 0.0)
            )
);
//double Wopt[maxT][maxD][maxS];          // fitness immediately after predator has/hasn't attacked, under optimal decision h
vector<                                   // vector 1: maxT
  vector<                                 // vector 2: maxD                                 
    vector<double>>>                         // vector 3: maxS
      Wopt(maxT,                                        
        vector<                                
          vector<double>> 
            (maxD, 
              vector<double> (maxS, 0.0)       
            )
);
//double W[maxT][maxD][maxH][maxS];       // expected fitness at start of time step, before predator does/doesn't attack
vector<                                   // vector 1: maxT
  vector<                                 // vector 2: maxD
    vector<                               // vector 3: maxH
      vector<double>>>>                   // vector 4: maxS
        W(maxT,
          vector<
            vector<
              vector<double>>>
                (maxD,
                  vector<
                    vector<double>>
                      (maxH,
                        vector<double> (maxS, 0.0)
                      )
                )
);
//double Wnext[maxT][maxD][maxH][maxS];   // expected fitness at start of next time step
vector<                                   // vector 1: maxT
  vector<                                 // vector 2: maxD
    vector<                               // vector 3: maxH
      vector<double>>>>                   // vector 4: maxS
        Wnext(maxT,
          vector<
            vector<
              vector<double>>>
                (maxD,
                  vector<
                    vector<double>>
                      (maxH,
                        vector<double> (maxS, 0.0)
                      )
                )
);
//double pPred[maxT];                     // probability that predator is present
vector<double> pPred(maxT, 0.0);
//double damage_new[maxD][maxH];          // array of current damage
vector<                                   // vector 1: maxD
  vector<double>>                         // vector 2: maxH
    damage_new(maxD,
      vector<double> (maxH, 0.0)
);
//double background_mortality[d]
vector<double> background_mortality(maxD, 0.0);

double totfitdiff;                        // fitness difference between optimal strategy in successive iterations
double maxfitdiff;
double tothormonediff;
double maxhormonediff;

int i;     // iteration

int c;     // count


/* SPECIFY FINAL FITNESS */
void FinalFit()
{
  int t,d,h,s;

  for (t=1;t<maxT;t++) // note that Wnext is undefined for t=0 because t=1 if predator has just attacked
  {
    for (d=0;d<maxD;d++)
    {
      for (h=0;h<maxH;h++)
      {
        for (s=0;s<maxS;s++)
          {
            Wnext[t][d][h][s] = 1.0;
          }
      }
    }
  }

}


/* CALCULATE PROBABILITY THAT PREDATOR IS PRESENT */
void Predator()
{
  int t;

  pPred[1] = 1.0-lambdaL; // if predator attacked in last time step
  for (t=2;t<maxT;t++) // if predator did NOT attack in last time step
  {
    pPred[t] = (pPred[t-1]*(1.0-pAtt)*(1.0-lambdaL)+(1.0-pPred[t-1])*lambdaA) / (1.0 - pPred[t-1]*pAtt);
  }

}


/* CALCULATE PROBABILITY OF BEING KILLED BY AN ATTACKING PREDATOR */
void Death()
{
  int h;

  for (h=0;h<maxH;h++)
  {
   // pKill[h] = 1.0 - pow(static_cast<double>(h)/static_cast<double>(maxH), alpha);
      pKill[h] = 1.0 - (alpha * (static_cast<double>(h)/static_cast<double>(maxH)));
  }

}


/* CALCULATE PROBABILITY OF REPRODUCING */
void Reproduction()
{
  int d,h,s;

  for (d=0;d<maxD;d++)
  {
    for (h=0;h<maxH;h++)
    {
        // EDIT: Tweaked reproduction function. It still includes affect of hormone AND of damage
        // could look at only includign affect of damage & not hormone level (avoids paying cost of hormones twice)
      repro[d][h][0] = exp(-(beta_b*(static_cast<double>(h)/static_cast<double>(maxH)) 
                            + gamma_g*(static_cast<double>(d)/static_cast<double>(maxD))));
      for (s=1;s<maxS;s++) 
      {
          repro[d][h][s] = 0.0;
      }
    }
  }

}


/* CALCULATE DAMAGE IN NEXT TIME STEP*/
// note d(t+1) is only a function of d(t) and h, which itself is a function of season
void Damage()
{
  int d,h;

  for (d=0;d<maxD;d++)
  {
    for (h=0;h<maxH;h++)
    {
      damage_new[d][h] = d + omega*pow(h - h0, 2) - rho;
      damage_new[d][h] = max(0.0, min(static_cast<double>(maxD-1), damage_new[d][h])); 
      //max function ensures damage level never goes below 0, min function ensures it never goes above maxD
    }
  }
}


/* CALCULATE MORTALITY BASED ON DAMAGE*/
void Mortality()
{
    int d;

    for (d=0;d<maxD;d++)
    {
        background_mortality[d] = min(1.0, (mu + (1.0 - exp(-kappa*(static_cast<double>(d)/static_cast<double>(maxD))))));
    }
}

/* CALCULATE OPTIMAL DECISION FOR EACH t */
void OptDec()
{
  int t,d,h,s;
  double fitness;
  
  int d1, d2;
  double ddif;

  // calculate optimal decision h given current t (N.B. t=0 if survived attack)
  for (t=0;t<maxT;t++)          // loop over all t
  {
    for (d=0;d<maxD;d++)        // loop over damage levels
    {
      for (s=0;s<maxS;s++)      // loopover breeding seasons
      {
        Wopt[t][d][s] = 0.0;    // initialise optimal fitness array with all 0s (all fitness values = 0)
        hormone[t][d][s] = 0;   // initialise all hormone strategies as 0
                                //
        for (h=0;h<maxH;h++)    // for every possible hormone level...
        {
          // ORIGINAL FITNESS CODE:
          // fitness = Wnext[min(maxT - 1, t + 1)]      // fitness as a function of t, d, s, h
          //                    [damage_new[d][h]]
          //                    [hormone[min(maxT - 1, t + 1)][damage_new[d][h]][(s + 1) % maxS]]
          //                    [(s + 1) % maxS];      // (s + 1) % maxS resets S back to 0 after it reaches maxS - 1

          // NEW FITNESS CODE:
            d1=floor(damage_new[d][h]); // lower integer value of damage
            d2=ceil(damage_new[d][h]); // top integer value
            ddif = damage_new[d][h] - static_cast<double>(d1); // calculate difference 

            // calculate fitness from W' as a function of t, d, s, h
            fitness = 
                (1.0 - ddif) * Wnext[min(maxT - 1, t + 1)][d1][h][(s + 1) % maxS] // deterministic rounding
                + 
                ddif * Wnext[min(maxT - 1, t + 1)][d2][h][(s + 1) % maxS];
	      

        // compare with current optimal fitness for this specific combination of t,d,s
          if (fitness>Wopt[t][d][s])                 
          {
            Wopt[t][d][s] = fitness;                 // overwrite Wopt, a particular hormone value yields better fitness for this combination of t, d, s
            hormone[t][d][s] = h;                    // store hormone value that provides this best fitness value for this combination of t,d,s
          }
        }
      }
    }
  }

  // calculate expected fitness as a function of t d, h, and s, before predator does/doesn't attack
  for (t=1;t<maxT;t++) // note that W is undefined for t=0 because nothing happens at t=0, predator is attacks. t=1 if predator has just attacked
  {
    for (d=0;d<maxD;d++)
    {
      for (h=0;h<maxH;h++)
      {
        for (s=0;s<maxS;s++)
        {
          W[t][d][h][s] = pPred[t]*pAtt*(1.0-pKill[h])*(1-(background_mortality[d]))*
                          (Wopt[0][d][s]+repro[d][h][s])     // survive attack
                        + (1.0-pPred[t]*pAtt)*(1-(background_mortality[d]))*
                          (Wopt[t][d][s]+repro[d][h][s]);  // no attack

          // Checking to see if fitness values are sensible:
          // cout << W[t][d][h][s] << ";" << t << ";" << d << ";" << h << ";" << s << endl;
        }
      }
    }
  }

}



/* OVERWRITE FITNESS ARRAY FROM PREVIOUS ITERATION */
void ReplaceFit()
{
  int t,d,h,s;

  maxfitdiff = 0.0;
  totfitdiff = 0.0;

  for (t=1;t<maxT;t++)
  {
    for (d=0;d<maxD;d++)
    {
      for (h=0;h<maxH;h++)
      {
        for (s=0;s<maxS;s++)
        {
          totfitdiff = totfitdiff + abs(Wnext[t][d][h][s]-W[t][d][h][s]);
	  maxfitdiff = max(maxfitdiff, abs(Wnext[t][d][h][s]-W[t][d][h][s]));
          Wnext[t][d][h][s] = W[t][d][h][s];
        }
      }
    }
  }

}

// comparing hormone values
void ReplaceHormone()
{
  int t,d,s;

  maxhormonediff = 0.0;
  tothormonediff = 0.0;

  for (t=1;t<maxT;t++)
  {
    for (d=0;d<maxD;d++)
      {
        for (s=0;s<maxS;s++)
        {
          tothormonediff = tothormonediff + abs(hormone_next[t][d][s]-hormone[t][d][s]);
	  maxhormonediff = max(maxhormonediff, static_cast<double>(abs(hormone_next[t][d][s]-hormone[t][d][s])));
          hormone_next[t][d][s] = hormone[t][d][s];
        }
      }
    }
  }

/* PRINT OUT OPTIMAL STRATEGY */
void PrintStrat()
{
  int t,d,s;

  outputfile << "t" << "\t" << "d" << "\t" << "s" << "\t" << "hormone" << endl;

  for (t=0;t<maxT;t++)
  {
    for (d=0;d<maxD;d++)
    {
      for (s=0;s<maxS;s++)
      {
        outputfile << t << "\t" << d << "\t" << s << "\t" << hormone[t][d][s] << endl;
      }
    }
  }
  outputfile << endl;
  outputfile << "nIterations" << "\t" << i << endl;
  outputfile << endl;
}




/* WRITE PARAMETER SETTINGS TO OUTPUT FILE */
void PrintParams()
{
  outputfile << endl << "PARAMETER VALUES" << endl
       << "lambdaL: " << "\t" << lambdaL << endl
       << "lambdaA: " << "\t" << lambdaA << endl
       << "pAtt: " << "\t" << pAtt << endl
       << "alpha: " << "\t" << alpha << endl
       << "beta: " << "\t" << beta_b << endl
       << "mu: " << "\t" << mu << endl
       << "rho: " << "\t" << rho << endl
       << "h0: " << "\t" << h0 << endl
       << "omega: " << "\t" << omega << endl
       << "gamma: " << "\t" << gamma_g << endl
       << "maxI: " << "\t" << maxI << endl
       << "maxT: " << "\t" << maxT << endl
       << "maxD: " << "\t" << maxD << endl
       << "maxH: " << "\t" << maxH << endl
       << "maxS: " << "\t" << maxS << endl;
}



void SimAcutePhases(const string &base_name) // Simulating predator attack at t=10
{
    // 1) Simulation parameters
    const int simTime  = 50;      // We'll simulate from time=0 to time=50
    const int N        = 1000;    // Total individuals
    // We'll assume N is evenly divisible by maxS:
    int nPerPhase      = N / maxS; // # individuals starting in each phase

    // 2) Create arrays to track sums, sums of squares, and counts
    //    We store them by time (0..simTime) and breeding phase (0..maxS-1).
    static double sumD[51][101];       // sumD[time][s]
    static double sumsqD[51][101];     // sum of (damage^2)
    static double sumH[51][101];       // sum of hormone (or proportion)
    static double sumsqH[51][101];     // sum of (hormone^2)
    static int    countInd[51][101];   // how many individuals are in phase s at time t

    // We also track if the predator attacks at a given time:
    static bool attack[51];

    // 3) Initialize all arrays to zero/false
    for(int t=0; t<=simTime; t++)
    {
        attack[t] = false;
        for(int s=0; s<maxS; s++)
        {
            sumD[t][s]    = 0.0;
            sumsqD[t][s]  = 0.0;
            sumH[t][s]    = 0.0;
            sumsqH[t][s]  = 0.0;
            countInd[t][s] = 0;
        }
    }

    // Let's define the predator to attack exactly once at time=10
    if(simTime >= 10)
    {
        attack[10] = true; // Predator attacks at t=10
    }

    // 4) Output file to store results
    string fname = "SimAcutePhases_" + base_name + ".txt";
    ofstream outFile(fname.c_str());
    if(!outFile)
    {
        cerr << "Error opening file: " << fname << endl;
        return;
    }

    // Write a header
    outFile << "ACUTE ATTACK SIMULATION (split pop by breeding phase)\n";
    outFile << "time\ts\tnInd\tmeanD\tsdD\tmeanH\tsdH\n";

    // 5) Simulate each individual
    //    We distribute individuals so that nPerPhase start in phase=0,
    //    another nPerPhase in phase=1, etc.
    uniform_real_distribution<double> Uniform(0.0, 1.0);

    // We'll loop over each phase sVal as the starting phase
    for(int sVal=0; sVal<maxS; sVal++)
    {
        // We'll create nPerPhase individuals that start in sVal
        for(int i=0; i<nPerPhase; i++)
        {
            // (a) Initialize state for this individual
            int tState = maxT - 1; // time-since-attack state for DP
            int d      = 0;       // damage
            int sInd   = sVal;    // breeding phase for this individual
            int time   = 0;       // simulation time steps (0..simTime)

            // (b) Step through time
            while(time <= simTime)
            {
                // If the predator attacks at time=10, we reset tState=0
                // else we increment tState
                if(time == 10)
                {
                    tState = 0;
                }
                else
                {
                    tState = min(tState + 1, maxT - 1);
                }

                // (c) Retrieve optimal hormone from DP
                int h = hormone[tState][d][sInd];

                // (d) Update damage according to damage_new
                double d_next = damage_new[d][h];
                int    d1     = floor(d_next);
                int    d2     = ceil(d_next);
                double frac   = d_next - double(d1);

                double randVal = Uniform(mt);
                if(randVal < frac) d = d2; else d = d1;

                // (e) Record stats: we store the individual's damage/hormone
                //     in the bin for (time, sInd).
                sumD[time][sInd]    += d;
                sumsqD[time][sInd]  += (double(d) * double(d));

                double hProp = double(h) / double(maxH);
                sumH[time][sInd]    += hProp;
                sumsqH[time][sInd]  += (hProp * hProp);

                // increment count of individuals in that phase
                countInd[time][sInd]++;

                // (f) Advance the breeding phase for next time step
                sInd = (sInd + 1) % maxS;

                // (g) Move forward in simulation time
                time++;
            }
        }
    }

    // 6) Now compute means and standard deviations
    //    and output them to the file for each time and breeding phase.
    for(int t=0; t <= simTime; t++)
    {
        for(int s=0; s < maxS; s++)
        {
            int n = countInd[t][s];
            if(n > 0)
            {
                double meanD_ = sumD[t][s] / double(n);
                double varD_  = (sumsqD[t][s]/double(n)) - (meanD_*meanD_);

                double meanH_ = sumH[t][s] / double(n);
                double varH_  = (sumsqH[t][s]/double(n)) - (meanH_*meanH_);

                double sdD = (varD_>0.0) ? sqrt(varD_) : 0.0;
                double sdH = (varH_>0.0) ? sqrt(varH_) : 0.0;

                // We'll print: time, breeding phase, # of individuals, meanD, sdD, meanH, sdH
                outFile << t << "\t" 
                        << s << "\t" 
                        << n << "\t" 
                        << meanD_ << "\t" 
                        << sdD << "\t"
                        << meanH_ << "\t"
                        << sdH << "\n";
            }
            else
            {
                // No individuals in this bin => just print zero or skip
                outFile << t << "\t" 
                        << s << "\t" 
                        << 0 << "\t" 
                        << 0.0 << "\t" 
                        << 0.0 << "\t"
                        << 0.0 << "\t"
                        << 0.0 << "\n";
            }
        }
    }

    outFile.close();
    cout << "Simulation complete! Results in: " << fname << endl;
}






/* MAIN PROGRAM */
int main()
{
	c = 0;

		///////////////////////////////////////////////////////
		outfile.str("");
		outfile << "stress.txt";
		string outputfilename = outfile.str();
		outputfile.open(outputfilename.c_str());
		///////////////////////////////////////////////////////

        outputfile << "Random seed: " << seed << endl; // write seed to output file

        FinalFit();
        Predator();
        Death();
        Reproduction();
        Damage();

        cout << "i" << "\t" << "totfitdiff" << "\t" << "maxfitdiff" << "\t" << "tothormonediff" << "\t" << "maxhormonediff" << "\t" << "c" << endl;
        for (i=1;i<=maxI;i++)
          {
          OptDec();
          ReplaceFit();
	  ReplaceHormone();

          if (maxfitdiff < 0.01) 
            {
               cout << "Converged at iteration: " << i << ", maxfitdiff: " << maxfitdiff << endl;
               break; // strategy has converged on optimal solution, so exit loop
            }
          if (i==maxI) 
            { 
              outputfile << "*** DID NOT CONVERGE WITHIN " << i << " ITERATIONS ***" << endl;
            }
	  if (maxhormonediff==0)
	    {
	      c++;  // Increment counter if maxfitdiff is still 0	    
	    }
	  else
	   {
	     c = 0;  
	   }
          if (c==150)
	   {
               cout << "Converged at iteration: " << i << ", maxfitdiff: " << maxfitdiff << endl;
               break; // strategy has converged on optimal solution, so exit loop
            }
 	  if (i%skip==0)
            {
              cout << i << "\t" << totfitdiff << "\t" << maxfitdiff << "\t" << tothormonediff << "\t" << maxhormonediff << "\t" << c << endl; // show fitness difference every 'skip' generations
            }
          }

        cout << endl;
        outputfile << endl;

        PrintStrat();
        PrintParams();
        outputfile.close(); //done with DP results, optimal strategy found

	 // 5) *** Now run your acute-attack simulation *** 
  	 // e.g. if you named it SimAcutePhases()
    	 SimAcutePhases("ShortCycle");

  return 0;
}
