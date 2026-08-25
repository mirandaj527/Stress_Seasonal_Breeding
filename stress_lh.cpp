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

double lambdaA  = 0.035;   // probability that predator arrives
double lambdaL  = 0.065;   // probability that predator leaves
const double pAtt     = 0.5;     // probability that predator attacks if present
double alpha    = 0.1;     // parameter controlling effect of hormone level on pKill
const double beta_b   = 0.0;     // parameter controlling effect of hormone level on reproductive rate
const double kappa    = 0.0;     // Parameter controlling affect of damage on mortality
double mu       = 0.02;   // background mortality (independent of hormone level and predation risk)
double rho      = 1.0;    // Fixed rate of repair
const double h0       = 20.0;   // Reference hormone level
double omega    = 0.1;   // Effect of deviations from h0 on damage build-up
double gamma_g  = 0.1;     // Effect of damage on reproductive output
double fec_scale = 1.0; 
double fec_baseline = 0.0;

int maxI        = 100000; // maximum number of iterations
int maxT        = 25;     // maximum number of time steps since last saw predator
int maxD        = 100;    // Number of discrete damage levels?
int maxH        = 100;    // maximum hormone level
int maxS        = 10;       // Length of the breeding cycle
int skip        = 10;       // interval between print-outs

bool stdoutput = true;

// Create a random engine with your chosen seed
std::random_device rd{}; 
unsigned seed = rd();
static std::mt19937 mt(seed);

std::ofstream outputfile; // output file
std::stringstream outfile; // for naming output file

//int hormone[maxT][maxD][maxS];          // hormone level (strategy)                               //int hormone_next[maxT][maxD][maxS];   // optimal hormone levels given next time step

std::vector<                                   // vector 1: maxT
  std::vector<                                 // vector 2: maxD                                 
    std::vector<int>>> hormone, hormone_next;                        // vector 3: maxS

//double pKill[maxH];                     // probability of being killed by an attacking predator
std::vector<double> pKill;

//double repro[maxD][maxH][maxS];         // reproductive output
std::vector<                                   // std::vector 1: maxD
  std::vector<                                 // std::vector 2: maxH
    std::vector<double>>>                      // std::vector 3: maxS
      repro;
//double Wopt[maxT][maxD][maxS];          // fitness immediately after predator has/hasn't attacked, under optimal decision h
std::vector<                                   // std::vector 1: maxT
  std::vector<                                 // std::vector 2: maxD                                 
    std::vector<double>>>                         // std::vector 3: maxS
      Wopt;

//double W[maxT][maxD][maxH][maxS];       // expected fitness at start of time step, before predator does/doesn't attack
std::vector<                                   // std::vector 1: maxT
  std::vector<                                 // std::vector 2: maxD
    std::vector<                               // std::vector 3: maxH
      std::vector<double>>>>                   // std::vector 4: maxS
        W;

//double Wnext[maxT][maxD][maxH][maxS];   // expected fitness at start of next time step
std::vector<                                   // std::vector 1: maxT
  std::vector<                                 // std::vector 2: maxD
    std::vector<                               // std::vector 3: maxH
      std::vector<double>>>>                   // std::vector 4: maxS
        Wnext;

//double pPred[maxT];                     // probability that predator is present
std::vector<double> pPred;
//double damage_new[maxD][maxH];          // array of current damage
std::vector<                                   // std::vector 1: maxD
  std::vector<double>>                         // std::vector 2: maxH
    damage_new(maxD,
      std::vector<double> (maxH, 0.0)
);
//double background_mortality[d]
std::vector<double> background_mortality(maxD, 0.0);

double totfitdiff;                        // fitness difference between optimal strategy in successive iterations
double maxfitdiff;
double tothormonediff;
double maxhormonediff;

int i;     // iteration

int c;     // count

/* initialize all vectors after we have
 * obtained dimensions via command line */
void InitVectors()
{
    // initialize hormone vector
    // which has dimensions
    //int hormone[maxT][maxD][maxS];          
    //
    // we do so by creating local vec
    // and assigning that to global
    std::vector <  // vector 1: maxT
        std::vector <  // vector 2: maxD
            std::vector < int >>> // vector 3: maxS
            local_hormone(maxT,
                    std::vector<                                
                      std::vector<int>> 
                        (maxD, 
                          std::vector<int> (maxS, 0)     // initialising internal vector of length maxS to contain 0s
                        )
            );

    hormone = local_hormone;
    hormone_next = local_hormone;

    std::vector <double> pKill_local(maxH, 0.0);
    pKill = pKill_local;

    //double repro[maxD][maxH][maxS];         // reproductive output
    std::vector < // maxD
        std::vector < // maxH
            std::vector < double >>> // maxS
                local_repro(
                        maxD,
                        std::vector <
                            std::vector < double >>
                                (maxH,
                                    std::vector < double > (maxS, 0.0)
                                )
                            );

    repro = local_repro;

    //double Wopt[maxT][maxD][maxS];          // fitness immediately after predator has/hasn't attacked, under optimal decision h
    std::vector<                                   // vector 1: maxT
        std::vector<                                 // vector 2: maxD                                          
            std::vector<double>>>                         // vector 3: maxS
                Wopt_local(maxT,                                        
                    std::vector<                                
                        std::vector<double>> 
                            (maxD, 
                              std::vector<double> (maxS)       
                            )
                        );

    Wopt = Wopt_local;

//double W[maxT][maxD][maxH][maxS];       // expected fitness at start of time step, before predator does/doesn't attack
    std::vector<                                   // vector 1: maxT
        std::vector<                                 // vector 2: maxD
            std::vector<                               // vector 3: maxH
                std::vector<double>>>>                   // vector 4: maxS
                    Wlocal(maxT,
                        std::vector<
                            std::vector<
                                std::vector<double>>>
                                    (maxD,
                                        std::vector<
                                            std::vector<double>>
                                                (maxH,
                                                    std::vector<double> (maxS)
                                                )
                                    )
                            );
    W = Wlocal;
//double Wnext[maxT][maxD][maxH][maxS];   // expected fitness at start of next time step

    Wnext = Wlocal;


    // initialize predation probabilities
    std::vector <double> pPred_local(maxT, 0.0);
    pPred = pPred_local;

    // initialize damage updating vector


    std::vector<                                   // vector 1: maxD
        std::vector<double>>                         // vector 2: maxH
            damage_new_local(maxD,
                std::vector<double> (
                    maxH, 0.0)
                            );

    damage_new = damage_new_local;

    std::vector<double> background_mortality_local(maxD, 0.0);

    background_mortality = background_mortality_local;

} // end InitVectors()

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
      repro[d][h][0] = fec_baseline + fec_scale * exp(-(beta_b*(static_cast<double>(h)/static_cast<double>(maxH)) 
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
      damage_new[d][h] = std::max(
              0.0, 
              std::min(static_cast<double>(maxD-1), damage_new[d][h])); 
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
        background_mortality[d] = std::min(
                1.0, 
                (mu + (1.0 - std::exp(-kappa*(static_cast<double>(d)/static_cast<double>(maxD))))));
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
            d1=std::floor(damage_new[d][h]); // lower integer value of damage
            d2=std::ceil(damage_new[d][h]); // top integer value
            ddif = damage_new[d][h] - static_cast<double>(d1); // calculate difference 

            // calculate fitness from W' as a function of t, d, s, h
            fitness = 
                (1.0 - ddif) * Wnext[std::min(maxT - 1, t + 1)][d1][h][(s + 1) % maxS] // deterministic rounding
                + 
                ddif * Wnext[std::min(maxT - 1, t + 1)][d2][h][(s + 1) % maxS];
	      

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
          // std::cout << W[t][d][h][s] << ";" << t << ";" << d << ";" << h << ";" << s << std::endl;
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
          totfitdiff = totfitdiff + std::abs(Wnext[t][d][h][s]-W[t][d][h][s]);
	  maxfitdiff = std::max(maxfitdiff, std::abs(Wnext[t][d][h][s]-W[t][d][h][s]));
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
          tothormonediff = tothormonediff + std::abs(hormone_next[t][d][s]-hormone[t][d][s]);
	  maxhormonediff = std::max(maxhormonediff, static_cast<double>(std::abs(hormone_next[t][d][s]-hormone[t][d][s])));
          hormone_next[t][d][s] = hormone[t][d][s];
        }
      }
    }
  }

/* PRINT OUT OPTIMAL STRATEGY */
void PrintStrat()
{
  int t,d,s;

  outputfile << "t" << "\t" << "d" << "\t" << "s" << "\t" << "hormone" << std::endl;

  for (t=0;t<maxT;t++)
  {
    for (d=0;d<maxD;d++)
    {
      for (s=0;s<maxS;s++)
      {
        outputfile << t << "\t" << d << "\t" << s << "\t" << hormone[t][d][s] << std::endl;
      }
    }
  }
  outputfile << std::endl;
  outputfile << "nIterations" << "\t" << i << std::endl;
  outputfile << std::endl;
}




/* WRITE PARAMETER SETTINGS TO OUTPUT FILE */
void PrintParams()
{
  outputfile << std::endl << "PARAMETER VALUES" << std::endl
       << "lambdaL: " << "\t" << lambdaL << std::endl
       << "lambdaA: " << "\t" << lambdaA << std::endl
       << "pAtt: " << "\t" << pAtt << std::endl
       << "alpha: " << "\t" << alpha << std::endl
       << "beta: " << "\t" << beta_b << std::endl
       << "mu: " << "\t" << mu << std::endl
       << "rho: " << "\t" << rho << std::endl
       << "fec_scale: " << "\t" << fec_scale << std::endl
       << "fec_baseline: " << "\t" << fec_baseline << std::endl
       << "h0: " << "\t" << h0 << std::endl
       << "omega: " << "\t" << omega << std::endl
       << "gamma: " << "\t" << gamma_g << std::endl
       << "maxI: " << "\t" << maxI << std::endl
       << "maxT: " << "\t" << maxT << std::endl
       << "maxD: " << "\t" << maxD << std::endl
       << "maxH: " << "\t" << maxH << std::endl
       << "maxS: " << "\t" << maxS << std::endl;
}



void SimAcutePhases(const std::string &base_name) // Simulating predator attack at t=10
{
    // 1) Simulation parameters
    const int simTime  = 50;      // We'll simulate from time=0 to time=50
    const int N        = 1000;    // Total individuals
    // We'll assume N is evenly divisible by maxS:
    int nPerPhase      = N / maxS; // # individuals starting in each phase

    // 2) Create arrays to track sums, sums of squares, and counts
    //    We store them by time (0..simTime) and breeding phase (0..maxS-1).
    std::vector < 
        std::vector < double >> sumD(simTime + 1, 
                                    std::vector< double > (maxS, 0.0)
                );

    std::vector < 
        std::vector < double >> sumsqD(simTime + 1, 
                                    std::vector< double > (maxS, 0.0)
                );
    
    std::vector < 
        std::vector < double >> sumH(simTime + 1, 
                                    std::vector< double > (maxS, 0.0)
                );
    
    std::vector < 
        std::vector < double >> sumsqH(simTime + 1, 
                                    std::vector< double > (maxS, 0.0)
                );
    
    std::vector < 
        std::vector < int >> countInd(simTime + 1, 
                                    std::vector< int > (maxS, 0.0)
                );

    // 3) Output file to store results
    std::string fname = "SimAttack_" + base_name + ".txt";
    std::ofstream outFile(fname.c_str());
    if(!outFile)
    {
        std::cerr << "Error opening file: " << fname << std::endl;
        return;
    }

    // Write a header
    outFile << "ACUTE ATTACK SIMULATION (split pop by breeding phase)\n";
    outFile << "time\ts\tnInd\tmeanD\tsdD\tmeanH\tsdH\n";

    // 4) Simulate each individual
    //    We distribute individuals so that nPerPhase start in phase=0,
    //    another nPerPhase in phase=1, etc.
    std::uniform_real_distribution<double> Uniform(0.0, 1.0);

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
                    tState = std::min(tState + 1, maxT - 1);
                }

                // (c) Retrieve optimal hormone from DP
                int h = hormone[tState][d][sInd];

                // (d) Update damage according to damage_new
                double d_next = damage_new[d][h];
                int    d1     = std::floor(d_next);
                int    d2     = std::ceil(d_next);
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

    if (stdoutput)
    {
        std::cout << "Simulation complete! Results in: " << fname << std::endl;
    }
}


/* MAIN PROGRAM */
int main(int argc, char* argv[])
{
    // ----------------------------------------------------
    // 1) Provide default filenames & parameters
    //    for DP output and simulation base name
    // ----------------------------------------------------
    std::string dpOutputFilename = "stress.txt";  
    std::string simOutputBase    = "Test";

    // ----------------------------------------------------
    // 2) Parse command-line arguments
    //    e.g. "./Program lambdaA=0.04 dpFile=dpOut.txt simBase=MySim"
    // ----------------------------------------------------
    for (int argIndex = 1; argIndex < argc; argIndex++)
    {
        std::string arg = argv[argIndex];

        if      (arg.rfind("lambdaA=", 0) == 0)  { lambdaA   = std::stod(arg.substr(8)); }
        else if (arg.rfind("lambdaL=", 0) == 0)  { lambdaL   = std::stod(arg.substr(8)); }
        else if (arg.rfind("alpha=",   0) == 0)  { alpha     = std::stod(arg.substr(6)); }
        else if (arg.rfind("mu=",      0) == 0)  { mu        = std::stod(arg.substr(3)); }
        else if (arg.rfind("rho=",     0) == 0)  { rho       = std::stod(arg.substr(4)); }
        else if (arg.rfind("omega=",   0) == 0)  { omega     = std::stod(arg.substr(6)); }
        else if (arg.rfind("gamma_g=", 0) == 0)  { gamma_g   = std::stod(arg.substr(8)); }
        else if (arg.rfind("stdoutput=", 0) == 0)  { stdoutput = static_cast<bool>(std::stoi(arg.substr(10))); }
        else if (arg.rfind("maxS=", 0) == 0)  { maxS = std::stoi(arg.substr(5)); }
        // else if (arg.rfind("maxS=",    0) == 0)  { maxS      = std::stoi(arg.substr(5)); }
        else if (arg.rfind("dpFile=",  0) == 0)  { dpOutputFilename = arg.substr(7); }
        else if (arg.rfind("simBase=", 0) == 0)  { simOutputBase    = arg.substr(8); }
        else {
            std::cerr << "Unrecognized argument: " << arg << std::endl;
        }
    }

    // ----------------------------------------------------
    // 3) Print out final parameter values 
    //    so we know what's being used
    // ----------------------------------------------------
    
    if (stdoutput)
    {
        std::cout << "Using parameters:\n";
        std::cout << "  lambdaA = " << lambdaA << "\n"
                  << "  lambdaL = " << lambdaL << "\n"
                  << "  alpha   = " << alpha << "\n"
                  << "  mu      = " << mu << "\n"
                  << "  rho     = " << rho << "\n"
                  << "  omega   = " << omega << "\n"
                  << "  gamma_g = " << gamma_g << "\n"
                  << "  maxS    = " << maxS << "\n\n";
        std::cout << "DP output file: " << dpOutputFilename << "\n";
        std::cout << "Simulation base: " << simOutputBase << "\n\n";
    }

    InitVectors();
    // ----------------------------------------------------
    // DP calculations to find optimal strat
    // ----------------------------------------------------
    c = 0;

    // For the DP output file:
    outfile.str("");
    outfile << dpOutputFilename;
    std::string outputfilename = outfile.str();
    outputfile.open(outputfilename.c_str());

    // Write seed
    outputfile << "Random seed: " << seed << std::endl;

    // Model init
    FinalFit();
    Predator();
    Death();
    Reproduction();
    Damage();

    // DP iteration
    if (stdoutput)
    {
        std::cout << "i" << "\t" << "totfitdiff" << "\t" << "maxfitdiff" << "\t"
             << "tothormonediff" << "\t" << "maxhormonediff" << "\t" << "c" << std::endl;
    }

    for (i=1; i <= maxI; i++)
    {
        OptDec();
        ReplaceFit();
        ReplaceHormone();

        if (maxfitdiff < 0.000001) 
        {
            if (stdoutput)
            {
               std::cout << "Converged at iteration: " << i 
                    << ", maxfitdiff: " << maxfitdiff << std::endl;
            }
           break; // strategy has converged, so exit loop
        }
        if (i == maxI) 
        { 
           outputfile << "*** DID NOT CONVERGE WITHIN " << i << " ITERATIONS ***" << std::endl;
        }
        if (maxhormonediff == 0)
        {
           c++;  // Increment counter if maxhormonediff is still 0
        }
        else
        {
           c = 0;
        }

        if (c >= 150)
        {
            if (stdoutput)
            {
               std::cout << "Converged at iteration: " << i 
                    << ", maxfitdiff: " << maxfitdiff << std::endl;
            }
           break; 
        }
        if (i % skip == 0)
        {
            if (stdoutput)
            {
               std::cout << i << "\t" << totfitdiff << "\t" << maxfitdiff << "\t"
                    << tothormonediff << "\t" << maxhormonediff << "\t" << c << std::endl;
            }
        }
    }

    if (stdoutput)
    {
        // Convergence done or maxI reached
        std::cout << std::endl;
    }
   
    outputfile << std::endl;

    // Print strategy & params
    PrintStrat();
    PrintParams();
    outputfile.close();

    // ----------------------------------------------------
    // 5) Finally, run the simulation 
    //    passing simOutputBase to the function
    // ----------------------------------------------------
    SimAcutePhases(simOutputBase);

    return 0;
}
