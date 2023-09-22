/*==========================================================
 * This file conducts the contraction mapping to obtain contract-month specific mean 
 * utilities. It takes into account consumer heterogeneity in fallback choice sets via 
 * advertising, endogenous PCW usage decisions, and switching costs.
 *
 * Inputs:
 *      tol = scalar indicating the contraction mapping tolerance
 *      maxiter = scalar indicating maximimum number of iterations
 *      expmvalold = vector of exponentiated mean utilities.
 *      expmu = matrix of individual specific (exponentiated) utility deviations
 *      cdindex = index vector as in Nevo's BLP-code that indicates last observation of each market.
 *      shares = vector of contract-month level market shares.
 *      shares_init = vector of type-specific contract-level market shares at beginning fo initial period.
 *      aware_adv = indicator vector for which contracts each consumer is informed about via advertising or previous choice.
 *      expPBeliefDiff = simulated price beliefs for products that a consumer is not perfectly informed about (already exponentiated and weighted with consumer-specific price coefficient)
 *      expEpsShocks = (exponentiated) simulated logit shocks for each consumer-market-contract observation (used for evaluating expected benefit from using PCW)
 *      exp_psi = (exponentiated) vector of consumer-specific switching cost parameters
 *      exp_kappa = vector of consumer-specific PCW usage costs.
 *      contrMap_dummy = binary indicator whether contraction mapping is run or whether we just solve for model predictions; 1 = solve contraction mapping, 0 = only model predictions for given value of mean utilities.
 *
 * Output:
 *      expmval = vector of updated (exponentiated) mean utilities.
 *      PCWUsagei = consumer-market specific vector of PCW usage probability
 *      pSharesi = predicted consumer-market specific vecotr of contract-level market shares
 *      pSharesiPCW =  predicted consumer-market specific vecotr of contract-level market shares conditional on using PCW
 *      pSharesiNoPCW =  predicted consumer-market specific vecotr of contract-level market shares conditional on not using PCW
 *      pChurni = predicted consumer-market specific churn rates 
 *      cSurplusPCWi = implied consumer-specific consumer surplus conditional on using PCW
 *      cSurplusNoPCWi = implied consumer-specific consumer surplus conditional on not using PCW
 *      ccpMat = predicted CCP matrices for each simulated consumer. 
 *      mktnorm = convergence criterion after concluding contraction mapping.   
 *      mktiter = number of iterations needed until convergence.
 *
 * Written by Stefan Weiergraeber, May 2019
 * Updated by Stefan Weiergraeber, August 2021
 *========================================================*/
#include "mex.h"
#include "math.h"
#include "stdio.h"

/*==========================================================
 * The printArray function makes it simpler to print arrays,
 * thereby faciliating the debugging process. 
 *========================================================*/
void printArray(int size, char *name, double *array)
{
    int j;
    for (j = 0; j < size; j++)
    {
        mexPrintf("%s[%d]=%f\n", name, j, array[j]);
    }
}

/*==========================================================
 * The contrMap function contains the computational routine.
 *========================================================*/
void contrMap(double tol, double maxiter,
              double *expmvalold, double *expmu,
              mwSize nJ, mwSize nI, mwSize nM, mwSize nsP,
              double *cdindex, double *shares, double *shares_init, double *aware_adv,
              double *expPBeliefDiff, double *expEpsShocks,
              double *exp_psi, double *exp_kappa, _Bool contrMap_dummy,
              double *mktnorm, double *mktiter,
              double *expmval, double *PCWUsagei, double *pSharesi, double *pSharesiPCW, double *pSharesiNoPCW,
              double *pChurni, double *cSurplusPCWi, double *cSurplusNoPCWi, double *ccpMat)
{

    const double LogitSmoother = 0.15;        // factor to smooth binary choice probabilities.
    const double invLS = 1.0 / LogitSmoother; // inverse of logit smoother above.
   
    /* Assigning variable names */
    mwSize j;      // product index
    mwSize i;      // consumer index
    mwSize mkt;    // market (i.e., month) indicator
    mwSize ii;     // loop index for simulated consumers
    mwSize k;      // loop index for current product
    mwSize k_lag;  // loop index for previous product
    mwSize l;      // loop index over simulated price beliefs
    mwSize nJM;    // number of market-specific products
    mwSize mStart; // Starting index of current month
    mwSize mEnd;   // End index of current month

    double invnI;  // inverse of number of consumers: 1/NS
    double invnsP; // inverse of number of simulated prices and epsilon shocks.
    // Inverse of number of consumers and simulation draws to compute averages over these dimensions.
    invnI = 1.0 / (double)nI;
    invnsP = 1.0 / (double)nsP;

    // Variables needed to measure convergence.
    double maxabsdev;
    double absdev;
    double norm;
    mwSize iter;
    // Some debugger variables.
    double sum_ind_shares;
    double sum_agg_shares;
    // Aux variable to hold CCP number.
    double ccpAux;

    // Declaring vector of market-specific search and switching costs parameters.
    mxArray *expPsiMx;   /* memory for market-specific vector of switching costs */
    mxArray *expKappaMx; /* memory for market-specific vector of PCW search costs */
    double *expPsiM;     /* market-specific vector of switching cost parameters */
    double *expKappaM;   /* market-specific vector of PCW search cost parameters */
    expPsiMx = mxCreateDoubleMatrix((mwSize)nI, 1, mxREAL);
    expPsiM = mxGetDoubles(expPsiMx);
    expKappaMx = mxCreateDoubleMatrix((mwSize)nI, 1, mxREAL);
    expKappaM = mxGetDoubles(expKappaMx);

    // Declaring a bunch of required vectors containing various utility components.
    // These are all pointers to one-dimensional arrays.
    mxArray *expmvaloldMx;
    double *expmvaloldM;
    mxArray *sharesMx;
    double *sharesM;
    mxArray *sharesiLagMx;
    double *sharesiLagM;
    mxArray *expmuMx;
    double *expmuM;

    // Declare market specific CCP matrix (requested by referee)
    mxArray *ccpiMx; // individual-specific CCPs for one market.
    double *ccpiM;
    mxArray *incValFBMx; // Denominator of choice probabilities for PCW non-users.
    double *incValFBM;
    mxArray *incValFAMx; // Denominator of choice probabilities for PCW users after having used PCW (choice under full information)
    double *incValFAM;
    mxArray *EUNoSearchMx; // Expected utility from not searching PCW.
    double *EUNoSearchM;
    mxArray *EUSearchNetMx; // Expected utility from searching PCW taking into account PCW search cost (=net).
    double *EUSearchNetM;
    mxArray *PCWUsageMx; // Market-specific PCW usage vector.
    double *PCWUsageM;
    mxArray *exputilMx; // market-specific exp(utility) as function of consumer type, previous contract, current contract
    double *exputilM;
    double exputilBelief; // aux variable for expected exp(utility) from using PCW (uses price beliefs instead of actual prices)

    // These objects are needed to simulate PCW usage decision.
    // Keep in mind that the following block is all in exponentiated utilities.
    // Simulated utilities for one consumer and lagged choice.
    // Dimensions: J x NS_p
    mxArray *USimNoPCWx; // Simulated utilities when not using PCW for a specific consumer: for each product and simulation draw for price and epsilon shock.
    double *USimNoPCW;
    mxArray *USimPCWx; // Simulated utilities when using the PCW for a specific consumer: for each product and simulation draw for price and epsilon shock.
    double *USimPCW;
    mxArray *UMaxSimNoPCWx; // simulated maximum utility for a given simulation draw (max over availabe contracts)
    double *UMaxSimNoPCW;
    mxArray *UMaxSimPCWx; // simulated maximum utility for a given simulation draw (max over availabe contracts)
    double *UMaxSimPCW;
    // KEEP THIS AS POINTERS FOR NOW BECAUSE THEY NEED TO OPERATE ON POINTER OBJECT!
    // Not sure this needs to be the case. Working, but could be made more elegant.
    mxArray *EUMaxNoPCWx; // expected maximum exp(utility) when NOT using PCW (scalar for given consumer and previous choice)
    double *EUMaxNoPCW;
    mxArray *EUMaxPCWx; // expected maximum exp(utility) when using PCW (scalar for given consumer and previous choice)
    double *EUMaxPCW;
    double PCWSmoothAux; // store aux object for smoothed PCW usage probability.
    mxArray *pSharesiMx; // predicted shares for market m (stacked for consumer types and current contracts)
    double *pSharesiM;
    // Analog conditional on using PCW and not using PCW.
    mxArray *pSharesiPCWMx;
    double *pSharesiPCWM;
    mxArray *pSharesiNoPCWMx;
    double *pSharesiNoPCWM;
    mxArray *pSharesMx; // aggregated market share predictions (averaged over consumers)
    double *pSharesM;
    mxArray *pChurniMx; // analogous prediction for market-specific churn rates for each simulated consumer.
    double *pChurniM;
    mxArray *awareiAdvMx; // market-specific indicator for consumer i subscribed to previous contract k_lag is aware of k  through advertising or previous subscription.
    double *awareiAdvM;
    mxArray *expPBeliefDiffMx; // aux object that captures difference between actual price and expected / simulated price. this is used to adjust realized utility computed by the model into utility belief.
    double *expPBeliefDiffM;
    mxArray *expEpsShocksMx; // exp(epsilon): logit shocks to utility are needed to compute expected utility from using PCW, but not for computation of CCPs.
    double *expEpsShocksM;
    mxArray *expmvalMx; // updated exp(delta) for market m
    double *expmvalM;

    ////////////////////////////////////////////////////////////////////////////
    // Declare variables for slotting output objects.
    ////////////////////////////////////////////////////////////////////////////
    // Detailed PCW usage vector for each consumer type, previous contract choice and market.
    mxArray *PCWUsagex;
    double *PCWUsage;
    PCWUsagex = mxCreateDoubleMatrix(nJ * nI, 1, mxREAL);
    PCWUsage = mxGetDoubles(PCWUsagex);
    // Averaged PCW usage vector over consumer types and previous choices.
    // This should me measured against the observed vector of PCW usage.
    mxArray *PCWUsageAvgx; // Total average PCW usage vector on consumer type level.
    double *PCWUsageAvg;
    PCWUsageAvgx = mxCreateDoubleMatrix(nM, 1, mxREAL);
    PCWUsageAvg = mxGetDoubles(PCWUsageAvgx);
    // Surplus depends on previous choice and consumer type and market.
    mxArray *ConsumerSurplusx;
    double *ConsumerSurplus;
    ConsumerSurplusx = mxCreateDoubleMatrix(nJ * nI, 1, mxREAL);
    ConsumerSurplus = mxGetDoubles(ConsumerSurplusx);
    
    /*==========================================================
    * Conducting the contraction mapping.
    * In our model markets are months. Because of state dependence model has to
    * be solved recursively.
    *========================================================================*/

    for (mkt = 0; mkt < nM; mkt++)
    {
        iter = 0;
        norm = 1.0;

        // Start/stop variables to help create market-specific objects.
        // These are technically not needed since we have the same number of products each month.
        if (mkt == 0)
        {
            mStart = 0;
            nJM = cdindex[mkt];
        }
        else
        {
            mStart = cdindex[mkt - 1];
            nJM = cdindex[mkt] - cdindex[mkt - 1];
        }
        mEnd = cdindex[mkt];

        /* Create month-specific objects: everything specifically for market m */
        // Creating these could be moved outside of the loop in order to save some time.
        expmvaloldMx = mxCreateDoubleMatrix((mwSize)nJM, 1, mxREAL); // exponentiated mean utility
        expmvaloldM = mxGetDoubles(expmvaloldMx);
        expmuMx = mxCreateDoubleMatrix(nJM * nI, 1, mxREAL); // exponentiated consumer-specific deviation from mean utility
        expmuM = mxGetDoubles(expmuMx);
        exputilMx = mxCreateDoubleMatrix(nJM * nJM * nI, 1, mxREAL); // utility from current contracts conditional on consumer type and lagged choice
        exputilM = mxGetDoubles(exputilMx);
        incValFBMx = mxCreateDoubleMatrix((mwSize)nJM * nI, 1, mxREAL); // denominator for CCPs for PCW non-users
        incValFBM = mxGetDoubles(incValFBMx);
        incValFAMx = mxCreateDoubleMatrix((mwSize)nJM * nI, 1, mxREAL); // denominator for CCPs for PCW users
        incValFAM = mxGetDoubles(incValFAMx);

        sharesMx = mxCreateDoubleMatrix((mwSize)nJM, 1, mxREAL); // aggregate market share predictions
        sharesM = mxGetDoubles(sharesMx);
        sharesiLagMx = mxCreateDoubleMatrix((mwSize)nJM * nI, 1, mxREAL); // consumer-specific lagged market shares
        sharesiLagM = mxGetDoubles(sharesiLagMx);

        expPBeliefDiffMx = mxCreateDoubleMatrix((mwSize)nJM * nI * nsP, 1, mxREAL); // difference between simulated and realized price draw
        expPBeliefDiffM = mxGetDoubles(expPBeliefDiffMx);
        expEpsShocksMx = mxCreateDoubleMatrix((mwSize)nJM * nI * nsP, 1, mxREAL); // logit shocks for simulating PCW usage decision
        expEpsShocksM = mxGetDoubles(expEpsShocksMx);

        EUSearchNetMx = mxCreateDoubleMatrix(nJM * nI, 1, mxREAL); // expected utility from using PCW conditional on consumer type and lagged contract.
        EUSearchNetM = mxGetDoubles(EUSearchNetMx);
        EUNoSearchMx = mxCreateDoubleMatrix(nJM * nI, 1, mxREAL); // expected utility from not using PCW conditional on consumer type and lagged contract.
        EUNoSearchM = mxGetDoubles(EUNoSearchMx);

        PCWUsageMx = mxCreateDoubleMatrix(nJM * nI, 1, mxREAL); // PCW usage share for each consumer and lagged contract
        PCWUsageM = mxGetDoubles(PCWUsageMx);
        pSharesiPCWMx = mxCreateDoubleMatrix(nJM * nI, 1, mxREAL); // consumer-specific market shares among PCW users
        pSharesiPCWM = mxGetDoubles(pSharesiPCWMx);
        pSharesiNoPCWMx = mxCreateDoubleMatrix(nJM * nI, 1, mxREAL); // consumer-specific market shares among PCW non-users
        pSharesiNoPCWM = mxGetDoubles(pSharesiNoPCWMx);
        pChurniMx = mxCreateDoubleMatrix(nJM * nI, 1, mxREAL); // consumer-specific churn rates
        pChurniM = mxGetDoubles(pChurniMx);
        pSharesMx = mxCreateDoubleMatrix(nJM, 1, mxREAL);
        pSharesM = mxGetDoubles(pSharesMx);
        // Market-specific and consumer-specific CCP matrix.
        ccpiMx = mxCreateDoubleMatrix(nJM * nJM * nI, 1, mxREAL);
        ccpiM = mxGetDoubles(ccpiMx);
        // Market-specific aggregate CCP matrix.
        // Exogenous awareness indicator through advertising and awareness
        // because it's a consumer's previous product.
        // for each combination of: current produt-simulated consumer
        // Its structure is analogous to that of exputiliM.
        awareiAdvMx = mxCreateDoubleMatrix(nJM * nJM * nI, 1, mxREAL);
        awareiAdvM = mxGetDoubles(awareiAdvMx);

        expmvalMx = mxCreateDoubleMatrix(nJM, 1, mxREAL); // updated exp(delta)
        expmvalM = mxGetDoubles(expmvalMx);

        // Consumer-lagged choice specific matrix for simulated utilities.
        // Dimension: considered product x simulation draw
        USimNoPCWx = mxCreateDoubleMatrix((mwSize)nJM * nsP, 1, mxREAL);
        USimNoPCW = mxGetDoubles(USimNoPCWx);
        USimPCWx = mxCreateDoubleMatrix((mwSize)nJM * nsP, 1, mxREAL);
        USimPCW = mxGetDoubles(USimPCWx);

        // Maximum utilities over products for one consumer and lagged choice when using and not using PCW, respectively.
        // Dimensions: NS_p (number of simulation draws)
        UMaxSimNoPCWx = mxCreateDoubleMatrix((mwSize)nsP, 1, mxREAL);
        UMaxSimNoPCW = mxGetDoubles(UMaxSimNoPCWx);
        UMaxSimPCWx = mxCreateDoubleMatrix((mwSize)nsP, 1, mxREAL);
        UMaxSimPCW = mxGetDoubles(UMaxSimPCWx);

        // Expected max utility for one consumer and lagged choice.
        // Dimensions: array of dimension 1 x 1
        EUMaxNoPCWx = mxCreateDoubleMatrix(1, 1, mxREAL);
        EUMaxNoPCW = mxGetDoubles(EUMaxNoPCWx);
        EUMaxPCWx = mxCreateDoubleMatrix(1, 1, mxREAL);
        EUMaxPCW = mxGetDoubles(EUMaxPCWx);

        // Fill month-specific objects with data for consumer-specific switching costs and search costs.
        for (i = 0; i < nI; i++) // loop over consumers.
        {
            expKappaM[i] = exp_kappa[mkt * nI + i];
            expPsiM[i] = exp_psi[mkt * nI + i];
            // Fill exogenous awareness component for each consumer-previous product -> current product awareness.
            // Loop over previous products
            for (k_lag = 0; k_lag < nJM; k_lag++)
            {
                for (j = 0; j < nJM; j++)
                {
                    awareiAdvM[i * nJM * nJM + k_lag * nJM + j] = aware_adv[mkt * nJM * nJM * nI + i * nJM * nJM + nJM * k_lag + j];
                    
                } // end loop over current products.
            }     // end loop over previous products.
        }         // end loop over consumers

        // Fill month-specific objects.
        for (j = mStart; j < mEnd; j++) // loop over current choices.
        {
            expmvaloldM[j - mStart] = expmvalold[j]; // exponentiated mean utilities
            // mexPrintf("expmvaloldM[%d]: %f \n",j-mStart,expmvaloldM[j-mStart]);
            sharesM[j - mStart] = shares[j]; // observed market shares.

            // Individual-specific utility components (mu)
            for (i = 0; i < nI; i++) // loop over consumers.
            {
                expmuM[nJM * i + j - mStart] = expmu[nJ * i + j];
                // mexPrintf("expmuM[%d]: %f \n",nJM*i+j-mStart,expmuM[nJM*i+j-mStart]);
                if (mkt == 0)
                {
                    // Lagged market shares for first month.
                    // If mkt = 0 -> mStart = 0
                    sharesiLagM[nJM * i + j] = shares_init[nJM * i + j];
                   
                }
                else
                // For all months except first one:
                {
                    // Fill with market shares predictions from previous months.
                    sharesiLagM[nJM * i + j - mStart] = pSharesiM[nJM * i + j - mStart];
                    
                }
                // Fill month-specific price belief data.
                for (l = 0; l < nsP; l++) // loop over p belief and epsilon simulations.
                {
                    // Filling needs to be done correctly in order of: consumers i - current choices k - simulation draws l.
                    expPBeliefDiffM[i * nJM * nsP + (j - mStart) * nsP + l] = expPBeliefDiff[mkt * nJM * nI * nsP + i * nJM * nsP + (j - mStart) * nsP + l];
                    // For epsilon shocks order shouldn't matter much anyway. For consistency still try to keep it analogous to price beliefs.
                    expEpsShocksM[i * nJM * nsP + (j - mStart) * nsP + l] = expEpsShocks[mkt * nJM * nI * nsP + i * nJM * nsP + (j - mStart) * nsP + l];
                    
                } // end loop over price belief and epsilon simulations.
            }     // end loop over consumers.
        }         // end loop over current choices.

        // Destroy previous period's market share predictions only here since they were needed for lagged shares in this period.
        if (mkt > 0)
        {
            // Destroy object from previous month and create new one for current month.
            mxDestroyArray(pSharesiMx);
        }
        // Predicted shares for each combination of: current produt-simulated consumer
        pSharesiMx = mxCreateDoubleMatrix(nJM * nI, 1, mxREAL);
        pSharesiM = mxGetDoubles(pSharesiMx);
        //}
        // FINISHED SLOTTING INPUT OBJECTS FOR MARKET m.
        ////////////////////////////////////////////////////////////////////////

        // Start contraction mapping code for market mkt.
        while (norm > tol && iter < maxiter)
        {
            // Initialize container for average PCW usage statistic.
            // Make sure to set this to zero here because the code uses +=.
            PCWUsageAvg[mkt] = 0.0;
            maxabsdev = 0.0;
            // Loop over individual consumers.
            for (ii = 0; ii < nI; ii++)
            {
                // Loop over previous period products.
                for (k_lag = 0; k_lag < nJM; k_lag++)
                {
                    // Initialize different types of inclusive value for consumer ii and previous choice k_lag.
                    // These have to be initialized to zero because we use += to sum over products!
                    incValFBM[(ii * nJM) + k_lag] = 0.0; // IncVal for PCW non-users.
                    incValFAM[(ii * nJM) + k_lag] = 0.0; // IncVal for PCW users ,i.e., from full choice set.

                    // Loop over current products.
                    for (k = 0; k < nJM; k++)
                    {
                        // Compute numerator of conditional choice probabilities.
                        if (k == k_lag) // condition if SC are contract-specific.
                        {
                            // Compute exponentiated utility when NOT paying switching costs.
                            exputilM[(ii * nJM * nJM) + (nJM * k_lag) + k] = expmvaloldM[k] * expmuM[(ii * nJM) + k];
                        }
                        else
                        {
                            // Compute exponentiated utility when paying switching costs.
                            exputilM[(ii * nJM * nJM) + (nJM * k_lag) + k] = expmvaloldM[k] * expmuM[(ii * nJM) + k] / expPsiM[ii];
                        }
                       
                        for (l = 0; l < nsP; l++)
                        {
                            // Subsitute actual price data with price belief data.
                            exputilBelief = exputilM[(ii * nJM * nJM) + (nJM * k_lag) + k] * expPBeliefDiffM[ii * nJM * nsP + k * nsP + l];

                            // Fill array of simulated utilities.
                            // Conditional on consumer type and previous contract choice
                            // Order: current product - simulation draw
                            // When using PCW.

                            // Utility expected from contract k with simulation draw l if PCW was used.
                            USimPCW[k * nsP + l] = (awareiAdvM[(ii * nJM * nJM) + (nJM * k_lag) + k] *
                                                    exputilM[(ii * nJM * nJM) + (nJM * k_lag) + k] *
                                                    expEpsShocksM[ii * nJM * nsP + k * nsP + l]) +
                                                   ((1.0 - awareiAdvM[(ii * nJM * nJM) + (nJM * k_lag) + k]) *
                                                    exputilBelief *
                                                    expEpsShocksM[ii * nJM * nsP + k * nsP + l]);
                            // Expected utility from contract k and simulation draw l When not using PCW.
                            USimNoPCW[k * nsP + l] = (awareiAdvM[(ii * nJM * nJM) + (nJM * k_lag) + k] *
                                                      exputilM[(ii * nJM * nJM) + (nJM * k_lag) + k] * expEpsShocksM[ii * nJM * nsP + k * nsP + l]) +
                                                     ((1.0 - awareiAdvM[(ii * nJM * nJM) + (nJM * k_lag) + k]) *
                                                      0.0); // 0.0 is exponentiated utility from negative infinity, i.e., choice is not available
                                                            // IMPORTANT: If choice is not available, it does not get an epsilon shock either!

                        } // end loop over simulation draws

                        // For each consumer and lagged contract choice:
                        // Updated inclusive value from fallback choice set (without PCW usage).
                        incValFBM[(ii * nJM) + k_lag] += awareiAdvM[(ii * nJM * nJM) + (nJM * k_lag) + k] *
                                                         exputilM[(ii * nJM * nJM) + (nJM * k_lag) + k];
                        // Updated inclusive value from full choice set (with PCW usage).
                        incValFAM[(ii * nJM) + k_lag] += exputilM[(ii * nJM * nJM) + (nJM * k_lag) + k];
                    } // END OF LOOP OVER CURRENT PRODUCTS.
                   
                    // Compute expected maximum utility from using PCW versus not using PCW.
                    // This is done for each consumer type and previous contract choice.
                    EUMaxNoPCW[0] = 0.0;
                    EUMaxPCW[0] = 0.0;
                    // Average over simulation draws.
                    for (l = 0; l < nsP; l++)
                    {
                        // To start, set maximum simulated utility to utility from first contract.
                        UMaxSimNoPCW[0] = USimNoPCW[l];
                        UMaxSimPCW[0] = USimPCW[l];
                        // Find maximum (over products) of simulated utilities
                        // by checking all other contracts simulated utilities.
                        for (k = 1; k < nJM; k++)
                        {
                            if (UMaxSimNoPCW[0] < USimNoPCW[k * nsP + l])
                            {
                                UMaxSimNoPCW[0] = USimNoPCW[k * nsP + l];
                            }
                            if (UMaxSimPCW[0] < USimPCW[k * nsP + l])
                            {
                                UMaxSimPCW[0] = USimPCW[k * nsP + l];
                            }
                        } // end looping over products, now max for each simulation draw is found.

                        // Update average of UMax for both PCW usage decisions.
                        // The division could be moved outside of loop over simulation draws, but that is more likely to run into overflow if utility scaling is unfortunate. Unlikely to happen here.
                        EUMaxNoPCW[0] += (UMaxSimNoPCW[0] * invnsP);
                        EUMaxPCW[0] += (UMaxSimPCW[0] * invnsP);
                    } // end of loop over simulation draws; now expected MaxU are found.
                    // PCW usage prediction with logit smoothing (should behave better in optimization since it's a smooth function now.) See, e.g., Train textbook, page 121.
                    PCWSmoothAux = pow((EUMaxPCW[0] / expKappaM[ii]) / EUMaxNoPCW[0], invLS);
                    PCWUsageM[(ii * nJM) + k_lag] = PCWSmoothAux / (1.0 + PCWSmoothAux);
                } // END OF LOOP OVER PREVIOUS PRODUCTS.

                // Compute predicted individual-specific market shares.
                for (k = 0; k < nJM; k++) // loop over current period's products.
                {
                    // Initialize predicted market share for each consumer type and product.
                    // Main market share prediction.
                    pSharesiM[(ii * nJM) + k] = 0.0;
                    // Share prediction conditional on using PCW and not using the PCW.
                    pSharesiPCWM[(ii * nJM) + k] = 0.0;
                    pSharesiNoPCWM[(ii * nJM) + k] = 0.0;
                    for (k_lag = 0; k_lag < nJM; k_lag++) // Loop over previous choices.
                    {
                        // Full share prediction assuming smoothed PCW usage prediction.
                        // Compute CCP for consumer ii choosing contract k when previous subscription was k_lag.
                        ccpAux = ((PCWUsageM[(ii * nJM) + k_lag] *                  // CCP contribution of PCW users
                                   exputilM[(ii * nJM * nJM) + (nJM * k_lag) + k] / // CCP of PCW users
                                   incValFAM[(ii * nJM) + k_lag]) +
                                  //  sharesiLagM[(ii * nJM) + k_lag] + // finally, multiply with previous period's market share of product k_lag

                                  ((1.0 - PCWUsageM[(ii * nJM) + k_lag]) *            // CCP contribution of PCW non-users
                                   awareiAdvM[(ii * nJM * nJM) + (nJM * k_lag) + k] * // whether ii is aware of k without PCW
                                   exputilM[(ii * nJM * nJM) + (nJM * k_lag) + k] /   // CCP of PCW non-users
                                   incValFBM[(ii * nJM) + k_lag]));

                        pSharesiM[(ii * nJM) + k] += // existing share accumulated from looping over previous products before k_lag.
                            ccpAux *
                            sharesiLagM[(ii * nJM) + k_lag]; // finally, multiply with previous period's market share of product k_lag

                        // Write CCP to market-specific vector.
                        // Order of assignments: (1) consumers (2) lagged choice (3) current choice
                        // Structure is analogous to exputilM.
                        ccpiM[ii * nJM * nJM + nJM * k_lag + k] = ccpAux;
                        // Share prediction conditional on using PCW.
                        pSharesiPCWM[(ii * nJM) + k] +=                      //existing share from looping over previous period's products before k_lag.
                            exputilM[(ii * nJM * nJM) + (nJM * k_lag) + k] / // CCP of PCW users
                            incValFAM[(ii * nJM) + k_lag] *
                            sharesiLagM[(ii * nJM) + k_lag]; // finally, multiply by previous period's market share of product k_lag
                        // Share prediction conditional on not using PCW.
                        pSharesiNoPCWM[(ii * nJM) + k] +=                      //existing share from looping over previous period's products before k_lag.
                            awareiAdvM[(ii * nJM * nJM) + (nJM * k_lag) + k] * // whether ii is aware of k without PCW
                            exputilM[(ii * nJM * nJM) + (nJM * k_lag) + k] /   // CCP of PCW non-users
                            incValFBM[(ii * nJM) + k_lag] *
                            sharesiLagM[(ii * nJM) + k_lag]; // finally, multiply by previous period's market share of product k_lag

                    } // end loop over previous products.

                    // PCW Usage depends on previous contract and consumer type.
                    // Need to sum over all previous products (k) and average over consumers ii.
                    // Careful: role of k and k_lag changes here!
                    PCWUsageAvg[mkt] += PCWUsageM[(ii * nJM) + k] * sharesiLagM[(ii * nJM) + k];
                } // end loop over current products.
            }     // end loop over consumers.
            // Average aggregated PCW usage over consumers.
            // Overflow is not an issue here, since PCW usage is probability between zero and one.
            PCWUsageAvg[mkt] = PCWUsageAvg[mkt] * invnI;

            /*=========================================================
            * Aggregate over consumers to obtain predicted shares for month mkt,
            * calculates expmval and deviations, and updates expmvalold
            *=========================================================*/
            // Only debug mode: Aggregate market share predictions should sum to one.
            // sum_agg_shares = 0.0;
            for (k = 0; k < nJM; k++) // loop over curent products.
            {
                pSharesM[k] = 0.0;
                for (ii = 0; ii < nI; ii++) // loop over consumers.
                {
                    // Main prediction for shares.
                    // In loop, only sum over individual market shares. Divide after loop.
                    pSharesM[k] += pSharesiM[(ii * nJM) + k]; // * invnI;
                    // Prediction for shares conditional on using PCW.
                    // pSharesPCWM[k] += pSharesiPCWM[(ii * nJM) + k] * invnI;

                } // end loop over consumers
                // Because of no risk of overflow here, decrease number of computational operations.
                pSharesM[k] *= invnI;
               
                if (k < nJM - 1 && contrMap_dummy == 1)
                {
                    // Classic version without any dampening.
                    // Update actual mean utilities for inside goods.
                    expmvalM[k] = expmvaloldM[k] * sharesM[k] / pSharesM[k]; // Update based on market share difference (data vs predicted)
                }
                else if (contrMap_dummy == 0 || k == nJM - 1)
                // Make sure not to update mean utility of outside good.
                {
                    expmvalM[k] = expmvaloldM[k];
                }

                // Test for convergence.
                absdev = fabs(expmvalM[k] - expmvaloldM[k]);
                if (absdev > maxabsdev)
                {
                    maxabsdev = absdev;
                }
                // Update mean utilities for next iteration of CM.
                expmvaloldM[k] = expmvalM[k];
                // Only diagnostics - debug mode: Sum shares to check that they sum to one.
                // sum_agg_shares += pSharesM[k];
            } // end loop over current products
            iter = iter + 1;
            norm = maxabsdev;
            // If we only want to predict market shares without updating mean utilities, set norm to zero
            // after first iteration.
            if (contrMap_dummy == 0)
            {
                //              mexPrintf("No CM performed, only prediction of market shares.\n");
                norm = 0;
            } // end of computations when contrMap_dummy==0.
            // mexPrintf("Indicator for contraction mapping: %d \n",contrMap_dummy);

        } // END OF CONTRACTION MAPPING WHILE LOOP.

        // Churn rate computation outside of contraction mapping is computationally less demanding.
        // Loop over consumers.
        for (ii = 0; ii < nI; ii++)
        { // loop over simulated consumers.
            for (k = 0; k < nJM; k++)
            { // loop over (current products).
                // Initialize predicted churn rates for each consumer type and product.
                // Below we subtract the CCP of sticking with firm (not contract),
                // so that churn rate = 1 - CCP(staying with firm).
                pChurniM[(ii * nJM) + k] = 1.0;
                // 3 case distinctions across products:
                // (a) first/conventional contract of two-contract firms: 0-2-5-7
                // (b) second/green contract of two-contract firms: 1-3-6-8
                // (c) single-contract firms (green only): 4-9-10
                // Record churn rate predictions for each individual consumer and product. Note that churn rate is firm, not contract specific.
                // if (k == k_lag)
                if (k == 0 || k == 2 || k == 5 || k == 7) // first contract of two-contract firms.
                {
                    // Version based on CCP matrix estimate.
                    pChurniM[(ii * nJM) + k] -=
                        ((ccpiM[(ii * nJM * nJM) + nJM * k + k] + ccpiM[(ii * nJM * nJM) + nJM * k + k + 1])               // CCP of going from first contract to any contract of firm
                             * sharesiLagM[(ii * nJM) + k]                                                                 // Pr of being with first contract of firm in previous period
                         + (ccpiM[(ii * nJM * nJM) + nJM * (k + 1) + k] + ccpiM[(ii * nJM * nJM) + nJM * (k + 1) + k + 1]) // CCP of going from first contract to any contract of firm
                               * sharesiLagM[(ii * nJM) + k + 1])                                                          // Pr og being with second contract of firm in previous period
                        / (sharesiLagM[(ii * nJM) + k] + sharesiLagM[(ii * nJM) + k + 1]);                                 // normalize by market sahre of firm in previous period
                }
                else if (k == 1 || k == 3 || k == 6 || k == 8) // second contract of two-contract firms
                {                                              // This is just the same churn rate as the one computed in the first case, since it is only at the firm level.
                    pChurniM[(ii * nJM) + k] = pChurniM[(ii * nJM) + k - 1];
                }
                else // single-contract firms.
                {
                    // Version based on CCP matrix estimate.
                    pChurniM[(ii * nJM) + k] -=
                        ccpiM[(ii * nJM * nJM) + nJM * k + k]; //simple case because these firms only have one product, so just subtract CCP.
                }                                              // end case distinction for churn rate computation.
            }                                                  // end loop over products for churn rate computation.
        }                                                      // end loop over consumer types for churn rate computation.

        // After convergence of market-specific mean utilities:
        // Slot output objects to global output arrays.
        // 1. Month-contract specific mean utilities.
        for (j = 0; j < nJM; j++) // For deltas, only loop over current contracts.
        {
            // Month-contract specific mean utilities to expmval array.
            expmval[j + mStart] = expmvalM[j];
        } // end slotting loop over products (for mean utilities)

        // 2. Consumer-specific output objects.
        for (ii = 0; ii < nI; ii++)   // First, loop over consumers.
        {                             // Note: This slotting assumes constant set of products over months.
            for (j = 0; j < nJM; j++) // Second, loop over current contracts.
            {
                // Slot consumer-contract specific predicted market shares.
                pSharesi[mkt * nJM * nI + ii * nJM + j] = pSharesiM[ii * nJM + j];
                // Analogous for shares conditional on PCW usage.
                pSharesiPCW[mkt * nJM * nI + ii * nJM + j] = pSharesiPCWM[ii * nJM + j];
                pSharesiNoPCW[mkt * nJM * nI + ii * nJM + j] = pSharesiNoPCWM[ii * nJM + j];
                // Slot consumer-contract specific churn rates.
                pChurni[mkt * nJM * nI + ii * nJM + j] = pChurniM[ii * nJM + j];
                // Slot consumer-previous contract specific PCW usage.
                PCWUsagei[mkt * nJM * nI + ii * nJM + j] = PCWUsageM[ii * nJM + j];
                // Slot consumer-previous contract specific surplus to global vector.
                // If consumer uses PCW.
                // Recall that incVal objects are exponentiated, as well as kappa object.
                cSurplusPCWi[mkt * nJM * nI + ii * nJM + j] = log(incValFAM[(ii * nJM) + j] / expKappaM[ii]);
                // If consumer does not use PCW.
                cSurplusNoPCWi[mkt * nJM * nI + ii * nJM + j] = log(incValFBM[(ii * nJM) + j]);
            } // end slotting loop over products.
        }     // end slotting loop over consumers.
        // Store month-specific convergence measures for diagnostics.
        mktnorm[mkt] = norm;
        mktiter[mkt] = iter;

        if (contrMap_dummy == 0)
        {
            // If no contraction mapping is performed, also slot CCP matrix elements.
            // Pay attention to slotting of elements here:
            // First, consumer types, then current contract, then lagged contract.
            // So, reshape in order to have current choice in columns and lagged choice in rows in model_myopic.m
            for (i = 0; i < nI; i++) // loop over consumer types.
            {
                for (k = 0; k < nJM; k++) // loop over current contracts.
                {
                    for (k_lag = 0; k_lag < nJM; k_lag++)
                    {
                        // ccpMat[mkt * nJM * nJM + k * nJM + k_lag] = ccpMatM[k * nJM + k_lag];
                        ccpMat[mkt * nJM * nJM * nI + i * nJM * nJM + k * nJM + k_lag] = ccpiM[i * nJM * nJM + k * nJM + k_lag];
                    } // end loop over previous choices.
                }     // end loop over current products.
            }         // end loop over consumer types.
        }             // end code when contrMap_dummy==0.

        /* Destroying market-specific data to prevent memory leaks */
        // Individual-specific market shares are destroyed only after
        // they have been used to create lagged market shares above.
        mxDestroyArray(expmvaloldMx);
        mxDestroyArray(expmuMx);
        // mxDestroyArray(incValMx);
        mxDestroyArray(incValFBMx);
        // mxDestroyArray(incValBeliefMx);
        mxDestroyArray(incValFAMx);
        mxDestroyArray(sharesMx);
        mxDestroyArray(awareiAdvMx);
        mxDestroyArray(expPBeliefDiffMx);
        mxDestroyArray(expEpsShocksMx);
        mxDestroyArray(exputilMx);
        mxDestroyArray(EUSearchNetMx);
        mxDestroyArray(EUNoSearchMx);
        mxDestroyArray(PCWUsageMx);
        mxDestroyArray(ccpiMx);
        // The individual predicted shares are needed or next month's prediction,
        // therefore, don't destroy it here!
        // mxDestroyArray(pSharesiMx);
        mxDestroyArray(pSharesiPCWMx);
        mxDestroyArray(pSharesiNoPCWMx);
        mxDestroyArray(pChurniMx);
        mxDestroyArray(pSharesMx);
        mxDestroyArray(expmvalMx);
        mxDestroyArray(USimNoPCWx);
        mxDestroyArray(USimPCWx);
        mxDestroyArray(UMaxSimNoPCWx);
        mxDestroyArray(UMaxSimPCWx);
        mxDestroyArray(EUMaxNoPCWx);
        mxDestroyArray(EUMaxPCWx);
    } // end of loop over markets
    if (contrMap_dummy == 0)
    {
        mexPrintf("No CM performed, only prediction of market shares.\n");
        mexPrintf("Computed matrix of CCPs for each period..\n");
    }
} // end of contrMap_ES

/*====================================================================
 * mexFunction is the gateway beween Matlab and the contrMap function.
 *===================================================================*/
void mexFunction(int nlhs, mxArray *plhs[],       /* Setup output arguments */
                 int nrhs, const mxArray *prhs[]) /* Setup input arguments */
{
    // Declare input arguments.
    double tol;             /* contraction mapping tolerance */
    double maxiter;         /* max iterations in contraction mapping */
    double *expmvalold;     /* initial vector of mean valuations */
    double *expmu;          /* matrix of individual-specific utility deviations */
    double *cdindex;        /* vector identifies last obs of each market */
    double *shares;         /* vector of observed market shares */
    double *shares_init;    /* vector of observed market shares in initial period for each consumer type */
    double *exp_psi;        /* vector of switching cost parameters */
    double *exp_kappa;      /* vector of PCW search cost parameters */
    double *aware_adv;      /* indicator vector of awareness through advertising */
    double *expPBeliefDiff; /* matrix of price belief data normalized for actual price and difference between inside and outside good */
    double *expEpsShocks;   /* vector of logit shocks for PCW usage simulation. */
    _Bool contrMap_dummy;   /* indicator whether contraction mapping is run or only market share predictions computed */
    // Extract several important dimensions from input arguments.
    size_t nJ;  /* number of product-mkt observations: J * T */
    size_t nI;  /* number of individual consumer draws  */
    size_t nM;  /* number of markets/months T  */
    size_t nsP; /* number of simulated price beliefs NS_p_eps */

    // Declare output arguments.
    double *expmval;        /* output vector of updated mean valuations */
    double *PCWUsagei;      /* output vector of predicted PCW usage */
    double *pSharesi;       /* output vector of predicted individual market shares */
    double *pSharesiPCW;    /* output vector of predicted individual market shares conditional on using PCW */
    double *pSharesiNoPCW;  /* output vector of predicted individual market shares conditional on not using PCW */
    double *pChurni;        /* output vector of predicted churn rates */
    double *cSurplusPCWi;   /* output vector of predicted consumer surplus for PCW users*/
    double *cSurplusNoPCWi; /* output vector of predicted consumer surplus for PCW non-users*/
    double *ccpMat;         /* output vector for CCP matrix (for all periods and each simulated consumers) */
    double *mktnorm;        /* output vector of realized tolerances for each market */
    double *mktiter;        /* output vector of # iterations for each market */

    // prhs[xxx]: xxx denotes the order of the input arguments.
    /* get the value of tol  */
    tol = mxGetScalar(prhs[0]);

    /* get the value of maxiter  */
    maxiter = mxGetScalar(prhs[1]);

    /* create a pointer to the real exmvalold data  */
    expmvalold = mxGetDoubles(prhs[2]);

    /* get dimensions of the inputs */
    nJ = mxGetM(prhs[3]);            // total number of product-market observations, here: T x J = 53 x 11 = 583
    nI = mxGetN(prhs[3]);            // number of simulated consumers, here: probably 1000 or 2000
    nM = mxGetM(prhs[4]);            // number of markets, for us: months T=53
    nsP = mxGetM(prhs[8]) / nJ / nI; // number of simulation draws for PCW decision, probably around 30 for us.
    //     mexPrintf("nJ=%d; nI=%d; nM=%d\n", nJ, nI, nM);

    /* create a pointer to the real expmu data  */
    expmu = mxGetDoubles(prhs[3]);
    /* create a pointer to the real cdindex data */
    cdindex = mxGetDoubles(prhs[4]);
    /* create a pointer to the real share data  */
    shares = mxGetDoubles(prhs[5]);
    /* create a pointer to the real initial share data  */
    shares_init = mxGetDoubles(prhs[6]);
    /* create a pointer to the awareness-through-advertising data  */
    aware_adv = mxGetDoubles(prhs[7]);
    // mexPrintf("aware_adv[0] = %f \n", aware_adv[0]);

    /* create a pointer to the price belief data  */
    // Product-month-inviduals in rows, simulation draws in columns.
    expPBeliefDiff = mxGetDoubles(prhs[8]);
    // mexPrintf("expPBeliefDiff[0] = %f \n", expPBeliefDiff[0]);
    // IID-Logit Shocks for PCW usage simulation.
    expEpsShocks = mxGetDoubles(prhs[9]);
    // mexPrintf("expEpsShocks[2] = %f \n", expEpsShocks[2]);
    /* create a pointer to the switching cost parameter vector */
    exp_psi = mxGetDoubles(prhs[10]);
    /* create a pointer to the PCW search cost parameter vector */
    exp_kappa = mxGetDoubles(prhs[11]);
    /* Dummy indicator for whether CM or market share prediction is run*/
    contrMap_dummy = mxGetScalar(prhs[12]);
    // mexPrintf("Value of contrMap dummy is %d.\n", contrMap_dummy);

    /* create the output data: expmval */
    plhs[0] = mxCreateDoubleMatrix((mwSize)nJ, 1, mxREAL);
    /* create the output data: PCWUsagei */
    // PCW usage decisions for each month, previous choice and simulated consumer type.
    plhs[1] = mxCreateDoubleMatrix((mwSize)nJ * nI, 1, mxREAL);
    /* create the output data: pSharesi */
    // Output predicted share for each product, month and consumer type!
    plhs[2] = mxCreateDoubleMatrix((mwSize)nJ * nI, 1, mxREAL);
    // Analogous for individual shares conditional on PCW and non PCW usage.
    plhs[3] = mxCreateDoubleMatrix((mwSize)nJ * nI, 1, mxREAL);
    plhs[4] = mxCreateDoubleMatrix((mwSize)nJ * nI, 1, mxREAL);

    /* create the output data: pChurni */
    // Churn rate depends on month, product and consumer type
    plhs[5] = mxCreateDoubleMatrix((mwSize)nJ * nI, 1, mxREAL);
    /* create the output data: cSurplusi for both PCW users and non-users */
    // Consumer surplus depends on month, previous choice and consumer type.
    plhs[6] = mxCreateDoubleMatrix((mwSize)nJ * nI, 1, mxREAL);
    plhs[7] = mxCreateDoubleMatrix((mwSize)nJ * nI, 1, mxREAL);
    // Create CCP matrix on individual level: J x J x T = 11 x 11 x 53.
    plhs[8] = mxCreateDoubleMatrix((mwSize)nJ * nJ / nM * nI, 1, mxREAL);
    /* create the output data: mktnorm */
    plhs[9] = mxCreateDoubleMatrix((mwSize)nM, 1, mxREAL);
    /* create the output data: mktiter */
    plhs[10] = mxCreateDoubleMatrix((mwSize)nM, 1, mxREAL);

    // plhs[xxx]: xxx denotes the order of the output! arguments.
    /* get a pointer to the real data in the output matrix */
    expmval = mxGetDoubles(plhs[0]);
    /* get a pointer to the real data in the output matrix */
    PCWUsagei = mxGetDoubles(plhs[1]);
    /* get a pointer to the real data in the output matrix */
    pSharesi = mxGetDoubles(plhs[2]);
    pSharesiPCW = mxGetDoubles(plhs[3]);
    pSharesiNoPCW = mxGetDoubles(plhs[4]);
    /* get a pointer to the real data in the output matrix */
    pChurni = mxGetDoubles(plhs[5]);
    /* get a pointer to the real data in the output matrix */
    cSurplusPCWi = mxGetDoubles(plhs[6]);
    cSurplusNoPCWi = mxGetDoubles(plhs[7]);
    /* get a pointer to CCP matrix */
    ccpMat = mxGetDoubles(plhs[8]);
    /* get a pointer to the realized tolerance */
    mktnorm = mxGetDoubles(plhs[9]);
    /* get a pointer to the realized tolerance */
    mktiter = mxGetDoubles(plhs[10]);

    /* call the computational routine */
    contrMap(tol, maxiter, expmvalold, expmu,
             (mwSize)nJ, (mwSize)nI, (mwSize)nM, (mwSize)nsP,
             cdindex, shares, shares_init, aware_adv, expPBeliefDiff, expEpsShocks, exp_psi, exp_kappa, contrMap_dummy,
             mktnorm, mktiter, expmval, PCWUsagei, pSharesi, pSharesiPCW, pSharesiNoPCW, pChurni, cSurplusPCWi, cSurplusNoPCWi, ccpMat);
}
