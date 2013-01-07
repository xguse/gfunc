"""
####################
stats.py
####################
Code supporting calculations of statistical probabilities for gfunc.
"""

from decimal import Decimal
import numpy as np
import operator as o
# see bottom for conditional import of "bestChoose" 



def basic_bootstrap_est(vec,reps=1000):
    """
    *GIVEN:*
        * ``vec`` = vector of sample values
        * ``reps`` = number of resampling reps
    
    *DOES:*
        * Resample w/ replacement ``reps`` times and record the medians
        * Calculate stdv of resampled medians which should approach
          the actual SE as ``reps`` approaches ``inf``.
        * Calculate the 95% CI bounds.
       
       
    
    *RETURNS:*
        * tuple([*median of resampled medians*, *SE est*, *loBound*, *hiBound*])
    """
    
    reSampledMeds = []
    sampleSize = len(vec)
    
    for rep in range(reps):
        sample = []
        for draw in range(sampleSize):
            randIndex = np.random.randint(0,sampleSize)
            sample.append(vec[randIndex])
        reSampledMeds.append(np.median(sample))
        
    seEstimate = np.std(reSampledMeds,ddof=1)
    medOfMeds = np.median(reSampledMeds)
    
    return (medOfMeds,
            seEstimate,
            np.percentile(reSampledMeds,2.5),
            np.percentile(reSampledMeds,97.5))

def benjHochFDR(table,pValColumn=-1):
    """
    *GIVEN:*
        * ``table``: 2D list(*hypothesis*,*p-value*) hypothesis could = *geneName* tested for enrichment
        * ``pValColumn``: integer of column index containing the *p-value*.
    *DOES:*
        * Calculates the Benjamini-Hochberg adjusted *p-values*
    *RETURNS:*
        * a new version of ``table`` with an extra column added to the end representing the BH corrected *p-values*
    
    """
    assert type(pValColumn) == type(1),\
           "ERROR: pValColumn must be int type!"
    # Sort table from highest to lowest after converting them all to floats.
    for i in range(len(table)):
        table[i][pValColumn] = float(table[i][pValColumn])
    table.sort(key=lambda x: x[pValColumn])
    table.reverse()
    
    n = len(table)
    
    lastPval = table[0][pValColumn]
    for i in range(len(table)):
        p    = table[i][pValColumn]
        adj  = (float(n)/(n-i))
        adjP = p*adj
        miN  = min(adjP,lastPval)
        table[i].append(miN)
        lastPval = table[i][-1]
    
    table.reverse()
    return table




def binComb(n, k):
    """
    *GIVEN:*
        * ``n``
        * ``k``
    *DOES:*
        * Computes ``n`` *choose* ``k``. 
    *RETURNS:*
        * The number of ways ``k`` objects can be sampled from a population of size ``n``.
    
    """
    if (k > n): return 0
    if (k < 0): return 0
    if (k > int(n/2)):
        k = n - k
    rv = 1
    for j in range(0, k):
        rv *= n - j
        rv /= j + 1
    return rv


def choose(n, k):
    if (k > n): return 0
    if (k < 0): return 0
    ntok = 1
    for t in xrange(min(k, n-k)):
        ntok = ntok*(n-t)//(t+1)
    return ntok

# determine whether we can use gmpy and set "bestChoose" to that or the standby.
try:
    from gmpy import comb as bestChoose #timeit says gmpy.comb is 10 times faster than choose!
    print "bestChoose is 'comb' from 'gmpy'."
except ImportError:
    bestChoose = choose # timeit says that choose is 7% faster than binComb
    print "bestChoose is 'choose' from 'rSeq'."
    
    

def hypergeoP(n,i,m,N):
    """
    | Calculates the non-cumulative hypergeometric *p-value* for variables:
    | 
    | ``n`` = number of positives in population
    | ``i`` = number of positives in sample
    | ``m`` = number of negatives in population
    | ``N`` = sample size
    | 
    | *P(x=i) = (choose(n,i)choose(m,N-i))/choose(n+m,N)*
    | 
    | For more details -> http://mathworld.wolfram.com/HypergeometricDistribution.html
    
    """
    return (bestChoose(n,i)*bestChoose(m,N-i))/float(bestChoose(n+m,N))




def cumHypergeoP(n,i,m,N):
    """
    | Calculates the cumulative hypergeometric *p-value* for variables:
    | 
    | ``n`` = number of positives in population
    | ``i`` = number of positives in sample
    | ``m`` = number of negatives in population
    | ``N`` = sample size
    | 
    | *P(i) = sum([as i->N] (choose(n,i)choose(m,N-i))/choose(n+m,N))*
    | 
    | For more details -> http://mathworld.wolfram.com/HypergeometricDistribution.html
    
    """

    cumPVal = 0

    for x in range(i,N+1):
        cumPVal = cumPVal + hypergeoP(n,x,m,N)

    return cumPVal


def binomialPval(n,k,p):
    """
    *RETURNS:* 
        * exact binomial P-value.
    
    | ``n`` = number of trials
    | ``k`` = number of successes
    | ``p`` = probability of a success
    | 
    | *P(k succeses in n trials) = choose(n,k) (p^k) ((1-p)^(n-k))*
    
    """
    
    return bestChoose(n,k) * p**k * (1-p)**(n-k)

def binomialPval_gte(n,k,p):
    """
    *RETURNS:* 
        * binomial *p-value* of ``k`` or greater successes in ``n`` trials with probability of success for each trial ``p``.
    
    | ``n`` = number of trials
    | ``k`` = number of successes
    | ``p`` = probability of a success
    | 
    | *sum( choose(n,k) (p^k) ( (1-p)^(n-k) ) )* as ``k`` goes from ``k`` to ``n``
    
    """    

    cumPVal = 0
    for k in range(k,n+1):
        cumPVal = cumPVal + binomialPval(n,k,p)
    return cumPVal