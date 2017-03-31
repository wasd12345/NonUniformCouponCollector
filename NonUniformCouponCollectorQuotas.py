# -*- coding: utf-8 -*-
"""
Created on Thu Mar 16 20:32:41 2017

@author: gk

Implements results from the following paper:
Coupon Collecting with Quotas; Russell May; 2008.
http://www.combinatorics.org/ojs/index.php/eljc/article/view/v15i1n31/pdf
    
"""


import numpy as np
import matplotlib.pyplot as plt
from itertools import chain, combinations
import time

from sympy import factorial
from sympy.functions.combinatorial.factorials import binomial



def NonUniformCouponCollectorQuotas__AnalyticalMean(quota_vector,probability_vector):
    """
    For the Coupon Collector Problem, with:
    Non-uniform coupon probabilities, and
    non-uniform coupon quotas.
    
    Calculates the expected number of trials until the collection is complete,
    meaning that the quota on each coupon is met.
    
    Implements the following paper:
    Coupon Collecting with Quotas; Russell May; 2008.
    http://www.combinatorics.org/ojs/index.php/eljc/article/view/v15i1n31/pdf
    [equation 6]
    
    
    Inputs:
    quota_vector - vector of length Ncoupons. Element i is the quota for coupon i
    probability_vector - vector of length Ncoupons. Element i is the probability of 
    getting coupon i. The probabilities are assumed fixed for every draw.
    
    Output:
    ExpectedValue - float. The expected value of the number of drwas until the 
    collection is complete.
    
    
    *** Is extremely slow and/or memory intensive for more than a few coupons (8-10),
    or when one or more of the coupons has a large quota
    """
    
    #Vectors must both be length Ncoupons:
    if quota_vector.size != probability_vector.size:
        raise Exception('Both vectors must both be length Ncoupons')
    else:
        Ncoupons = quota_vector.size
        indices = np.arange(Ncoupons)
        
    #Make sure quota_vector is integers:
#    if type(quota_vector)!=int:
#        print 'quota_vector converted to ints'
#        quota_vector = quota_vector.astype(int)
    
    #Generate M, the set of all combinations of non-empty subsets of letters from 
    #the dictionary L (the set of all letters or coupons).
    #combinations(indices,1), combinations(indices,2), ..., combinations(indices,Ncoupons)
    #e.g. when Ncoupons=3: is: [(0,), (1,), (2,), (0, 1), (0, 2), (1, 2), (0, 1, 2)]
    #is the set of all combinations of indices, for sets of length 1 to Ncoupons
    
    #If generate in advance, even for small Ncoupons (e.g. >25) starts to take up too much memory, 
    #even though there are only < 10M tuples in list (although maybe ~20 elements per tuple average so is big in memory)
    #The cumulative number of elements over all tuples in this list is the nth Mersenne number.
    #M = chain(*[combinations(indices,i) for i in np.arange(1,Ncoupons+1)])])
    #M_cardinalities = [len(kk) for kk in M]
    
    #So instead, just step through iterable
    #But this will have over a billion terms after only ~ Ncoupons=30. So not very helpful either.
    M = chain(*[combinations(indices,i) for i in np.arange(1,Ncoupons+1)])
    
    
#    quota_vector = list(quota_vector.astype(int)) #To make iterable
    
    
    ExpectedValue = 0.
    
    #Iterate through each set M
    for mm in M:
#        print 'mm', mm
        
        factor = (-1.)**(len(mm) + 1)
#        print 'factor', factor
        
        #Terms that are constant for a given mm:
        q = quota_vector[np.array(mm)]
        p = probability_vector[np.array(mm)]
        p_sum = np.sum(p)
#        print 'q', q
#        print 'p', p
#        print 'p_sum', p_sum
        
        #Generate all possible r vector values:
        # 0 <= r_l < q_l
        qvals = [np.arange(i) for i in q]
#        print qvals
        
    
        r_grid = np.meshgrid(*qvals)
        #To keep proper ordering, need to swap axes 0 and 1 because of way meshgrid works
        if len(r_grid[0].shape)>=2:
            r_grid = [np.swapaxes(i,0,1) for i in r_grid]
        FLAT = np.array([i.flatten() for i in r_grid])
        FLAT_sum = np.sum(FLAT,axis=0)
        Nr = FLAT.shape[1]
#        print r_grid
#        print len(r_grid)
#        print r_grid[0].shape
#        print np.sum(r_grid,axis=0)
#        print FLAT
        
        




        #product of p^r values
        p = np.expand_dims(p,axis=0).T
        p_r = np.prod(np.repeat(p,Nr,axis=1)**FLAT,axis=0)


        #Once Sympy developers fix factorial ufuncify, can easily vectorize this step:
        multichoose = np.array([factorial(i) for i in FLAT_sum]) / np.prod(np.array([[factorial(y) for y in x] for x in FLAT]),axis=0)
        
        numerator = multichoose * p_r    
        denominator = p_sum **(1+FLAT_sum)
##            print 'multichoose', multichoose
##            print 'factorial(r_sum)', factorial(r_sum)
##            print 'np.prod([factorial(i) for i in r])', np.prod([factorial(i) for i in r])
##            print 'numerator', numerator
##            print 'denominator', denominator
##            print 'numerator/denominator', numerator/denominator
##            print
#
#            
#            SUM += numerator/denominator
    
#        print 'SUM', SUM

        


        ExpectedValue += factor*np.sum(numerator/denominator)
    
    return ExpectedValue



    
    
    
    
    
def UniformCouponCollectorQuotas__AnalyticalMean(d,n):
    """
    For UNIFORM quotas and probabilities, since non-uniform is very 
    computationally intensive.
    
    Implements the following paper:
    Coupon Collecting with Quotas; Russell May; 2008.
    http://www.combinatorics.org/ojs/index.php/eljc/article/view/v15i1n31/pdf
    [equation 8]

    
    Inputs:
    
    d - integer. The number of coupons
    
    n - integer. The quota (assumed same for all coupons)
    
    
    Output:
    
    ExpectedValue - float. The EV for the number of draws needed to complete 
    the collection.
    
    
    **For trivial case of d=1, n=n, will give error. But obviously the exact 
    value of the mean in this case is just n.
    """
    
    #Because of notation in paper:
    n = n - 1    
    
    #n = Symbol('n', integer=True)
    m1_term = d*(n+1.) / (factorial(d)**n)
#    print m1_term
    
    #For r = 2, 3, ..., d; 
    #get the m_d summation, and append to list of values
    M = [m1_term]
    for r in xrange(2,d+1):
#        print r
        
        #term1 = [n + m_(r-1)] choose [n]
        #term2 = ((r-1)/r)^m_r
    
        #Vector of all possible m_r values to be summed over:
        m_r = np.arange(0, ((r-1)*n+1) + 1).astype(int)
        
        #term1 is the binomial coefficient term.
        #If sympy ufuncify is working, might be faster to vectorize it...
        #Also, this actually only needs to be generated once for the biggest 
        #value of r. For smaller values of r, it just uses a subset of the 
        #first few terms. So could save computation if this step becomes problem.
        term1 =  np.array([binomial(n+i, n) for i in m_r])
        
        #Term2 is the ratio term        
        term2 = ((r-1.)/r)**m_r
        
        #Cumulative sum of the elementwise product
        M += [term1*term2]
        
#        print term1
#        print term2
#        print
        
#    print M
#    print


#==============================================================================
# 
#==============================================================================
    FULL = M[-1]
    for r in xrange(1,d):
        #The cumulative sum
        FULL = np.cumsum(FULL)[n:]
        FULL *= M[-r-1]
    ExpectedValue = FULL[-1]    
    
#    print ExpectedValue
    return ExpectedValue







def NonUniformCouponCollectorQuotas__Simulation(quota_vector,probability_vector,Nsimulations,batch_method):
    """
    For the Coupon Collector Problem, with:
    Non-uniform coupon probabilities, and
    non-uniform coupon quotas.
    

    Empircally simulates:
    
    Calculates the expected number of trials until the collection is complete,
    meaning that the quota on each coupon is met.
    
    And calculates the variance in the number of trials until completion.
    
    
    Inputs:
    
    quota_vector - vector of length Ncoupons. Element i is the quota for coupon i

    probability_vector - vector of length Ncoupons. Element i is the probability of 
    getting coupon i. The probabilities are assumed fixed for every draw.    
    
    Nsimulations - int. Number of simulations to  run 
    [when not too many coupons or large quotas can do few hundred simulations pretty quickly]
    
    batch_method - True / False. Whether or not to process in batches.
    Default to True sice seems to be faster in batches.
    """
    
    #Vectors must both be length Ncoupons:
    if quota_vector.size != probability_vector.size:
        raise Exception('Both vectors must both be length N')
    else:
        Ncoupons = quota_vector.size
        indices = np.arange(Ncoupons)  
        
        
        
    #Container variable:
    #For each simulation, the number of trials needed to compelte the coupon collection
    Vals = np.zeros(Nsimulations).astype(int)
    
    for i in xrange(Nsimulations):


        if batch_method == False:
            #Iterative one by one way:
            counts_vector = 0.*quota_vector
            #Until the collection is complete, take random draws according to probability vector        
            while ~np.alltrue(counts_vector >= quota_vector):
                counts_vector[np.random.choice(indices,p=probability_vector)] += 1
            
            #The total number of trials required to complete the coupon collection
            #is the sum of counts of each coupon type.
            Vals[i] = counts_vector.sum()
        
        
        
        
        
        if batch_method == True:
            counts_vector = 0.*quota_vector
            batch_size = 200
            Nbatches = 0
            rep = np.repeat(np.expand_dims(quota_vector,axis=0),batch_size,axis=0)
            #Until the collection is complete, take random draws according to probability vector        
            while ~np.alltrue(counts_vector >= quota_vector):
    
                coupons = np.random.choice(indices,size=batch_size,p=probability_vector)
    
                #Make a (batch_size x Ncoupons) 0-1 matrix. Element (i,j) is 1 if simulation i drew coupon j.
                #So each row has exactly 1 nonzero entry.
                batch = np.zeros((batch_size,Ncoupons))#.astype(int)
                batch[np.arange(batch_size),coupons] = 1
    
                #Get cumulative sum over simulations.
                #The collection is completed at the first row at which the cumulative count for all coupons is >= quota for that coupon.
                cumu_batch = np.cumsum(batch,axis=0)
    
                
                #Find first instance where cumulative counts >= quota vector at every element
                #Subtract the quota vector from each row. The first row witch all nonnegative numbers is the row.
                cumu_batch2 = cumu_batch + np.repeat(np.expand_dims(counts_vector,axis=0),batch_size,axis=0) - rep
                check = np.sum(cumu_batch2>=0,axis=1)
                
                #If a row has all nonnegative elements, then all quotas have been met at that row so the collection is complete.
                cond = np.where(check==Ncoupons)[0]
    #            print cond
    #            print cond.size>0
                if cond.size>0:
                    Vals[i] = cond[0] + 1 + Nbatches*batch_size
    #                counts_vector += cumu_batch[cond[0]]
    #                print 'done'
    #                print counts_vector - quota_vector
    #                print 'done'
                    break
                #Otherwise, if collection still not complete, add the counts vector for this batch to the running counts and keep going
                else:
                    counts_vector += cumu_batch[-1]
    #                print counts_vector
                    Nbatches += 1
                
    #        print '\n'*5
            
            
            
    
    
    #Once the simulations are done, get basic statistics:
    MEAN = np.mean(Vals)
    MEDIAN = np.median(Vals)
    STD = np.std(Vals)
    #quartiles ...
    
    #Plot PMF of Vals:
    y = np.bincount(Vals)
    plt.figure(figsize=(18,12))
    plt.title('Unnormalized PMF of #Draws to Complete Collection',fontsize=30)
    plt.xlabel('# Draws',fontsize=30)
    plt.ylabel('# Occurrences',fontsize=30)
    plt.plot(y,marker='o',color='b')
    
    #Plot Histogram of Vals:
    plt.figure(figsize=(18,12))
    plt.title('Histogram of #Draws to Complete Collection',fontsize=30)
    plt.xlabel('# Draws',fontsize=30)
    plt.ylabel('# Occurrences',fontsize=30)
    plt.hist(Vals,bins=100)

    print 'Nsimulations', Nsimulations
    print 'Ncoupons', Ncoupons
    print 'MEAN', MEAN
    print 'MEDIAN', MEDIAN
    print 'STD', STD
    print 'Vals', Vals
    
    return Vals, MEAN, MEDIAN, STD
    
    
    
    
    
    
    
    
    
if __name__ == "__main__":
    
    
#==============================================================================
#     Testing analytical mean for NON-uniform case
#==============================================================================
#    quota_vector = np.array([1,1,2,3,4])
#    probability_vector = np.array([.2, .2, .3, .1, .2])
    
    #To match paper Dr. Pepper example:
    probability_vector = np.array([.25, .25, .15, .35])
    quota_vector = np.array([1,2,3,2])
    #For 2nd example of win "twice":
    #quota_vector *= 2
    
    NonUniExpectedValue = NonUniformCouponCollectorQuotas__AnalyticalMean(quota_vector,probability_vector)
    print 'NonUniExpectedValue', NonUniExpectedValue





#==============================================================================
#     Testing analytical mean for uniform case
#==============================================================================
#    #100 and 100 to match original paper
#    d = 100
#    n = 100
    d = 43
    n = 12
    UNIExpectedValue = UniformCouponCollectorQuotas__AnalyticalMean(d,n)    
    print 'UNIExpectedValue', UNIExpectedValue
    
    
    
    
    
    
    
    
    
#==============================================================================
#     Calculating mean and variance by simulation
#==============================================================================
    Nsimulations = 20000#0
    batch_method = True#False
    
    #Coupon skew example:
    from scipy.stats import halflogistic
    N = 50#200 #Number of Coupons
    quota_vector = np.ones(N)
    x = np.linspace(halflogistic.ppf(0.01), halflogistic.ppf(0.99), N)
    probability_vector = halflogistic.pdf(x)
    probability_vector /= probability_vector.sum()
    plt.figure()
    plt.title('Skew Coupon PDF')
    plt.plot(probability_vector)
    
    t0 = time.clock()
    print 'Doing Simulation...'
    Vals, MEAN, MEDIAN, STD = NonUniformCouponCollectorQuotas__Simulation(quota_vector,probability_vector,Nsimulations,batch_method)
    t1 = time.clock()
    print 'Time Elapsed: ', t1-t0