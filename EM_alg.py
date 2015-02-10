"""
A python implementation of the expectation - maximization algorithm (EM)
example in the Primer (Nature Biotechnology)

What is the expectation maximization algorithm?
Chuong B Do & Serafim Batzoglou

Link: http://www.nature.com/nbt/journal/v26/n8/full/nbt1406.html

@original author: Saeed Al Turki (salturki at gmail.com)
@modified by: Yusan Lin (yusan@psu.edu)
@modifications: 1) changed from terminal-friendly to Python IDLE-friendly
                2) changed from fixed number of iterations to detect convergence
"""
import math
import random # for flipping coins

def countHeads(experiment):
    return experiment.count('H')

def binom(trials, x, p):
    return (math.factorial(trials) / (math.factorial(x) * math.factorial(trials - x))) * (p ** x) * ((1-p) ** (trials - x))

def stepE(experiments, theta_A, theta_B):
    """
    Compute a probability distribution over possible completions using the current parameters
    """
    completions = {}
    for nbr in sorted(experiments):
        experiment  = experiments[nbr]
        trials      = len(experiment)
        if not completions.has_key(nbr):
            completions[nbr] = {'theta_A':[], 'theta_B':[]}
        headsCount  = countHeads(experiment)
        tailsCount  = trials - headsCount
        binom_A     = binom(trials, headsCount, theta_A)
        binom_B     = binom(trials, headsCount, theta_B)
        binom_total = binom_A + binom_B
        norm_A      = binom_A / binom_total
        norm_B      = binom_B / binom_total
        completions[nbr]['theta_A'] = [norm_A * headsCount, norm_A * tailsCount] 
        completions[nbr]['theta_B'] = [norm_B * headsCount, norm_B * tailsCount]
    return completions


def stepM(completions):
    """
    new parameters are determined using the current completions
    """
    updated_theta_A = 0.
    updated_theta_B = 0.
    heads_A = 0.
    tails_A = 0.
    heads_B = 0.
    tails_B = 0.
    for nbr in completions:
        heads_A += completions[nbr]['theta_A'][0]
        tails_A += completions[nbr]['theta_A'][1]
        heads_B += completions[nbr]['theta_B'][0]
        tails_B += completions[nbr]['theta_B'][1]

    updated_theta_A = heads_A / (heads_A+tails_A)
    updated_theta_B = heads_B / (heads_B+tails_B)
    return updated_theta_A, updated_theta_B


def generateData(theta_A, theta_B, n_observations, n_flips, pA):
    experiments = dict()
        
    for n in range(n_observations):
        
        # pick coin
        coin = pickCoin(pA)
        if coin == 'A':
            theta = theta_A
        else:
            theta = theta_B

        # flip coin            
        flips = ""
        for i in range(n_flips):
            flips += flipCoin(theta)
        experiments[n] = flips
        
    return experiments


def flipCoin(prob_H):
    randRoll = random.random() # random number [0,1)
    if randRoll < prob_H:
        return 'H'
    else:
        return 'T'


def pickCoin(pA):
    randPick = random.random()
    if randPick < pA:
        return 'A'
    else:
        return 'B'


def run(pA = 0.5, theta_A_real = 0.5, theta_B_real = 0.5, n_observations = 1000, n_flips = 100):

    # generate data
    experiments = generateData(theta_A_real, theta_B_real, n_observations, n_flips, pA)

    # initialization
    theta_A = random.random()
    theta_B = random.random()

    # iterations and convergence
    iterations = 0
    max_iterations = 50
    eps = 0.0001
    theta_A_old, theta_B_old = float("inf"), float("inf")

    header = ['t','Theta_A', 'Theta_B']
    print "\t".join(header)

    # test for convergence and limit the maximum number of convergence
    while ((abs(theta_A - theta_A_old) > eps) or (abs(theta_B - theta_B_old) > eps))\
          and (iterations < max_iterations):
        
        theta_A_old = theta_A
        theta_B_old = theta_B
        
        print "%d\t%.3f\t%.3f" % (iterations, theta_A, theta_B)
        completions      = stepE(experiments, theta_A, theta_B)
        theta_A, theta_B = stepM(completions)
        iterations += 1
