#!/usr/bin/env python
# -*- coding: utf-8 -*-
import numpy as np
import matplotlib as mat
mat.use("agg")
import matplotlib.pyplot as plt


def powerLaw(y, x):
    """
    'When the frequency of an event varies as power of some attribute of that
    event the frequency is said to follow a power law.' (wikipedia)
    This is represented by the following equation, where c and alpha are
    constants:
    y = c . x ^ alpha
    Args
    --------
    y: array with frequency of events >0
    x: numpy array with attribute of events >0
    Output
    --------
    (c, alpha)
    c: the maximum frequency of any event
    alpha: defined by (Newman, 2005 for details):
        alpha = 1 + n * sum(ln( xi / xmin )) ^ -1
    """
    c = 0
    alpha = .0

    if len(y) and len(y)==len(x):
        c = max(y)
        xmin = float(min(x))
        alpha = 1 + len(x) * pow(sum(np.log(x/xmin)),-1)

    return (c, alpha)


def plotPowerLaws(y, x, c=[], alpha=[]):
    """
    Plots the relationship between x and y and a fitted power law on LogLog
    scale.
    Args
    --------
    y: array with frequency of events >0
    x: array with attribute of events >0
    c: array of cs for various power laws
    alpha: array of alphas for various power laws
    """
    plt.figure()
    plt.loglog()
    plt.plot(x,
             y,
             'r+')
    for _c, _alpha in zip(c,alpha):
        plt.plot( (1, max(x)),
                  (_c, _c * pow(max(x), _alpha)),
                  label='~x^%.2f' % _alpha)
        plt.legend()
    plt.savefig("test.powerlaw.png")




if __name__ == '__main__':
    """
    Checking Zipfs law, where the frequency and rank of a word follow a
    specific power law, using the nltk genesis text in english.
    """
    from collections import defaultdict
    import numpy as np

    from nltk.corpus import genesis
    words = genesis.raw('english-web.txt').split(' ')
    y = defaultdict(int)
    for i in words:
        y[i]+=1
    y = sorted(y.values(),reverse=True)
    x = np.array(xrange(1,len(y)+1))

    c, alpha = powerLaw(y, x)
    print 'According to Zipfs law %.2f should be close to 1.' % alpha
plotPowerLaws(y, x, [c,c], [-1,-alpha])
