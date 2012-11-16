#!/usr/bin/env python
"""
Utilities for SHA python version 
""" 
# Needed libraries:
import os, sys
import numpy as np 

from time import sleep 
import MySQLdb as mdb

from my_util.functools import * 
from my_util.numer import * 

from pynga.utils import * 
from pynga import * 


# commonly used functions
def GaussianDistribution(x,mu=0,sigma=1):
    A = x-mu
    B = 1./(sigma*np.sqrt(2.*np.pi))
    C = np.exp( -0.5*(A/sigma)**2. )
    return B*C

def GaussianCumulative(LowerLimit, UpperLimit, mu,std, N = 100, met='Simpson'): 
    h = (UpperLimit-LowerLimit) / (N-1)
    x = np.arange( N ) * h + LowerLimit
    ft = GaussianDistribution(x, mu=mu, sigma=std) 
    integral = time_integral( ft, LowerLimit, UpperLimit, h, met=met )
    return integral 





