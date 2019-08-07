import numpy as np
import time
from params import input_params
# import matplotlib.pyplot as plt

class Random_DNA_Strand_Displacement_Circuit(object):
    '''
        Random DNA Strand Displacement Circuit class to generate a random DNA strand diaplacement circuit
        Take parameters (input, etc.) from params.py
            Optimal random DNA chemistry using results from optimization methods (turning point measure and genetic search)
            (nL, nU, k) = (3, 4, 3)
            Global ordered type
    '''

    def __init__(self, input_params={}):
        self.input_params = input_params.copy() # copy the input_params dict from params.py
        nU, nL, U, L = self.initialize_single_strands(self.input_params) # run initialize_single_strands() method
        nF, nP, F, P = self.initialize_double_strands(input_params=self.input_params, nU=nU, nL=nL, U=U, L=L) # run initialize_double_strands() method

    def initialize_single_strands(self, input_params):
        '''
            Method to initialize the set of single DNA strands
            Args:
                input_params (dict): input parameters
            Returns:
                nU (int): number of upper single strands
                nL (int): number of lower single strands
                U (tuple): set of upper single strands
                L (tuple): set of lower single strands
        '''

        nL = int(round(self.input_params['n'] / (1 + self.input_params['p'])))
        nU = self.input_params['n'] - nL
        U = ('U0', 'U1', 'U2', 'U3', 'U4', 'U5', 'U6', 'U7', 'U8', 'U9')
        L = ('L0', 'L1', 'L2', 'L3', 'L4', 'L5', 'L6', 'L7', 'L8', 'L9')
        U = U[0:nU]
        L = L[0:nL]
        return nU, nL, U, L

    def initialize_double_strands(self, input_params, nU, nL, U, L):
        '''
            Method to initialize the set of single DNA strands
            Args:
                input_params (dict): input parameters
                nU (int): number of upper single strands
                nL (int): number of lower single strands
                U (tuple): set of upper single strands
                L (tuple): set of lower single strands
            Returns:
                nF (int): number of full double strands
                nP (int): number of partial double strands
                F (tuple): set of full double strands
                P (tuple): set of partial double strands
        '''
        min_single = min([nL, nU])
        nF = int(self.input_params['y'] * min_single)
        index_full = np.random.choice(a=min_single, size=nF, replace=False)
        F = []
        for i in index_full:
            F.append(U[i] + L[i])

        # positive normal distribution of partial double strands per upper strand
        norm_dist = np.random.normal(loc=self.input_params['phi']['mean'], 
                                     scale=self.input_params['phi']['variance'],
                                     size=nU)
        self.input_params['phi'].update({'norm_dist': norm_dist})

        num_double_per_upper = []
        for i in norm_dist:
            num_double_per_upper.append(int(round(i)))
        if len(set(num_double_per_upper)) == 1:
            num_double_per_upper = num_double_per_upper[0] # k (int): number of double strands per upper strands
        else:
            raise BaseException

        P = F
        while len(set(F + P)) != len(F + P):
            P = []
            for upper in U: # for each upper strand
                # Choose randomly number of lower strand counterparts without repetitions
                lower_per_upper = np.random.choice(a=L, size=num_double_per_upper, replace=False)
                for lower in lower_per_upper:
                    partial = upper + lower
                    P.append(partial)
        nP = len(P)
        F = tuple(F)
        P = tuple(P)
        return nF, nP, F, P

network = Random_DNA_Strand_Displacement_Circuit(input_params)
# print network.input_params['n']