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
        self.impose_order_partial(P=P, nP=nP, F=F, nF=nF)
        nS, nSS, nDS, S, SS, DS = self.create_species_set(U=U, nU=nU, L=L, nL=nL, F=F, nF=nF, P=P, nP=nP)
        nI, nO, I, O = self.create_influx_efflux(input_params=self.input_params, S=S, nS=nS)

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

    def mirror_selection(self):
        pass

    def impose_order_partial(self, P, nP, F, nF):
        '''
            Method to create a lookup dictionary to show ordering of the double strands
            Args:
                P (tuple): set of partial double strands
                nP (int): number of partial double strands
                F (tuple): set of full double strands
                nF (int): number of full double strands
            Returns:
                order_lookup (dict): lookup dictionary to show double strands' ordering
        '''

        order_lookup = {}
        for p_index in range(nP):
            partial_with_order = np.random.choice(a=P)
            while partial_with_order in order_lookup.values():
                partial_with_order = np.random.choice(a=P)
            order_lookup.update({'{}'.format(p_index): partial_with_order})
        if nF <= 1:
            order_lookup.update({'max': F[0]})
        else:
            order_lookup.update({'max': F})
        return order_lookup

    def create_species_set(self, U, nU, L, nL, F, nF, P, nP):
        '''
            Method to initialize the set of single DNA strands
            Args:
                nU (int): number of upper single strands
                nL (int): number of lower single strands
                U (tuple): set of upper single strands
                L (tuple): set of lower single strands
                nF (int): number of full double strands
                nP (int): number of partial double strands
                F (tuple): set of full double strands
                P (tuple): set of partial double strands
            Returns:
                nS (int): number of all DNA strands in network
                nSS (int): number of all single DNA strands in network
                nDS (int): number of all double DNA strands in network
                S (tuple): set of all DNA strands in network
                SS (tuple): set of all single DNA strands in network
                DS (tuple): set of all double DNA strands in network
        '''

        S = U + L + F + P
        nS = len(S)
        SS = U + L
        nSS = len(SS)
        DS = F + P
        nDS = len(DS)

        if nS != (nU + nL + nF + nP):
            raise BaseException
        if nSS != (nU + nL):
            raise BaseException
        if nDS != (nF + nP):
            raise BaseException

        return nS, nSS, nDS, S, SS, DS

    def create_influx_efflux(self, input_params, S, nS):
        nI = int(round(self.input_params['a_in'] * nS))
        nO = int(round(self.input_params['a_out'] * nS))
        I = np.random.choice(a=S, size=nI, replace=False)
        O = np.random.choice(a=S, size=nO, replace=False)
        return nI, nO, I, O

    def create_reaction_set(self):
        pass

network = Random_DNA_Strand_Displacement_Circuit(input_params)