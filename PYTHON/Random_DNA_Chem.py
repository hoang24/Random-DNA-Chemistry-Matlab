import numpy as np
import time


class Random_DNA_Strand_Displacement_Circuit(object):
    '''
        Random DNA Strand Displacement Circuit class to generate a random DNA strand diaplacement circuit
        Take parameters (input, etc.) from params.py
            Optimal random DNA chemistry using results from optimization methods (turning point measure and genetic search)
            (nL, nU, k) = (3, 4, 3)
            Global ordered type
    '''

    def __init__(self, input_params={}):
        '''
            Args:
                input_params (dict): input parameters
        '''

        self.input_params = input_params.copy() # copy the input_params dict from params.py

        nU, nL, U, L = self.create_species_single() # run create_species_single() method
        
        nF, nP, F, P = self.create_species_double(nU=nU, nL=nL, U=U, L=L) # run create_species_double() method

        self.mirror_selection()

        order_lookup = self.impose_order_partial(P=P, nP=nP, F=F, nF=nF)
        self.order_lookup = order_lookup

        nS, nSS, nDS, S, SS, DS = self.create_all_species_sets(U=U, nU=nU, L=L, nL=nL, F=F, nF=nF, P=P, nP=nP)

        nI, nO, I, O = self.create_influx_efflux(S=S, nS=nS)
        self.species_lookup = {
            'U': U,
            'L': L,
            'F': F,
            'P': P,
            'SS': SS,
            'DS': DS,
            'S': S,
            'I': I,
            'O': O,
            'nU': nU,
            'nL': nL,
            'nF': nF,
            'nP': nP,
            'nSS': nSS,
            'nDS': nDS,
            'nS': nS,
            'nI': nI,
            'nO': nO
        }

        R_DISPLACE, nR_DISPLACE = self.create_reaction_displacement(SS=SS, order_lookup=order_lookup)
        self.R_DISPLACE = R_DISPLACE
        R_BIND, nR_BIND = self.create_reaction_binding(U=U, L=L, DS=DS)
        self.R_BIND = R_BIND
        R_IN, nR_IN = self.create_reaction_influx(I=I)
        R_OUT, nR_OUT = self.create_reaction_efflux(O=O)
        R, nR = self.create_all_reaction_set(R_BIND=R_BIND, R_DISPLACE=R_DISPLACE, R_IN=R_IN, R_OUT=R_OUT, 
                                             nR_BIND=nR_BIND, nR_DISPLACE=nR_DISPLACE, nR_IN=nR_IN, nR_OUT=nR_OUT)
        self.reaction_lookup = {
            'R_DISPLACE': R_DISPLACE,
            'R_BIND': R_BIND,
            'R_IN': R_IN,
            'R_OUT': R_OUT,
            'R': R,
            'nR_DISPLACE': nR_DISPLACE,
            'nR_BIND': nR_BIND,
            'nR_IN': nR_IN,
            'nR_OUT': nR_OUT,
            'nR': nR
        }

        ic_U, ic_L, ic_F, ic_P, IC = self.create_initial_concentration(nU=nU, nL=nL, nF=nF, nP=nP)
        self.concentration_lookup = IC

        k_BIND = self.create_rateConst_binding(nR_BIND)
        k_DISPLACE = self.create_rateConst_displacement(nR_DISPLACE)
        k_IN = self.create_rateConst_influx(nR_IN)
        k_OUT = self.create_rateConst_efflux(nR_OUT)
        self.rateConst_lookup = {
            'k_BIND': k_BIND,
            'k_DISPLACE': k_DISPLACE,
            'k_IN': k_IN,
            'k_OUT': k_OUT
        }

    def create_species_single(self):
        '''
            Method to initialize the set of single DNA strands
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

    def create_species_double(self, nU, nL, U, L):
        '''
            Method to initialize the set of single DNA strands
            Args:
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

        min_single = np.min([nL, nU])
        nF = int(self.input_params['y'] * min_single)
        index_full_list = np.random.choice(a=min_single, size=nF, replace=False)
        F = []
        for index_full in index_full_list:
            F.append(U[index_full] + L[index_full])

        # positive normal distribution of partial double strands per upper strand
        norm_dist = np.random.normal(loc=self.input_params['phi']['mean'], 
                                     scale=self.input_params['phi']['variance'],
                                     size=nU)
        self.input_params['phi'].update({'norm_dist': norm_dist})

        num_double_per_upper = []
        for norm_dist_value in norm_dist:
            num_double_per_upper.append(int(round(norm_dist_value)))
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

    def create_all_species_sets(self, U, nU, L, nL, F, nF, P, nP):
        '''
            Method to initialize the all-species sets
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

    def create_influx_efflux(self, S, nS):
        '''
            Methods to create set of influx and efflux
            Args:
                S (tuple): set of all species in network
                nS (int): number of all species in network
            Returns:
                nI (int): number of influx species
                nO (int): number of efflux species
                I (tuple): set of influx species
                O (tuple): set of efflux species
        '''
        nI = int(round(self.input_params['a_in'] * nS))
        nO = int(round(self.input_params['a_out'] * nS))
        I = np.random.choice(a=S, size=nI, replace=False)
        O = np.random.choice(a=S, size=nO, replace=False)
        return nI, nO, I, O

    def _get_key_by_value_from_dict(self, dictionary, value):
        '''
            Method to retrieve a key given its value in a dictionary. Apply to dictionary with non-repeating values.
            Args:
                dictionary (dict): a dictionary
                value (int, float, string, etc): value in dictionary
            Returns:
                key (string): key of the value in the dictionary
        '''
        item_list = dictionary.items()
        key = -1
        for item in item_list:
            if value == item[1]:
                key = item[0]
        if key == -1:
            raise Exception('Value given does not exist in given dictionary. Value = {}'.format(value))
        return key

    def create_reaction_binding(self, U, L, DS):
        '''
            Method to create a set of binding reactions
            Args:
                U (tuple): set of upper strand species
                L (tuple): set of lower strand species
                DS (tuple): set of all double strand species
            Returns:
                R_BIND (tuple): set of binding reactions
                nR_BIND (int): number of binding reactions
        '''

        R_BIND = []
        for u in U:
            for l in L:
                ds_result = u + l
                if ds_result in DS:
                    R_bind = u + ' + ' + l + ' --> ' + ds_result
                    R_BIND.append(R_bind)
        R_BIND = tuple(R_BIND)
        nR_BIND = len(R_BIND)
        return R_BIND, nR_BIND

    def create_reaction_displacement(self, SS, order_lookup):
        '''
            Method to create a set of displacement reactions
            Args:
                SS (tuple): set of all single strand species
                order_lookup (dict): lookup dictionary showing all double strand species and their ordering
            Returns:
                R_DISPLACE (tuple): set of displacement reactions
                nR_DISPLACE (int): number of displacement reactions
        '''

        R_DISPLACE = []
        for reactant_ss in SS: # for each single strand
            for reactant_ds_order, reactant_ds in order_lookup.items(): # for each double strand in order
                if reactant_ss not in reactant_ds: # if single strand is not part of double strand
                    if 'U' in reactant_ss: # if single strand is upper strand
                        result_ds = reactant_ss + reactant_ds[-2:]
                        result_ss = reactant_ds[:2]
                    elif 'L' in reactant_ss: # if single strand is lower strand
                        result_ds = reactant_ds[:2] + reactant_ss
                        result_ss = reactant_ds[-2:]
                    else: # invalid single strand notation
                        raise Exception('Unknown DNA species type. Defined type: U for Upper strand, L for Lower strand.')
                    if result_ds in order_lookup.values(): # if the resulted double strand made out of the single strand and the reactant double strand exists
                        result_ds_order = self._get_key_by_value_from_dict(order_lookup, result_ds)
                        if result_ds_order > reactant_ds_order: # if the resulted double strand has higher order than the reactant double strand
                            R_displace = reactant_ss + ' + ' + reactant_ds + ' --> ' + result_ds + ' + ' + result_ss
                            R_DISPLACE.append(R_displace)
        R_DISPLACE = tuple(R_DISPLACE)
        nR_DISPLACE = len(R_DISPLACE)
        return R_DISPLACE, nR_DISPLACE

    def create_reaction_influx(self, I):
        '''
            Method to create a set of influx reactions
            Args:
                I (tuple): set of influx species
            Returns:
                R_IN (tuple): set of influx reactions
                nR_IN (int): number of influx reactions
        '''
        R_IN = []
        for i in I:
            R_in = '0 --> ' + i
            R_IN.append(R_in)
        R_IN = tuple(R_IN)
        nR_IN = len(R_IN)
        return R_IN, nR_IN

    def create_reaction_efflux(self, O):
        '''
            Method to create a set of efflux reactions
            Args:
                O (tuple): set of efflux species
            Returns:
                R_OUT (tuple): set of efflux reactions
                nR_OUT (int): number of efflux reactions
        '''
        R_OUT = []
        for o in O:
            R_out = o + ' --> 0'
            R_OUT.append(R_out)
        R_OUT = tuple(R_OUT)
        nR_OUT = len(R_OUT)
        return R_OUT, nR_OUT

    def create_all_reaction_set(self, R_BIND, R_DISPLACE, R_IN, R_OUT, nR_BIND, nR_DISPLACE, nR_IN, nR_OUT):
        '''
            Method to create set of all reactions
            Args:
                R_BIND (tuple): set of binding reactions
                nR_BIND (int): number of binding reactions
                R_DISPLACE (tuple): set of displacement reactions
                nR_DISPLACE (int): number of displacement reactions
                R_IN (tuple): set of influx reactions
                nR_IN (int): number of influx reactions
                R_OUT (tuple): set of efflux reactions
                nR_OUT (int): number of efflux reactions
            Returns:
                R (tuple): set of all reactions
                nR (int): total number of reactions in network
        '''
        R = R_BIND + R_DISPLACE + R_IN + R_OUT
        nR = len(R)
        if nR != (nR_BIND + nR_DISPLACE + nR_IN + nR_OUT):
            raise BaseException
        return R, nR

    def create_initial_concentration(self, nU, nL, nF, nP):
        '''
            Method to create sets of initial concentration for each DNA strand types.
            Args:
                nU (int): number of upper single strands
                nL (int): number of lower single strands
                nF (int): number of full double strands
                nP (int): number of partial double strands
            Returns:
                ic_U (list): initial concentrations for all upper single strands
                ic_L (list): initial concentrations for all lower single strands
                ic_F (list): initial concentrations for all full double strands
                ic_P (list): initial concentration for all partial double strands
                IC (list): initial concentrations for all species, with this order: upper --> lower --> full --> partial
        '''
        ic_U = list(np.random.randint(0, 1000, nU))
        ic_L = list(np.random.randint(0, 1000, nL))
        ic_F = list(np.random.randint(0, 1000, nF))
        ic_P = list(np.random.randint(0, 1000, nP))
        IC = ic_U + ic_L + ic_F + ic_P
        ic_U = tuple(ic_U)
        ic_L = tuple(ic_L)
        ic_F = tuple(ic_F)
        ic_P = tuple(ic_P)
        IC = tuple(IC)
        return ic_U, ic_L, ic_F, ic_P, IC

    def create_rateConst_binding(self, nR_BIND):
        '''
            Method to create sets reaction rate constants for all binding reactions.
            Args:
                nR_BIND (int): number of binding reactions
            Returns:
                k_BIND (list): list of reaction rate constants for all binding reactions
        '''
        norm_dist_bind = np.random.normal(loc=self.input_params['theta']['mean'], 
                                          scale=np.sqrt(self.input_params['theta']['variance']),
                                          size=nR_BIND)
        self.input_params['theta'].update({'norm_dist_bind': norm_dist_bind})
        k_BIND= []
        for norm_dist_bind_value in norm_dist_bind:
            k_BIND.append(np.abs(norm_dist_bind_value))
        k_BIND = tuple(k_BIND)
        return k_BIND

    def create_rateConst_displacement(self, nR_DISPLACE):
        '''
            Method to create sets reaction rate constants for all displacement reactions.
            Args:
                nR_DISPLACE (int): number of displacement reactions
            Returns:
                k_DISPLACE (list): list of reaction rate constants for all displacement reactions
        '''
        norm_dist_displace = np.random.normal(loc=self.input_params['theta']['mean'], 
                                              scale=np.sqrt(self.input_params['theta']['variance']),
                                              size=nR_DISPLACE)
        self.input_params['theta'].update({'norm_dist_displace': norm_dist_displace})
        k_DISPLACE = []
        for norm_dist_displace_value in norm_dist_displace:
            k_DISPLACE.append(np.abs(norm_dist_displace_value))
        k_DISPLACE = tuple(k_DISPLACE)
        return k_DISPLACE

    def create_rateConst_influx(self, nR_IN):
        '''
            Method to create sets reaction rate constants for all influx reactions.
            Args:
                nR_IN (int): number of influx reactions
            Returns:
                k_IN (list): list of reaction rate constants for all influx reactions
        '''
        norm_dist_in = np.random.normal(loc=self.input_params['theta_in']['mean'], 
                                        scale=np.sqrt(self.input_params['theta_in']['variance']),
                                        size=nR_IN)
        self.input_params['theta_in'].update({'norm_dist_in': norm_dist_in})
        k_IN = []
        for norm_dist_in_value in norm_dist_in:
            k_IN.append(np.abs(norm_dist_in_value))
        k_IN = tuple(k_IN)
        return k_IN

    def create_rateConst_efflux(self, nR_OUT):
        '''
            Method to create sets reaction rate constants for all efflux reactions.
            Args:
                nR_OUT (int): number of efflux reactions
            Returns:
                k_OUT (list): list of reaction rate constants for all efflux reactions
        '''
        norm_dist_out = np.random.normal(loc=self.input_params['theta_out']['mean'], 
                                         scale=np.sqrt(self.input_params['theta_out']['variance']),
                                         size=nR_OUT)
        self.input_params['theta_out'].update({'norm_dist_out': norm_dist_out})
        k_OUT = []
        for norm_dist_out_value in norm_dist_out:
            k_OUT.append(np.abs(norm_dist_out_value))
        k_OUT = tuple(k_OUT)
        return k_OUT

    def run_chem(self):
        '''
            Method to run the random DNA strand displacement circuit chemistry
            Args:
                
            Returns:
                
        '''
        pass

    def Gillespie_initialization(self):
        '''
            Method for Gillespie algorithm's initialization step
            Args:
                
            Returns:
                
        '''
        pass

    def Gilespie_monte_carlo(self):
        '''
            Method for Gillespie algorithm's Monte Carlo step
            Args:
                
            Returns:
                
        '''
        pass

    def Gillespie_update(self):
        '''
            Method for Gillespie algorithm's update step
            Args:
                
            Returns:
                
        '''
        pass

    def Gillespie_iterate(self):
        '''
            Method for Gillespie algorithm's iterate step
            Args:
                
            Returns:
                
        '''
        pass
