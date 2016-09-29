import numpy as np
import scipy.constants as const



class PlasmaCalculator():
    '''
    PlasmaCalculator provides means to evaluate PIC simulation parameters for a given plasma
    '''

    coulomb_log = 10

    species_dict = {}

    def __init__(self, m_i, T_i, n_i, m_e = const.m_e, T_e = None, n_e = None, B=None):

        '''

        Initializes Plasma Calculator for electron -ion plasma

        :param m_i:
        :param m_e:
        :param T_i:
        :param T_e:
        '''

        if T_e is None:
            T_e = T_i

        if n_e is None:
            n_e = n_i

        #initialize species dictionary
        ions = {'temperature': T_i, 'mass' : m_i, 'density': n_i}
        electrons = {'temperature': T_e, 'mass': m_e, 'density': n_e}

        self.species_dict = {'electrons': electrons, 'ions': ions}

        #calculate plasma paramaters
        self.calc_plasma_potential()

    def calc_plasma_potential(self):

        '''
        Returns plasma potential for a two species electron ion plasma
        :return:
        '''

        m_i = self.species_dict['ions']['mass']
        T_i_K = self.species_dict['ions']['temperature']

        m_e = self.species_dict['electrons']['mass']
        T_e_K = self.species_dict['electrons']['temperature']

        #calculate logarithmic factor for plasma potential
        log_factor = np.log(1. / 2. / np.pi * m_i / const.electron_mass / (1 + T_i_K / T_e_K))

        self.plasma_potentical = 0.5 * const.k * T_e_K / const.e * log_factor

        return self.plasma_potentical

    def calc_debye_lenght(self):

        '''
        Calculate Debye length of two species plasma
        :return:
        '''

        m_i = self.species_dict['ions']['mass']
        T_i_K = self.species_dict['ions']['temperature']

        m_e = self.species_dict['electrons']['mass']
        T_e_K = self.species_dict['electrons']['temperature']
        n_e = self.species_dict['electrons']['density']

        self.debye_lenght = np.sqrt(const.k * T_e_K * const.epsilon_0 / const.e ** 2 / n_e)

        return self.debye_lenght

    def collision_factor(self, m, T, n):
        return const.epsilon_0 ** 2 * const.k * T * np.sqrt(m * const.k * T) / n / const.e ** 4 / self.coulomb_log




