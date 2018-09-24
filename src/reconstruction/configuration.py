"""
 Created by Jorge Gomes on 21/05/2018
 MReconstruction
 configuration
 
"""


class Configuration:
    """
    Configuration class to provide parameters for different metabolic model reconstruction algorithms
    """

    def __init__(self, lb=-1000, ub=1000, flux_threshold=1*10**-4, solver='cplex', ):
        """
        :param lb: numeric, lower bound
        :param ub: numeric, upper bound
        :param flux_threshold: numeric, threshold to be considered by optimization solver
        :param solver: string, solver to be used
        """
        self._lb = lb
        self._ub = ub
        self._flux_threshold = flux_threshold
        self._solver = solver
        if not self._validateConfig():
            raise TypeError('Configuration parameters are not valid')

    def lb(self):
        return self._lb

    def ub(self):
        return self._ub

    def flux_threshold(self):
        return self._flux_threshold

    def solver(self):
        return self._solver

    def set_lb(self, lb):
        if self._validateConfig():
            self._lb = lb

    def set_ub(self, ub):
        if self._validateConfig():
            self._ub = ub

    def set_fluxThreshold(self, ft):
        if self._validateConfig():
            self._flux_threshold = ft

    def set_solver(self, solver):
        if self._validateConfig():
            self._solver = solver

    def _validateConfig(self):
        if type(self._lb) not in (int, float):
            print('Lower bound must be numeric')
            return False

        if type(self._ub) not in (int, float):
            print('Upper bound must be numeric')
            return False

        if self._ub <= self._lb:
            print('Upper bound must always be greater than the lower bound')
            return False

        if type(self._flux_threshold) not in (int, float):
            print('Flux Threshold must be numeric')
            return False

        if self._solver.lower() not in ('cplex', 'glpk'):
            print('Solver must be either \'cplex\' or \'glpk\'')
            return False

        return True

if __name__ == '__main__':
    c = Configuration('ss')
    print(c.lb())