"""
 Created by Jorge Gomes on 02/05/2018
 MReconstruction
 fastcore
 
"""
from optlang.cplex_interface import Variable, Constraint, Model, Objective
from framed import simplify
from framed import load_cbmodel
from reconstruction.configuration import Configuration
from copy import deepcopy
from reconstruction.reconstruction_base import Reconstruction

"""
FASTCORE Algorithm for context specific metabolic network reconstruction
by: Nikos Vlassis, Maria Pires Pacheco, Thomas Sauter

Inputs: 
    - A consistent metabolic network model {N, Sn} where:
        -N is the reaction set {1,2,3, ...} 
        -Sn is the stoichiometric matrix
        * To be consistent a metabolic network model must have no blocked reactions. To accomplish this simplify (model)
        from framed package can be used
    - Core reaction set C that are supported by strong evidence to be active in the desired context such that
        - C contains some but not all of N, thus, a strict subset of N

Outputs:
    - A consistent induced subnetwork {A, Sa} of {N, Sn} such that all of C is included in A (subset)
    
Notes:
    * A mode is a feasible flux vector that satisfies a group of steady-state constraints namely :
        - Sv = 0
        - v (belongs to) B, with B being the lower and upper bound per reaction , i.e. B = [lb, ub]
    * In a mode, the set of nonzero entries are known as the support of that node.
    * An elementary mode is a feasible flux vector v!=0 with minimal support. That is, there is no other mode which has 
    a smaller support set.
    * A reaction i is called blocked if it cannot be active under any mode, that is, there is no mode v, where vi != 0 
    * A dense mode is a mode that contains as many nonzero fluxes as possible.
    * A sparse mode is a mode that contains as many zero fluxes as possible.
Further documentation can be found:
    https://opencobra.github.io/cobratoolbox/stable/modules/dataIntegration/transcriptomics/FASTCORE/index.html
"""


class Fastcore(Reconstruction):

    def __init__(self, model, core_reactions, config, simplified=False, scaling_factor=10**5):
        """
        :param model: framed model object, metabolic model
        :param core_reactions: set, core reactions used as input to the FASTCORE algorithm
        :param config: Configuration object, contains parameters necessary to the algorithm
        :param simplified: boolean, whether the model is consistent (all reactions are able to carry flux -> True)
                           or still needs to be checked for consistency (False). When set to True erroneously results
                           will not be trustworthy.
        :param scaling_factor: numeric, scaling factor to be considered during the algorithm
        """
        Reconstruction.__init__(self, model, core_reactions, config, simplified)
        self._scalingFactor = scaling_factor
        self._swap = False

    def generateModel(self):
        """
        Generates a context specific metabolic network reconstruction based on the FASTCORE algorithm.
        singleton is boolean, a flag indicating whether to check the reversible reactions in the model.
        flipped is another boolean flag that ensures that reversible reactions are checked for both negative and forward
        directions. N is a set, the reaction set that composes a metabolic model. C: set, the core reaction set
        supported by strong evidence. Letters kept as in the paper to ease comparison between pseudo-code and the actual
        python implementation.

        :return: A consistent induced subnetwork {A, Sa} of {N, Sn} such that all of C is included in A (subset), in the
                 shape of a framed object metabolic model.
        """
        print("Consistent input network has N=", len(self.model.reactions.keys()))
        N = set(self.model.reactions.keys())
        C = set(self.core)
        Irrev = self.get_IrrevReactions(N)  # I from pseudo-code
        visited = set()
        J = C & Irrev  # IrrevCore no código da Sara
        P = N - C
        flipped, singleton = False, False
        print('N:{0}, J:{1}, P:{2}'.format(len(N), len(J), len(P)))
        A = self.FindSparseMode(J, P, singleton)
        J = C - A

        while len(J) > 0:
            P = P - A
            A = A | self.FindSparseMode(J, P, singleton)
            self._swap = False
            # print(visited)
            print('J:', len(J), 'A:', len(A))
            if len(J & A) != 0:
                J = J - A
                flipped = False
            else:
                if flipped:
                    flipped = False
                    singleton = True
                else:
                    flipped = True

                    J_line = set()  # J' I_Line no código da Sara
                    if singleton:
                        it = iter(J)
                        elem = ""
                        iter_count = 1
                        while elem == "" and iter_count <= len(J):
                            elem = next(it)
                            if elem in visited:
                                elem = ""
                            iter_count += 1

                        J_line.add(elem)
                        visited.add(elem)

                    else:
                        J_line = J_line | J

                    J_line = J_line - Irrev

                    for reaction in J_line:
                        if len(J_line) == 1 and reaction == "":
                            print('escape route')
                            return self._buildTissueSpecificModel(A)
                        if self.model.reactions[reaction].reversible:
                            self.swapDirectionMatrix(reaction)
                    self._swap = True

        return self._buildTissueSpecificModel(A)

    def FindSparseMode(self, J, P, singleton):
        """
        Called with each iteration of the FASTCORE algorithm this function adds to the set A the support of a mode that
        dense in J and sparse in P.

        :param J: set, a subset of the C (core) reaction set containing only the irreversible reactions in C (irrevCore)
        :param P: set, a penalty set containing the reactions in N that are not included in C
        :param singleton: boolean, a flag indicating whether to check each reaction individually
        :return: set, a set of reactions that compose a sparse mode
        """
        res = set()

        if len(J) > 0:
            if self._swap:  # if reaction signs have been flipped new sv=0 constraints must be created
                self._createSVO(reset=True)
            if singleton:  # if singleton is true only the first reaction of J is checked
                lp7_res = self.MaxNumberReactions(next(iter(J)))
            else:
                lp7_res = self.MaxNumberReactions(J)

            K = self.get_ActiveFluxes(lp7_res, J, bothDir=False)

            if len(K) != 0:
                lp10_res = self.MinimizeFluxPenaltySet(K, P)
                res = self.get_ActiveFluxes(lp10_res, self.model.reactions.keys(), True)

        return res

    def swapDirectionMatrix(self, reaction):
        """
        Flips the signs of all metabolites coefficients on corresponding reaction column, while swapping its lb and ub
        to the inverse values.

        :param reaction: str, reaction id
        """
        reac = self.model.reactions[reaction]  # reduce verbose
        lb = reac.lb
        ub = reac.ub

        reac.lb = ub  # swap lower bound
        if reac.lb is not None:
            reac.lb *= -1
            lb = reac.lb

        reac.ub = lb  # swap upper bound
        if reac.ub is not None:
            reac.ub *= -1
            ub = reac.ub

        # changing the stoichiometric matrix
        for m_id, coef in reac.stoichiometry.items():
            reac.stoichiometry[m_id] = coef * -1

        # changing v_i variable
        var = self._VI[reaction]
        var.set_bounds(lb, ub)

    def MaxNumberReactions(self, J):
        """
        LP-7 formulation from FASTCORE algorithm paper

        Maximize the number of reaction, from a set of reactions, with positive flux rate, i.e., this LP tries to
        maximize the number of feasible fluxes in J (R & I) whose value is at least the flux threshold.

        OF:     Max SUM z_i

        s.t:    Sv = 0
                v_i >= z_i , i in Reacs
                z_i in [0,e] , e threshold
        """

        # CREATE EMPTY PROBLEM
        lp7 = Model(name='LP-7')

        # VARIABLE DEFINITION : z_i and v_i vectors
        v_i = self._createVI(scaled=True, scaling_factor=self._scalingFactor)[0]  # creates the v_i vector of variables

        z_i = {}
        for r_id, r in v_i.items():
            if r_id in J:
                name = 'z_' + r_id
                temp_ub = self.config.flux_threshold()
                z = Variable(name, lb=0, ub=temp_ub)
                z_i[name] = z  # dicionário com as variáveis z_i
                # lp7.add(z_i.values())

        # CONSTRAINT DEFINITION

        # S.V = 0
        for metabolite_con in self._createSVO():
            lp7.add(metabolite_con)

        # v_i - z_i  >= 0
        # contains only z_i for those i contained in J

        for z_id, z in z_i.items():  # variable name (str), variable object (optlang)
            r = v_i[z_id.replace('z_', '')]
            c = Constraint(r - z, lb=0)
            lp7.add(c)

        # OBJECTIVE FUNCTION DEFINITION
        obj = Objective(sum(z_i.values()), direction='max')
        lp7.objective = obj

        # OPTIMIZATION
        lp7.optimize()
        # print('lp7 obj:', lp7.objective.value)
        return lp7.variables.items()

    def MinimizeFluxPenaltySet(self, K, PenaltySet):
        """
        LP-10 formulation from FASTCORE algorithm paper

        Minimize the number of reaction from penalty set, when the K set must have flux.

        OF:     Min SUM Z_i

        s.t:    Sv = 0
                v_i >= e , i in K
                vi in [-zi, zi], i in PenaltySet

        :param K: Active subset of J reaction set , computed by the LP-7 formulation.
        :param PenaltySet: A set of reactions that contains all reactions outside of C (core), that have not been added
                           yet to the set A. Using Python's set syntax : P = (N - C) - A
        :return:
        """

        # CREATE EMPTY PROBLEM
        lp10 = Model(name='LP-10')

        # VARIABLE DEFINITION : z_i and v_i vectors
        v_i = self._createVI(scaling_factor=self._scalingFactor, scaled=True)[0]  # creates v_i vector of variables

        z_i = {}  # for reactions in the penalty set

        for r_id, r in v_i.items():
            # z_i
            if r_id in PenaltySet:
                # lower bound for Z_i in matlab = max(abs(model.ub(P)),abs(model.lb(P)))
                reac = self.model.reactions[r_id]
                name = 'z_' + r_id
                if reac.lb is not None and reac.ub is not None:
                    Zub = max(abs(reac.lb), abs(reac.ub)) * self._scalingFactor
                else:
                    Zub = None
                z = Variable(name, lb=0, ub=Zub)  # z belongs to R+, TRY NONE MAYBE?
                z_i[name] = z  # dicionário com as variáveis z_i
                # lp7.add(z_i.values())

        # CONSTRAINT DEFINITION
        # S.v = 0
        for metabolite_con in self._createSVO():
            lp10.add(metabolite_con)

        # v_i in [-z_i, z_i], for i in P
        for z_id, z in z_i.items():  # z_i variables only exists for i's included in P
            r = v_i[z_id.replace('z_', '')]
            con1 = Constraint(r - z, ub=0)
            con2 = Constraint(r + z, lb=0)
            lp10.add([con1, con2])

        # v_i >= threshold , for i in K
        for reaction in K:
            r = v_i[reaction]
            con = Constraint(r - (self._scalingFactor * self.config.flux_threshold()), lb=0)  # scaled threshold
            lp10.add(con)

        # OBJECTIVE DEFINITION
        obj = Objective(sum(z_i.values()), direction='min')
        lp10.objective = obj

        # TODO: Create additional v_is for lp10 since all bounds are scaled according to matlabs implementation

        # OPTIMIZATION
        lp10.optimize()
        # print('lp10 obj:', lp10.objective.value)
        return lp10.variables.items()

    def get_ActiveFluxes(self, optResults, reacs, bothDir=True):
        res = set()
        optRes = {x: y.primal for x, y in optResults}

        for r in reacs:
            if optRes[r] >= 0.99 * self.config.flux_threshold() or \
               (bothDir and optRes[r] <= -1.0 * 0.99 * self.config.flux_threshold()):
                res.add(r)
        return res


if __name__ == '__main__':
    example = 'C:/Users/Tese_Avoid_Namespaces/Tese/MReconstruction/models/Ec_iAF1260_flux1.xml'
    print('loading model')
    model1 = load_cbmodel(example)
    config1 = Configuration()
    fast = Fastcore(model1, [list(model1.reactions.keys())[x] for x in range(len(model1.reactions.keys())) if x < 100], config1)
    fast.swapDirectionMatrix('R_12DGR120tipp')

