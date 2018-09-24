"""
 Created by Jorge Gomes on 01/08/2018
 MReconstruction
 reconstruction
 
"""
from framed import simplify
from copy import deepcopy
from optlang.cplex_interface import Variable, Constraint, Model, Objective
from reconstruction.configuration import Configuration


class Reconstruction:
    """
    This class is intended to work as a base class (super) for any network reconstruction algorithm. In it
    common attributes and methods to all network reconstruction algorithm classes' are stored. This is done to simplify
    any updates/ upgrades to the current framework. Some methods are just placeholders to illustrate the desirable
    design of the subclasses.
    """

    def __init__(self, model, core_reactions, config, simplified=False):
        """
        Args:
            model: framed model, metabolic network model loaded with the framed package
            core_reactions: set, set or list containing the core reactions that will be used for model reconstruction
            config: Configuration object, contains parameters necessary to the algorithm
            simplified: boolean, whether the model is consistent (all reactions are able to carry flux -> True)
                        or still needs to be checked for consistency (False). When set to True erroneously results will
                        not be trustworthy.
        """
        if simplified:
            self.model = model
        else:
            self.model = simplify(model, inplace=True)
        self.core = self.get_ConsistentCore(core_reactions)
        self.config = config
        self._SVO = None
        self._VI = None
        self._scaledVI = None

    def generateModel(self):
        """
        @Placeholder method, to be overwritten by the actual algorithm of each method.

        Generates a context specific metabolic network reconstruction based on a specific reconstruction algorithm.

        Returns: a context specific metabolic network reconstruction stored on a model object from framed.
        """
        reactionSet = set()

        return self._buildTissueSpecificModel(reactionSet)

    def _buildTissueSpecificModel(self, reactionSet):
        """
        Builds a tissue-specific metabolic model from the reaction set that is outputted from a given algorithm

        :param reactionSet: set, a set of reactions based on which the original model will be trimmed
        :return: model: framed model object, tissue specific metabolic model
        """
        tsModel = deepcopy(self.model)
        toRemove = list(set(self.model.reactions.keys()) - reactionSet)

        tsModel.remove_reactions(toRemove)

        return tsModel

    def get_IrrevReactions(self, reactions=None):
        """
        Returns from a set of reactions those which are not reversible, or if no set is prompted all irreversible
        reactions contained in the metabolic model

        :param reactions: a set of reactions, with both reversible and irreversible reactions
        :return: a filtered set of reactions containing only the irreversible reactions
        """
        if reactions is None:
            reactions = self.model.reactions.values()
        return set([r for r in reactions if not self.model.reactions[r].reversible])

    def get_ConsistentCore(self, core):
        """
        Removes from the core set reactions that are unable to carry flux.

        :param core: set, core reactions that are input to a given algorithm
        :return: set, a filtered reaction set from the core that are able to carry flux
        """
        return set([reaction for reaction in core if reaction in self.model.reactions.keys()])

    def _createVI(self, scaled=False, scaling_factor=1):
        """
        Create the V_i variables for each reaction for later addition to the optimization problem. The storage
        of these variables as an attribute avoids its creation with every iteration of the algorithm, thus
        improving its performance.

        scaled: boolean, whether or not to use a scaling factor, which is useful with some algorithms
        scaling factor:

        Returns: dictionary, variable name(string): variable object(optlang)

        """
        if self._VI is None:
            v_i = {}
            for r_id, r in self.model.reactions.items():  # be it reaction_id(string), reaction object(framed reaction)
                v = Variable(r_id, lb=r.lb, ub=r.ub)
                v_i[r_id] = v

            self._VI = v_i

            if scaled:
                scaledV_i = deepcopy(v_i)
                for var in scaledV_i.values():
                    if var.lb is not None:
                        nlb = var.lb * scaling_factor
                    else:
                        nlb = None
                    if var.ub is not None:
                        nub = var.ub * scaling_factor
                    else:
                        nub = None
                    var.set_bounds(lb=nlb, ub=nub)
                self._scaledVI = scaledV_i

                return self._VI, self._scaledVI

        else:
            if scaled:
                return self._VI, self._scaledVI
            else:
                return self._VI

    def _createSVO(self, reset=False):
        """
        Create the S.V = 0 restrictions for each metabolite for later addition to the optimization problem. The storage
        of these constraints as an attribute avoids its creation every iteration of the algorithm, thus
        improving its performance.

        reset: boolean, indicates if the restrictions should be created again. Useful when reactions are changed.

        Returns: list, list of Restriction objects(optlang)

        """

        if self._SVO is None or reset:
            SVO = []
            for m_id in self.model.metabolites:
                row = []
                for r_id, reaction in self.model.reactions.items():
                    if m_id in reaction.stoichiometry:
                        row.append(reaction.stoichiometry[m_id] * self._VI[r_id])
                con_row = Constraint(sum(row), lb=0, ub=0)
                SVO.append(con_row)

            self._SVO = SVO
            return self._SVO

        else:
            return self._SVO