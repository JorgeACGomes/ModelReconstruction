"""
 Created by Jorge Gomes on 05/04/2018
 TsmRec
 integrate
 
"""
from framed.omics.simulation import gene_to_reaction_expression as gene2reaction
from omics.omics_data_map import OmicsDataMap
from omics.id_converter import searchNomenclature

"""
 Function responsible for the integration of different omics data with a metabolic model loaded with framed package.
 see help(integrateOmics)
"""

# Nested function implementation to ease variable scoping


def integrateOmics(CBModel, omicsContainer, modelField='id', and_func=min, or_func=max):
    """
    Function responsible for the integration of different omics data with a metabolic model loaded with framed package.
    Matches model ids for gene_ids, metabolites or reaction ids with those present in the omicsContainer object.

    :param CBModel: (obj) a metabolic model object previously loaded with framed package
    :param omicsContainer: (obj) an omics container object previously created using OmicsContainer class.
    :param modelField: (str) the model field where ids shall be retrieved and used for the integration.(must be either
                        "id" or "name")
    :param and_func:(func) the mathematical function to replace the "AND" operator present in the Gene-Protein-Rules
    :param or_func:(func) the mathematical function to replace the "OR" operator present in the Gene-Protein-Rules


    :return m: (obj) an OmicsDataMap object which contains the mapping between reactions/metabolites and its fluxes
    based on the supplied omics data.
   """

    def suffixAndPrefix():
        """
        Since in metabolic models loaded with "framed" package genes come listed in the shape of G_XXXXX_Y ,
        this function is responsible for converting the data in the omics_container to correctly match the genes
        contained in the metabolic model.
        """

        res = {}
        genes = set()  # set of genes listed in every reaction gpr

        # Check each reaction in the model for its gene protein rule, to register which genes can be useful
        for r_id, reaction in CBModel.reactions.items():
            if reaction.gpr is not None:
                genes.update(reaction.gpr.get_genes())

        psid = {x.split('_')[1]: x for x in genes}  # {9122 : G_9122_2}

        # Converting the ids in the omicsContainer object to the nomenclature used in the metabolic model genes
        model_nomenclature = searchNomenclature(list(psid.keys()))
        if model_nomenclature != omicsContainer.nomenclature:
            omicsContainer.convertIds(model_nomenclature)

        for geneId, exp in omicsContainer.data.items():
            if geneId in psid.keys():
                res[psid[geneId]] = exp
        omicsContainer.set_data(res)

    def g2rIntegrate():
        """
        Handles integration of both proteomics and transcriptomics expression data relying on framed's gene2reaction
        """
        suffixAndPrefix()
        d = gene2reaction(CBModel, omicsContainer.get_Data(), or_func=or_func, and_func=and_func)
        return aux_createMap(d, 'ReactionDataMap')

    def dIntegrate():
        """
        Handles integration of both metabolomics and fluxomics data, which does not require gene to reaction score
        conversion.
        """
        dataTypes = {'fluxomics': 'ReactionsDataMap',
                     'metabolomics': 'MetabolitesDataMap'}
        try:
            target = dataTypes[omicsType.lower()].replace('DataMap', '').lower()  # which entity r we 4 looking in the model?
        except KeyError:
            print('OmicsContainer object has an invalid omicsType')
            return

        skipVal = False  # flag to validate if modelField is correct (i.e. one of Id or Name)

        d = {}  # {entry:val for each entry both in model and omics container, and respective val in container}

        for entry in getattr(CBModel, target).values():  # accessing values from a dic of {'EntityID':EntityObject}'
            if not skipVal:
                try:
                    field = getattr(entry, modelField.lower())
                    skipVal = True
                except AttributeError:
                    print('User defined modelField is not present in the Metabolic Model. Must be "id" or "name"')

            else:
                field = getattr(entry, modelField.lower())  # returns a string with an identifier ('Id' or 'Name')

                if field in omicsContainer.data.keys():
                    d[field] = omicsContainer.data[field]

        return aux_createMap(d, dataTypes[omicsType.lower()])

    def aux_createMap(mMap, mapType):
        m = OmicsDataMap(mMap, mapType)
        return m

    # execution commands

    omicsType = omicsContainer.otype.lower()

    if omicsType.lower() in ['proteomics', 'transcriptomics']:
        return g2rIntegrate()
    else:
        return dIntegrate()








