import mackinac
from Gapfilling import *
import pandas as pd
import re
from contextlib import contextmanager

# the lines below must be written in the notebook which will allow mackinac to work
"""
mackinac.modelseed.ms_client.url = 'https://p3.theseed.org/services/ProbModelSEED/'
mackinac.workspace.ws_client.url = 'https://p3.theseed.org/services/Workspace'
mackinac.genome.patric_url = 'https://www.patricbrc.org/api/'

# password: ASafSLUqQc@7zyP
mackinac.get_token("fcomnozz")
"""

# Function for avoiding stdout
@contextmanager
def suppress_stdout():
    with open(os.devnull, "w") as devnull:
        old_stdout = sys.stdout
        sys.stdout = devnull
        try:  
            yield
        finally:
            sys.stdout = old_stdout

def get_models(ID_list):
    """
    This function receives a list with PATRIC IDs, reconstructs those models and converts them into COBRA models.
    It returns a dictionary where keys are model IDs and values are the models.
    """
    # Reconstructing models; this step will likely take a long time
    for ID in ID_list:
        # RAST can be used too (source = "rast)
        mackinac.reconstruct_modelseed_model(ID)
    # Creating COBRA models
    d = {}
    for ID in ID_list:
        d[ID] = mackinac.create_cobra_model_from_modelseed_model(model_id=ID)
    return d

def modelseed_gapfilling(model_dictionary):
    """
    Performs ModelSEED gapfilling against SEED reaction database.
    """
    for x in model_dictionary.values():
        with suppress_stdout():
            mackinac.gapfill_modelseed_model(x)

def seed_to_bigg(model, m_table, r_table):
    """
    Translates metabolites and reactions of the model to BiGG nomenclature.
    """
    new_model = cobra.Model(str(model) + " (translated)")
    # metabolites
    for metabolite in model.metabolites:
        m = metabolite.copy()
        sid = m.id
        # assumption -> nx_c; maybe will need to change it
        c = sid[-2:]
        sid = sid[:-2]
        try:
            bid = m_table[m_table["SEED"] == sid]["BiGG"].values[0]
        except IndexError:
            new_model.add_metabolites(m)
            continue
        ID = str(bid)+str(c)
        m.id = ID
        new_model.add_metabolites(m)
    # reactions
    for reaction in model.reactions:
        r = reaction.copy()
        sid = r.id
        # We will use metabolites table for exchange reactions
        # ID
        if sid.startswith('EX_'):
            c = sid[-2:]
            sid = sid[3:-2]
            try:
                bid = m_table[m_table["SEED"] == sid]["BiGG"].values[0]
                bid = "EX_" + str(bid) + str(c)
                r.id = bid
            except IndexError:
                pass
        else:
            c = sid[-2:]
            sid = sid[:-2]
            try:
                bid = r_table[r_table["SEED"] == sid]["BiGG"].values[0]
                # BiGG reactions' IDs don't have compartment termination
                r.id = bid
            except IndexError:
                pass
        # METABOLITES WITHIN REACTION
        # metabolite objects
        metab = [x.copy() for x in r.metabolites]
        # SEED IDs
        met_id = [str(x) for x in r.metabolites.keys()]
        # stoichiometric coefficients
        coef = [float(x) for x in r.metabolites.values()]
        # managing compartment terminations
        met_c = [x[-2:] for x in met_id]
        met_id = [x[:-2] for x in met_id]
        # getting BiGG ids
        new_ids = []
        for x in met_id:
            try:
                new_id = m_table[m_table["SEED"] == x]["BiGG"].values[0]
                new_ids.append(new_id)
            except IndexError:
                new_ids.append(x)
        # dict
        m_dict = {}
        # adding compartment terminations
        for i in range(len(new_ids)):
            new_ids[i] = str(new_ids[i]) + str(met_c[i])
            # changing metabolites' id
            metab[i].id = new_ids[i]
            # making the dict
            m_dict[metab[i]] = coef[i]
        # removing old metabolites
        r.subtract_metabolites(r.metabolites)
        # adding translated metabolites
        r.add_metabolites(m_dict)
        # ADDING REACTION
        new_model.add_reaction(r)
    return new_model

def modelseed(id_list, gapfilling = True, tobigg = True, M_table, R_table):
    """
    Function containing the three above.
    """
    # models is a dict
    models = get_models(ID_list = ID_list)
    if gapfilling == True:
        modelseed_gapfilling(models)
    if tobigg == True:
        bigg_models = {}
        for i in models:
            with suppress_stdout():
                bigg_models[i] = seed_to_bigg(models[i], m_table = M_table, r_table = R_table)
    return bigg_models