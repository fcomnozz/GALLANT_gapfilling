import cobra
import re
from cobra.flux_analysis.gapfilling import GapFiller
from cobra.io import read_sbml_model, write_sbml_model
import gurobipy
from cobra import exceptions
import os, sys

def fake_copy(model):
    write_sbml_model(model, filename = "tmp.xml")
    copied_model = read_sbml_model("tmp.xml")
    os.remove("tmp.xml")
    return copied_model

def load_query_model(model, obj = None):
    """
    Loads a model. Objective can be changed to either "biomass" (regex
    will be used) or a valid reaction ID.
    Print is used temporary. A general log will be generated in the 
    final version.
    """
    Model = model
    if obj == "biomass":
        b = re.compile("biomass", re.IGNORECASE)
        bc = re.compile("(biomass){1}.*(core){1}", re.IGNORECASE)
        reactions = [reaction.id for reaction in Model.reactions]
        # searching for a biomass core reaction
        core = list(filter(bc.match, reactions))
        if core:
            Model.objective = core[0]
            return Model
        # searching for a non core biomass reaction
        biomass = list(filter(b.match, reactions))
        if biomass:
            Model.objective = biomass[0]
            return Model
        print("biomass not found in the query model")
        return Model
    if obj != "biomass" and obj != None:
        try:
            Model.objective = obj
        except ValueError:
            print(str(obj) + " not found in the query model")
            return Model
        return Model
    if obj == None:
        return Model

def load_template_models(template_list, obj = None):
    """
    Takes a list of template models and changes objective if specified.
    Objective can be either "biomass" or a specific reaction ID.
    """
    templates = template_list
    failures = []
    if obj == None:
        return templates
    if obj == "biomass":
        b = re.compile("biomass", re.IGNORECASE)
        bc = re.compile("(biomass){1}.*(core){1}", re.IGNORECASE) 
        for model in templates:
            reactions = [reaction.id for reaction in model.reactions]
            # Searching for a biomass_core reaction
            core = list(filter(bc.match, reactions))
            if core:
                model.objective = core[0]
                continue
            # Searching for a non core biomass reaction
            biomass = list(filter(b.match, reactions))
            if biomass:
                model.objective = biomass[0]
            # If biomass reactions are not found, the model name is stored into "failures" list.
            # The model won't change its objective, but it will be used anyway for gap filling
            else:
                failures.append(str(model))
            print(failures)
            return templates
    if obj != None and obj != "biomass":
        for model in templates:
            try:
                model.objective = obj
            except ValueError:
                failures.append(str(model))
        print(failures)
        return templates
    
def gapfilling(model, template, integer_threshold, it = 1):
    """
    Calls GapFiller class from cobrapy.
    """
    gapfiller = GapFiller(model, template, integer_threshold = integer_threshold,
                         demand_reactions = False)
    return gapfiller.fill(iterations = it)

def add_exchange_reactions(model, template):
    """
    Adds all exchange reactions from a template matching model metabolites. 
    """
    EX = [template.reactions.get_by_id(i.id) for i in template.reactions if i.id.startswith("EX_")]
    added_exchange = []
    for reaction in EX:
        if reaction not in model.reactions and str(list(reaction.metabolites.keys())[0]) in [i.id for i in model.metabolites]:
            model.add_reaction(reaction.copy())
            added_exchange.append(reaction.id)
    return model, added_exchange

def is_transport(reaction, compartments_list, all_compounds = False, ignore_h = False):
    """
    Takes a cobra model reaction as input and determines if it is a transport reaction.
    all_compounds = True -> both sides of the reaction must be the same
    all_compounds = False -> at least one compound must be in both sides
    ignore_h = True -> does not consider proton transference as transport
    """
    # left part of the reaction
    r = list(str(x) for x in reaction.reactants)
    # right part
    p = list(str(x) for x in reaction.products)
    #if len(r) != len(p):
    #   return False
    c = compartments_list
    # removing terminations
    R = []
    for i in r:
        for e in c:
            i = i.replace(e, "") 
        R.append(i)
    P = []
    for i in p:
        for e in c:
            i = i.replace(e, "")
        P.append(i)
    # we sort the lists to avoid missing transport reactions where the sequence of the compounds is not maintained
    R = sorted(R)
    P = sorted(P)
    if all_compounds == False:
        if ignore_h == True:
            if "h" in R:
                R.remove("h")
            if "h" in P:
                P.remove("h")
        for i in R:
            if i in P:
                return True
            else:
                return False
    if R == P:
        return True
    else:
        return False
    
def add_transport(model, template, all_compounds = False, ignore_h = False):
    """
    Adds transport reactions from a template which metabolites are present in the model.
    """
    # PREPARATION
    # added transport reactions register
    added_transport = []
    # template compartments
    t_compartments = list(template.compartments.keys())
    t_Compartments = []
    for i in t_compartments:
        t_Compartments.append("_" + str(i))
    # getting compartments and metabolites from query model for further use
    m_compartments = list(model.compartments.keys())
    m_Compartments = []
    for i in m_compartments:
        m_Compartments.append("_" + str(i))
    m_metabolites = list(str(x) for x in model.metabolites)
    # removing suffixes from metabolites
    m_Metabolites = []
    for i in m_metabolites:
        for e in m_Compartments:
            i = i.replace(e, "")
        m_Metabolites.append(i)
    # sorting and removing duplicates
    m_Metabolites = list(dict.fromkeys(m_Metabolites))
    # REACTION ADDING
    for reaction in template.reactions:
        if is_transport(reaction, t_Compartments, all_compounds = all_compounds, ignore_h = ignore_h) and reaction not in model.reactions:
            # we will only use reactants as they'll be the same as products (apart from location)
            m = list(str(x) for x in reaction.reactants)
            M = []
            for i in m:
                for e in t_Compartments:
                    i = i.replace(e, "")
                M.append(i)
            if all(x in m_Metabolites for x in M):
                model.add_reaction(reaction.copy())
                added_transport.append(('Transport reaction: ' + str(reaction.id), [str(list(reaction.genes)[i]) for i in range(len(list(reaction.genes)))]))
    return model, added_transport

def homology_gapfilling(model, templates, model_obj = None, template_obj = None, use_all_templates = False,
                       integer_threshold = 1e-6, force_exchange = False, force_transport = False, t_all_compounds = False,
                       t_ignore_h = False, value_fraction = 0.8):
    """
    Performs gap filling on a model using homology models as templates.
    """
    Model = load_query_model(model, obj = model_obj)
    Model.solver = 'gurobi'
    Templates = load_template_models(template_list = templates, obj = template_obj)
    # this dict will store used models, genes and reactions
    added_reactions = {}
    # initial flux value
    value = Model.optimize().objective_value
    if value == None:
        value = 0.0
    if use_all_templates == False:
        for template in Templates:
            # reactions "log"
            log = []
            # adding exchange reactions
            if force_exchange == True:
                Model, added_exchange = add_exchange_reactions(Model, template)
                for reaction in added_exchange:
                    log.append((reaction, "Exchange reaction"))
            # adding transport reactions
            if force_transport == True:  
                Model, added_transport = add_transport(Model, template, all_compounds = t_all_compounds, ignore_h = t_ignore_h)
                for reaction in added_transport:
                    log.append(reaction)
            template.solver = 'gurobi'
            try:
                # result variable will store the reactions ids
                result = gapfilling(Model, template, integer_threshold = integer_threshold)
                for reaction in result[0]:
                    if reaction.id.startswith("EX_"):
                        log.append((reaction.id, "Exchange reaction"))
                    else:
                        log.append((reaction.id, [str(list(reaction.genes)[i]) 
                                              for i in range(len(list(reaction.genes)))]))
                added_reactions[str(template)] = log
                # Adding reactions to the model
                [Model.add_reaction(reaction.copy()) for reaction in result[0]]
                # Flux will be evaluated here
                new_value = Model.optimize().objective_value
                if new_value != None and new_value > value:
                    value = new_value
                elif new_value == None:
                    continue
                elif new_value != None and new_value == value:
                    break
                elif new_value < value:
                    if new_value >= value * value_fraction:
                        for i in range(len(log)):
                            Model.reactions.get_by_id(log[i][0]).lower_bound = 0.
                            Model.reactions.get_by_id(log[i][0]).upper_bound = 0.
                        new_value = Model.optimize().objective_value
                        value = new_value
                    else:
                        [Model.remove_reactions(Model.reactions.log[i][0]) for i in range(len(log))]
                        del added_reactions[str(template)]
                        break
                write_sbml_model(Model, filename = "model.tmp")
                Model = read_sbml_model("model.tmp")
                os.remove("model.tmp")
            except RuntimeError:
                print("\n" + str(template) + ": failed to validate gapfilled model, try lowering the integer_threshold")
            except exceptions.Infeasible:
                print("\n" + str(template) + ": gapfilling optimization failed (infeasible)")
        return Model, added_reactions
    else:
        for template in Templates:
            log = []
            # adding exchange reactions
            if force_exchange == True:
                Model, added_exchange = add_exchange_reactions(Model, template)
                for reaction in added_exchange:
                    log.append((reaction, "Exchange reaction"))
            # adding transport reactions
            if force_transport == True:
                Model, added_transport = add_transport(Model, template, all_compounds = t_all_compounds, ignore_h = t_ignore_h)
                for reaction in added_transport:
                    log.append(reaction)
            template.solver = 'gurobi'
            try:
                result = gapfilling(Model, template, integer_threshold = integer_threshold)   
                for reaction in result[0]:
                    if reaction.id.startswith("EX_"):
                        log.append((reaction.id, "Exchange reaction"))
                    else:
                        log.append((reaction.id, [str(list(reaction.genes)[i]) 
                                              for i in range(len(list(reaction.genes)))]))
                added_reactions[str(template)] = log
                [Model.add_reaction(reaction.copy()) for reaction in result[0]]
                write_sbml_model(Model, filename = "model.tmp")
                Model = read_sbml_model("model.tmp")
                os.remove("model.tmp")
            except RuntimeError:
                print("\n" + str(template) + ": failed to validate gapfilled model, try lowering the integer_threshold")
            except exceptions.Infeasible:
                print("\n" + str(template) + ": gapfilling optimization failed (infeasible)") 
        return Model, added_reactions
