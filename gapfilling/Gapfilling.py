import cobra
import re
from cobra.flux_analysis.gapfilling import GapFiller
import gurobipy
from cobra import exceptions

def load_query_model(model, obj = None):
    """
    Loads a model. Objective can be changed to either "biomass" (regex
    will be used) or a valid reaction ID.
    Print is used temporary. A general log will be generated in the 
    final version.
    """
    model = model
    if obj == "biomass":
        b = re.compile("biomass", re.IGNORECASE)
        bc = re.compile("(biomass){1}.*(core){1}", re.IGNORECASE)
        reactions = [reaction.id for reaction in model.reactions]
        # searching for a biomass core reaction
        core = list(filter(bc.match, reactions))
        if core:
            model.objective = core[0]
            return model
        # searching for a non core biomass reaction
        biomass = list(filter(b.match, reactions))
        if biomass:
            model.objective = biomass[0]
            return model
        print("biomass not found in the query model")
        return model
    if obj != "biomass" and obj != None:
        try:
            model.objective = obj
        except ValueError:
            print(str(obj) + " not found in the query model")
            return model
        return model
    if obj == None:
        return model

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
    for reaction in EX:
        if reaction not in model.reactions and str(list(reaction.metabolites.keys())[0]) in [i.id for i in model.metabolites]:
            model.add_reaction(reaction.copy())
    return model

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
    if len(r) != len(p):
        return False
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
    return model

def homology_gapfilling(model, templates, model_obj = None, template_obj = None, use_all_templates = False,
                       integer_threshold = 1e-6, force_exchange = False, force_transport = False, t_all_compounds = False,
                       t_ignore_h = False):
    """
    Performs gap filling on a model using homology models as templates.
    """
    model = load_query_model(model, obj = model_obj)
    model.solver = 'gurobi'
    templates = load_template_models(templates, obj = template_obj)
    # this dict will store used models, genes and reactions
    added_reactions = {}
    # initial flux value
    value = model.optimize().objective_value
    if value == None:
        value = 0.0
    if use_all_templates == False:
        for template in templates:
            # adding exchange reactions
            if force_exchange == True:
                add_exchange_reactions(model, template)
            # adding transport reactions
            if force_transport == True:
                add_transport(model, template, all_compounds = t_all_compounds, ignore_h = t_ignore_h)
            template.solver = 'gurobi'
            try:
                # result variable will store the reactions ids
                result = gapfilling(model, template, integer_threshold = integer_threshold)
                # reactions "log"
                log = []
                for reaction in result[0]:
                    if reaction.id.startswith("EX_"):
                        log.append((reaction.id, "Exchange reaction"))
                    else:
                        log.append((reaction.id, [str(list(reaction.genes)[i]) 
                                              for i in range(len(list(reaction.genes)))]))
                added_reactions[str(template)] = log
                # Adding reactions to the model
                [model.add_reaction(reaction.copy()) for reaction in result[0]]
                # Flux will be evaluated here
                new_value = model.optimize().objective_value
                if new_value != None and new_value > value:
                    value = new_value
                elif new_value == None:
                    continue
                elif new_value != None and new_value == value:
                    break
            except RuntimeError:
                print("\n" + str(template) + ": failed to validate gapfilled model, try lowering the integer_threshold")
            except exceptions.Infeasible:
                print("\n" + str(template) + ": gapfilling optimization failed (infeasible)")
        return model, added_reactions
    else:
        for template in templates:
            # adding exchange reactions
            if force_exchange == True:
                add_exchange_reactions(model, template)
            # adding transport reactions
            if force_transport == True:
                add_transport(model, template, all_compounds = t_all_compounds, ignore_h = t_ignore_h)
            template.solver = 'gurobi'
            try:
                result = gapfilling(model, template, integer_threshold = integer_threshold)
                log = []
                for reaction in result[0]:
                    if reaction.id.startswith("EX_"):
                        log.append((reaction.id, "Exchange reaction"))
                    else:
                        log.append((reaction.id, [str(list(reaction.genes)[i]) 
                                              for i in range(len(list(reaction.genes)))]))
                added_reactions[str(template)] = log
                [model.add_reaction(reaction.copy()) for reaction in result[0]]
            except RuntimeError:
                print("\n" + str(template) + ": failed to validate gapfilled model, try lowering the integer_threshold")
            except exceptions.Infeasible:
                print("\n" + str(template) + ": gapfilling optimization failed (infeasible)") 
        return model, added_reactions
