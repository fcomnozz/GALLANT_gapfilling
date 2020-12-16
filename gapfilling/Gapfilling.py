import cobra
import re
from cobra.flux_analysis import gapfill
import gurobipy

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
            print(str(obj) + "not found in the query model")
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

def homology_gapfilling(model, templates, model_obj = None, template_obj = None):
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
    for template in templates:
        template.solver = 'gurobi'
        # result variable will store the reactions ids
        result = gapfill(model, template, demand_reactions = False)
        # log
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
    return model, added_reactions
