import re 

class Parameter:
    def __init__(self):
        # spec is ID, value, constant or not, units
        self.Id = None
        self.val = None
        self.cts = None
        self.units = None
        self.rxn_ind = None

    def __str__(self):
        if self.units is not None and self.units != "":
            return "{} {} #{}".format(self.Id, self.val, self.units)
        else:
            return "{} {}".format(self.Id, self.val)

    def __repr__(self):
        return str(self)

class Compartment:
    def __init__(self):
        self.Id = None
        self.dim = None
        self.size = None
        self.cmt = None 
        self.unit = None

    def __str__(self):
        if cmt != '':
            txt = "{} {} {} #{}".format(self.Id, self.dim, self.size)
        else:
            txt = "{} {} {}".format(self.Id, self.dim, self.size)
        return txt

    def __repr__(self):
        return str(self)

class Molecule:
    def __init__(self):
        self.translator = {}

    def parse_raw(self, raw):
        self.raw = raw
        self.Id = raw['returnID']
        self.initConc = raw['initialConcentration']
        self.initAmount = raw['initialAmount']
        self.isConstant = raw['isConstant']
        self.isBoundary = raw['isBoundary']
        self.compartment = raw['compartment']
        self.name = raw['name']
        self.identifier = raw['identifier']

    def __str__(self):
        if self.Id in self.translator:
            # str2 is molecule types?
            txt = "{}".format(self.translator[self.Id].str2())
        else:
            txt = "{}()".format(self.Id)
        return txt

    def __repr__(self):
        return str(self)

class Species:
    def __init__(self):
        self.noCompartment = False
        self.translator = {}

    def parse_raw(self, raw):
        self.raw = raw
        self.Id = raw['returnID']
        self.initConc = raw['initialConcentration']
        self.initAmount = raw['initialAmount']
        self.isConstant = raw['isConstant']
        self.isBoundary = raw['isBoundary']
        self.compartment = raw['compartment']
        self.name = raw['name']
        self.identifier = raw['identifier']
        if self.initAmount > 0:
            self.val = self.initAmount
        elif self.initConc > 0:
            # TODO: Figure out what to do w/ conc
            self.val = self.initConc
            
    def __str__(self):
        trans_id = self.translator[self.Id] if self.Id in self.translator else self.Id+"()"
        mod = "$" if self.isConstant else ""
        if self.noCompartment or self.compartment == "":
            txt = "{}{} {} #{} #{}".format(mod, trans_id, self.val, self.raw['returnID'], self.raw['identifier'])
        else:
            txt = "@{}:{}{} {} #{} #{}".format(self.compartment, mod, trans_id, self.val, self.raw['returnID'], self.raw['identifier'])
        return txt

    def __repr__(self):
        return str(self)

class Observable:
    def __init__(self):
        self.Id = None
        self.type = "Species"
        self.compartment = None
        self.noCompartment = False
        self.translator = {}
        
    def parse_raw(self, raw):
        self.raw = raw
        self.Id = raw['returnID']
        self.initConc = raw['initialConcentration']
        self.initAmount = raw['initialAmount']
        self.isConstant = raw['isConstant']
        self.isBoundary = raw['isBoundary']
        self.compartment = raw['compartment']
        self.name = raw['name']
        self.identifier = raw['identifier']

    def __str__(self):
        txt = self.type
        pattern = self.translator[self.raw['returnID']] if self.Id in self.translator else self.raw['returnID']+"()"
        if self.noCompartment or self.compartment == "":
            txt += " {0} {1} #{2}".format(self.Id, pattern, self.name)
        else:
            txt += " {0}_{2} @{2}:{1} #{3}".format(self.Id, pattern, self.compartment, self.name)
        return txt

    def __repr__(self):
        return str(self)

class Function:
    def __init__(self):
        self.Id= None
        self.definition = None
        self.rule_ptr = None
        self.local_dict = None
        self.replaceLocParams = False

    def replaceLoc(self, func_def, pdict):
        for parameter in pdict:
            func_def = re.sub(r'(\W|^)({0})(\W|$)'.format(parameter),r'\g<1>{0}\g<3>'.format(pdict[parameter]),func_def)
        return func_def

    def renameLoc(self, pname, rind):
        return "r{}_{}".format(rind, pname)

    def __str__(self):
        # TODO: implement parameter replacement here
        # now we have a pointer to the rule itself, we can 
        # use that to resolve the parameter values, we can 
        # link to the overall model to get the values there
        # or we can just re-generate the local parameter names
        # and use those
        fdef = self.definition
        if self.replaceLocParams:
            # check possible places, local dict first
            if self.local_dict is not None:
                fdef = self.replaceLoc(self.definition, self.local_dict)
            # or pull from the pointer to the rule itself
            elif self.rule_ptr is not None:
                if len(self.rule_ptr.raw_param) > 0:
                    rule_dict = dict([(i[0],i[1]) for i in self.rule_ptr.raw_param])
                    fdef = self.replaceLoc(self.definition, rule_dict)
        # if we are not replacing, we need to rename local parameters
        # to the correct index if the function is related to a rule 
        else:
            if self.rule_ptr is not None:
                # this is a fRate, check for local parameters
                if len(self.rule_ptr.raw_param) > 0:
                    # gotta rename these in the function
                    rule_dict = dict([(i[0],self.renameLoc(i[0],self.rule_ptr.rule_ind)) for i in self.rule_ptr.raw_param])
                    fdef = self.replaceLoc(self.definition, rule_dict)
        fdef = self.adjust_func_def(fdef)
        return "{} = {}".format(self.Id, fdef)

    def __repr__(self):
        return str(self)

    def adjust_func_def(self,fdef):
        # TODO: All function definitions has to be adjusted to 
        # work with BNGL. E.g. switching log to ln and many more
        # adjustments 
        return fdef

class Rule:
    def __init__(self):
        self.Id= ""
        self.reactants = []
        self.products = []
        self.rate_cts = (None,)
        self.comment = ""
        self.reversible = False
        self.translator = {}

    def parse_raw(self, raw):
        self.raw = raw
        self.raw_react = raw['reactants']
        self.raw_prod = raw['products']
        self.raw_param = raw['parameters']
        self.raw_rates = raw['rates']
        self.raw_num = raw['numbers']
        self.raw_splt = raw['split_rxn']
        self.reversible = raw['reversible']
        self.Id = raw['reactionID']

    def __str__(self):
        if self.Id != "":
            txt = "{}: ".format(self.Id)
        else:
            txt = ""
        # reactants 
        if len(self.reactants) == 0:
            txt += "0"
        else:
            for ir,react in enumerate(self.reactants):
                if ir != 0:
                    txt += " + "
                txt += str(self.translator[react[0]]) if react[0] in self.translator else str(react[0])+"()"
        # correct rxn arrow
        if self.reversible and len(self.rate_cts) == 2:
            txt += " <-> "
        else:
            txt += " -> "
        # products
        if len(self.products) == 0:
            txt += "0"
        else:
            for ip,prod in enumerate(self.products):
                if ip != 0:
                    txt += " + "
                txt += str(self.translator[prod[0]]) if prod[0] in self.translator else str(prod[0])+"()"
        # rate cts
        if self.reversible and len(self.rate_cts) == 2:
            txt += " {},{}".format(self.rate_cts[0], self.rate_cts[1])
        else:
            txt += " {}".format(self.rate_cts[0])
        
        comment = 'Modifiers({0})'.format(', '.join(self.raw['modifiers'])) if self.raw['modifiers'] else ''
        if comment != "":
            txt += " #{}".format(comment)
        return txt

    def __repr__(self):
        return str(self)

class ARule:
    def __init__(self):
        self.Id = None
        self.rates = None
        self.isAssignment = None
        self.isRate = None

    def parse_raw(self, raw):
        self.Id = raw[0]
        self.rates = raw[1]
        self.isAssignment = raw[2]
        self.isRate = raw[3]

    def __str__(self):
        return "{} {}".format(self.Id, self.rates)

    def __repr__(self):
        return str(self)

# Model objects done 

class bngModel:
    '''
    Takes in atomizer stuff and turns everything 
    into objects which can be used to print the 
    final model
    '''
    def __init__(self):
        self.parameters = {}
        self.compartments = {}
        self.molecules = {}
        self.molecule_ids = {}
        self.species = {}
        self.observables = {}
        self.rules = {}
        self.arules = {}
        self.functions = {}
        self.metaString = ""
        self.translator = {}
        self.molecule_mod_dict = {}
        self.noCompartment = None
        self.useID = False
        self.replaceLocParams = False

    def __str__(self):
        txt = self.metaString

        txt += "begin model\n"

        txt += "begin parameters\n"
        for param in self.parameters.values():
            txt += "  " + str(param) + "\n"
        txt += "end parameters\n"

        if not self.noCompartment:
            txt += "begin compartments\n"
            for comp in self.compartments.values():
                txt += "  " + str(comp) + "\n"
            txt += "end compartments\n"

        txt += "begin molecule types\n"
        for molec in self.molecules.values():
            molec.translator = self.translator
            txt += "  " + str(molec) + "\n"
        txt += "end molecule types\n"

        txt += "begin seed species\n"
        for spec in self.species.values():
            spec.translator = self.translator
            if spec.val > 0:
                spec.noCompartment = self.noCompartment
                txt += "  " + str(spec) + "\n"
        txt += "end seed species\n"

        txt += "begin observables\n"
        for obs in self.observables.values():
            obs.translator = self.translator
            obs.noCompartment = self.noCompartment
            txt += "  " + str(obs) + "\n"
        txt += "end observables\n"

        if len(self.functions) > 0:
            txt += "begin functions\n"
            for func in self.functions.values():
                func.replaceLocParams = self.replaceLocParams
                txt += "  " + str(func) + "\n"
            txt += "end functions\n"

        txt += "begin reaction rules\n"
        for rule in self.rules.values():
            rule.translator = self.translator
            txt += "  " + str(rule) + "\n"
        txt += "end reaction rules\n"

        txt += "end model"

        return txt

    def __repr__(self):
        return str((self.parameters, self.molecules))

    def consolidate_arules(self):
        '''
        this figures out what to do with particular
        assignment rules pulled from SBML. 
        a) A non-constant parameter can be turned
           into an ODE basically
        b) Any species in the system can be modified
           by an assignment rule. This turns the species
           into a function which also requires a modification
           of any reaction rules the species is associated with
        '''
        for arule in self.arules.values():
            # first one is to check parameters
            # import IPython;IPython.embed()
            if arule.isRate:
                # this is a rate rule, it'll be turned into a 
                # reaction/function. 
                pass
            elif arule.isAssignment:
                # rule is an assignment rule
                # let's first check parameters
                if arule.Id in self.parameters:
                    a_param = self.parameters[arule.Id]
                    if not a_param.cts:
                        # this means that one of our parameters 
                        # is _not_ a constant and is modified by 
                        # an assignment rule
                        # TODO: Not sure if anything else 
                        # can happen here. Confirm via SBML spec
                        a_param = self.parameters.pop(arule.Id)
                        # TODO: check if an initial value to 
                        # a non-constant parameter is relevant?
                        # I think the only thing we need is to 
                        # turn this into a function
                        fobj = self.make_function()
                        fobj.Id = arule.Id
                        fobj.definition = arule.rates[0]
                        self.add_function(fobj)
                elif arule.Id in self.molecule_ids:
                    # we are an assignment rule that modifies 
                    # a molecule, this will be converted to 
                    # a function if true
                    mname = self.molecule_ids[arule.Id]
                    molec = self.molecules[mname]
                    # We can't have the molecule be _constant_
                    # at which point it's supposed to be encoded
                    # with "$" in BNGL
                    if not molec.isConstant:
                        # we can have it be boundary or not, doesn't 
                        # matter since we know an assignment rule is
                        # modifying it and it will take over reactions

                        # this should be guarantee
                        molec = self.molecules.pop(mname)

                        # we should also remove this from species
                        # and/or observables, this checks for 
                        # namespace collisions. 
                        # TODO: We might want to 
                        # remove parameters as well 
                        if molec.name in self.observables:
                            obs = self.observables.pop(molec.name)
                        if molec.name in self.species:
                            spec = self.species.pop(molec.name)

                        # this will be a function
                        fobj = self.make_function()
                        fobj.Id = mname + "()"
                        fobj.definition = arule.rates[0]
                        self.add_function(fobj)

                        # TODO: This molecule should be
                        # converted into a funcion and a function only
                        # and the reactions it's a part of should be 
                        # modified appropriately 
                        if mname in self.molecule_mod_dict:
                            for rule in self.molecule_mod_dict[mname]:
                                print("{} is being turned to a function, change rule: {}".format(mname, rule))
                else:
                    # this is just a simple assignment (hopefully)
                    # just convert to a function
                    fobj = self.make_function()
                    fobj.Id = arule.Id + "()"
                    fobj.definition = arule.rates[0]
                    self.add_function(fobj)
            else:
                # not sure what this means, read SBML spec more
                pass

    def consolidate_molecules(self):
        to_remove = []
        for molec in self.molecules:
            if molec not in self.molecule_mod_dict:
                to_remove.append(molec)
        for molec in to_remove:
            m = self.molecules.pop(molec)

    def consolidate(self):
        self.consolidate_arules()
        # self.consolidate_molecules()

    def reorder_functions(self):
        '''
        this one is to make sure the functions are reordered
        correctly, should be ported from the original codebase
        '''
        pass

    # model object creator and adder methods 
    def add_parameter(self, param):
        self.parameters[param.Id] = param

    def make_parameter(self):
        return Parameter()

    def add_compartment(self, comp):
        self.compartments[comp.Id] = comp

    def make_compartment(self):
        return Compartment()

    def add_molecule(self, molec):
        # we might have to add molecules that
        # didn't have rawSpecies associated with 
        if hasattr(molec, "raw"):
            self.molecule_ids[molec.raw['identifier']] = molec.name
        self.molecules[molec.name] = molec

    def make_molecule(self):
        return Molecule()

    def add_species(self, sspec):
        self.species[sspec.name] = sspec

    def make_species(self):
        return Species()

    def add_observable(self, obs):
        self.observables[obs.Id] = obs

    def make_observable(self):
        return Observable()

    def make_function(self):
        return Function()

    def add_function(self, func):
        self.functions[func.Id] = func

    def make_rule(self):
        return Rule()

    def add_rule(self, rule):
        # add this to keep track of molecule modifications
        # this will allow us to quickly loop over reactons 
        # a molecule is a part of if the molecule gets 
        # turned to a function
        for react in rule.reactants:
            if react[0] not in self.molecule_mod_dict:
                self.molecule_mod_dict[react[0]] = []
            self.molecule_mod_dict[react[0]].append(rule)
        for prod in rule.products:
            if prod[0] not in self.molecule_mod_dict:
                self.molecule_mod_dict[prod[0]] = []
            self.molecule_mod_dict[prod[0]].append(rule)
        self.rules[rule.Id] = rule

    def make_arule(self):
        return ARule()

    def add_arule(self, arule):
        self.arules[arule.Id] = arule
