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
            txt = "{}".format(self.translator[self.Id])
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

    def __str__(self):
        return "{} = {}".format(self.Id, self.definition)

    def __repr__(self):
        return str(self)

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
        pass

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
        self.species = {}
        self.observables = {}
        self.rules = {}
        self.arules = {}
        self.functions = {}
        self.metaString = ""
        self.noCompartment = None
        self.useID = False
        self.translator = {}

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

        txt += "begin functions\n"
        for func in self.functions.values():
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
        self.rules[rule.Id] = rule

    def make_arule(self):
        return ARule()

    def add_arule(self, arule):
        self.arules[arule.Id] = arule
