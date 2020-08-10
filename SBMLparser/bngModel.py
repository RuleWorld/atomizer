class Parameter:
    def __init__(self, spec):
        # spec is ID, value, constant or not, units
        if len(spec) == 4:
            pid, pval, cts, units = spec
            self.cts = cts
            self.units = units
        elif len(spec) == 2:
            pid, pval = spec
        self.Id = pid
        self.val = pval

    def __str__(self):
        return "{} {} #{}".format(self.Id, self.val, self.units)

    def __repr__(self):
        return str(self)

class Compartment:
    def __init__(self, spec):
        cid, cdim, csize, cmt, unit = spec
        self.Id = cid
        self.dim = cdim
        self.size = csize
        self.cmt = cmt 
        self.unit = unit

    def __str__(self):
        if cmt != '':
            txt = "{} {} {} #{}".format(self.Id, self.dim, self.size)
        else:
            txt = "{} {} {}".format(self.Id, self.dim, self.size)
        return txt

    def __repr__(self):
        return str(self)

class Molecule:
    def __init__(self, spec):
        self.Id = spec['returnID']
        self.initConc = spec['initialConcentration']
        self.initAmount = spec['initialAmount']
        self.isConstant = spec['isConstant']
        self.isBoundary = spec['isBoundary']
        self.compartment = spec['compartment']
        self.name = spec['name']
        self.identifier = spec['identifier']
        self.atomizer_text = spec['atomizer_text']

    def __str__(self):
        return "{}".format(self.atomizer_text)

    def __repr__(self):
        return str(self)

class bngModel:
    '''
    Takes in atomizer stuff and turns everything 
    into objects which can be used to print the 
    final model
    '''
    def __init__(self):
        self.parameters = []
        self.compartments = []
        self.molecules = []
        self.species = []
        self.observables = []
        self.rules = []
        self.metaString = None
        self.noCompartment = None

    def __str__(self):
        txt = ""
        if self.metaString is not None:
            txt += metaString 

        txt += "begin model\n"

        txt += "begin parameters\n"
        for param in self.parameters:
            txt += "  " + str(param) + "\n"
        txt += "end parameters\n"

        if not self.noCompartment:
            txt += "begin compartments\n"
            for comp in self.compartments:
                txt += "  " + str(comp) + "\n"
            txt += "end compartments\n"

        txt += "begin molecule types\n"
        for molec in self.molecules:
            txt += "  " + str(molec) + "\n"
        txt += "end molecule types\n"

        txt += "end model"

        return txt


    def __repr__(self):
        return str((self.parameters, self.compartments, self.molecules))

    def add_parameter(self, param_spec):
        self.parameters.append(Parameter(param_spec))

    def add_compartment(self, comp_spec):
        self.compartments.append(Compartment(comp_spec))

    def add_molecule(self, molec_spec):
        self.molecules.append(Molecule(molec_spec))

    def add_species(self, spec):
        pass

    def add_obs(self, obs):
        pass

