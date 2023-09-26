import Config.Core as cepgen

class Process(object):
    class Variable:
        linear = 'linear'
        exponential = 'exponential'
        square = 'square'
        powerLaw = 'power_law'

    def __init__(self, outgoing: list[int]):
        print('base object created')
        self._variables = cepgen.Parameters()
        self._outgoing = outgoing

        self.variable = cepgen.Parameters()

    def __eval__(self):
        yield

    def addVariable(self, name, range, type):
        self._variables[name] = cepgen.Parameters(
            range = range,
            type = type
        )
        self.variable[name] = 0.
