import Config.Core as cepgen


class PPtoFF(cepgen.Process):
    def __init__(self):
        super().__init__([13, 13])
        print('PPtoFF process initialised')
        self.addVariable('q2', range=(1.e-5, 10.), type=cepgen.Process.Variable.exponential)

    def __repr__(self):
        return 'PPtoFF process'

    def __eval__(self):
        q2 = self.variable['q2']
        return 1./q2
