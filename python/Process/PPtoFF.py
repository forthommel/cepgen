import Config.Core as cepgen


class PPtoFF(cepgen.Process):
    def __init__(self):
        super().__init__()
        print('PPtoFF process initialised')
        self.addVariable('q2', range=(1.e-5, 10.), type=cepgen.Process.Variable.exponential)

        self.mf2 = 0.10566**2

    def __repr__(self):
        return 'PPtoFF process'

    def __eval__(self):
        q2 = self.variable['q2']
        return 1./q2
