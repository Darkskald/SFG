class InvalidSpectrumError(Exception):
    def __init__(self, *args, **kwargs):
        Exception.__init__(self, *args, **kwargs)


class IntegrationError(Exception):
    def __init__(self,*args,**kwargs):
        Exception.__init__(self, *args, **kwargs)


class CoverageCalculationImpossibleError(Exception):
    pass