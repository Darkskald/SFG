import json

import numpy as np
import pandas as pd
from scipy.integrate import simps as sp
from scipy.integrate import trapz as tp

from SFG.spectrum.exceptions import InvalidSpectrumError, IntegrationError


class MetaSpectrum(type):

    def __call__(cls, *args, **kwargs):
        temp = super().__call__(*args, **kwargs)

        for attr in ("x", "y", "x_unit", "y_unit", "name"):
            if not hasattr(temp, attr):
                raise InvalidSpectrumError(f'Tried to instantiate spectrum object without suitable attribute {attr}!')

        return temp


class BaseSpectrum(metaclass=MetaSpectrum):

    def __init__(self, name=None, x=None, y=None, x_unit=None, y_unit=None, timestamp=None):
        self.name = name
        self.x = x
        self.y = y
        self.x_unit = x_unit
        self.y_unit = y_unit
        self.timestamp = timestamp

    # magic methods

    def __repr__(self):
        return f'{type(self).__name__} Object with name "{self.name}"'

    # working on the data
    def yield_spectral_range(self):
        """returns a list containing maximum and minimum wavenumer and the number of data points"""
        return [min(self.x), max(self.x), len(self.x)]

    def get_xrange_indices(self, lower, upper):
        """Takes a high (upper) and a low (lower) target x value as argument. Returns
        the indices of the wavenumber array of the spectrum that are the borders of this interval."""
        lower_index = np.argmax(self.x >= lower)
        upper_index = np.argmax(self.x >= upper)
        return int(lower_index), int(upper_index)

    def get_xrange(self):
        # todo: ensure this functions work as well for y_values
        """Returns the slice of the x values in the borders of lower to upper"""
        lower, upper = self.get_xrange_indices()
        return self.x[lower, upper + 1]

    def normalize(self, external=None):
        """Normalize the spectrum's y data either to the maximum of the y values or an
        external factor"""
        if external is None:
            return np.max(self.y)
        else:
            return self.y/external

    def integrate_slice(self, x_array, y_array):
        """Integrates the y_array which has a spacing given by the x_array. First it tries to apply
        simpson rule integration rule, but if this fails the function invokes integration via
        trapeziodal rule"""
        try:
            area = sp(y_array, x_array)
            if np.isnan(area):
                print(f"""Integration failed in spectrum {self.name} using Simpson's rule. 
                Falling back to trapezoidal rule.""")
                area = tp(y_array, x_array)
            return area
        except:
            raise IntegrationError(f'Integration not possible for {self.name}')

    # export functions
    def properties_to_dict(self):
        temp = {
                "name": self.name,
                "x_unit": self.x_unit,
                "y_unit": self.y_unit,
                "x": self.x,
                "y": self.y,
                "timestamp": self.timestamp
                }
        return temp

    def to_pandas_dataframe(self) -> pd.DataFrame:
       pd.DataFrame(data=self.properties_to_dict())

    def to_csv(self) -> None:
        self.to_pandas_dataframe().to_csv(self.name + ".csv", index=False, sep=";")

    def to_json(self) -> str:
        temp = self.properties_to_dict()
        return json.dumps(temp)




