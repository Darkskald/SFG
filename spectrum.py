from abc import ABC, abstractmethod
class AbstractSpectrum(ABC):
    """The baseclass which ensures that all necessary properties for plotting
    are available as well as some basic functionality all spectra have in common."""

    @property
    @abstractmethod
    def x(self):
        return self._x
    @x.setter
    def set_x(self, value):
        self._x = value

    @property
    @abstractmethod
    def y(self):
        return self._y
    @y.setter
    def set_y(self, value):
        self._x = value

    @property
    @abstractmethod
    def name(self):
        return self._name

    @x.setter
    def set_name(self, value):
        self._name = value


#base class

# sfg_spectrum

# lt_isotherm

# ir_raman

# uv_vis