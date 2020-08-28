from math import sqrt, sin, pi
import warnings
from typing import *


class CircuitError(Exception):
    pass


class Switch:
    def __init__(self):
        self.__closed = True

    @property
    def is_closed(self) -> bool:
        return self.__closed

    @property
    def is_open(self) -> bool:
        return not self.is_closed

    def open(self):
        self.__closed = False

    def close(self):
        self.__closed = True

    def toggle(self):
        self.__closed = not self.__closed

    def __repr__(self) -> str:
        return f"{type(self)} ({id(self)}): is_closed = {self.is_closed}"


number = Union[int, float, Callable[[float], Union[float, int]]]


class Battery:
    def __init__(self, emf: number):
        if type(emf) in (int, float):
            self.__emf = lambda t: emf
        elif hasattr(emf, '__call__'):
            self.__emf = emf

    def emf(self, t: float) -> float:
        return self.__emf(t)

    def set(self, emf: number):
        if type(emf) in (int, float):
            self.__emf = lambda t: emf
        elif hasattr(emf, '__call__'):
            self.__emf = emf

    def __repr__(self) -> str:
        return f"{type(self)} ({id(self)})"


class ACBattery(Battery):
    def __init__(self, emf_rms: float, omega: float, phi: float):
        if omega == 0.0:
            raise CircuitError("omega cannot be zero for ac battery")

        super().__init__(lambda t: sqrt(2) * emf_rms * sin(omega * t + phi))
        self.__emf_rms = emf_rms
        self.__omega = omega
        self.__phi = phi

    @property
    def emf_rms(self) -> float:
        return self.__emf_rms

    @property
    def omega(self) -> float:
        return self.__omega

    @property
    def time_period(self) -> float:
        return 2 * pi / self.__omega

    @property
    def frequency(self) -> float:
        return 1 / self.time_period

    @property
    def phi(self) -> float:
        return self.__phi

    def __repr__(self) -> str:
        return f"{type(self)} ({id(self)}): emf_rms = {self.__emf_rms}, " \
               f"omega = {self.__omega}, " \
               f"phi = {self.__phi}"


class Resistor:
    def __init__(self, resistance: number):
        if type(resistance) in (int, float):
            self.__resistance = lambda t: resistance
        elif hasattr(resistance, '__call__'):
            self.__resistance = resistance

    def resistance(self, t: float) -> float:
        res = self.__resistance(t)
        if res == 0.0:
            pass  # warnings.warn(f"0 resistance was encountered with {self}", RuntimeWarning)
        return res

    def potential_drop(self, t: float, current: float) -> float:
        return self.__resistance(t) * current

    def __repr__(self) -> str:
        return f"{type(self)} ({id(self)})"


class Capacitor:
    def __init__(self, capacitance: number, charge=0.0):
        if type(capacitance) in (int, float):
            self.__capacitance = lambda t: capacitance
        elif hasattr(capacitance, "__call__"):
            self.__capacitance = capacitance
        self.__init_charge = charge

    def capacitance(self, t: float) -> float:
        res = self.__capacitance(t)
        if res == 0.0:
            warnings.warn("0 capacitance was encountered", RuntimeWarning)
        return res

    @property
    def init_charge(self) -> float:
        return self.__init_charge

    def potential_drop(self, t: float, charge: float) -> float:
        return charge / self.__capacitance(t) if self.__capacitance(t) != 0.0 else 0.0

    def __repr__(self) -> str:
        return f"{type(self)} ({id(self)})"


class Inductor:
    def __init__(self, inductance: number):
        self.__isnull = False
        if type(inductance) in (int, float):
            if inductance == 0.0:
                self.__isnull = True
            else:
                self.__inductance = lambda t: inductance
        elif hasattr(inductance, "__call__"):
            self.__inductance = inductance

    def is_not_null(self) -> bool:
        return not self.__isnull

    def inductance(self, t: float) -> float:
        res = self.__inductance(t)
        if res == 0.0:
            warnings.warn("0 inductance was encountered", RuntimeWarning)
        return res

    def potential_drop(self, t: float, di_dt: float) -> float:
        return self.__inductance(t) * di_dt

    def __repr__(self) -> str:
        return f"{type(self)} ({id(self)}): inductance = {self.inductance}"


class Galvanometer(Resistor):
    def __init__(self, resistance: number, max_current: float):
        super().__init__(resistance)
        if max_current < 0.0:
            raise CircuitError(f"{type(self)}.max_current must be non-negative, cannot be: {max_current}")
        self.__max_current = max_current

    @property
    def max_current(self) -> float:
        return self.__max_current

    def get_reading(self, current: float) -> float:
        return current if abs(current) < self.__max_current else self.__max_current

    def __repr__(self) -> str:
        return super().__repr__() + f", max_current = {self.__max_current}"


Ammeter = Galvanometer


class Voltmeter(Resistor):
    def __init__(self, resistance: number, max_voltage: float):
        super().__init__(resistance)
        if max_voltage < 0.0:
            raise CircuitError(f"{type(self)}.max_voltage must be non-negative, cannot be: {max_voltage}")
        self.__max_voltage = max_voltage

    @property
    def max_voltage(self) -> float:
        return self.__max_voltage

    def get_reading(self, t: float, current: float) -> float:
        voltage = self.potential_drop(t, current)
        return voltage if abs(voltage) < self.max_voltage else self.max_voltage

    def __repr__(self) -> str:
        return super().__repr__() + f", max_voltage = {self.max_voltage}"
