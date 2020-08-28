from __future__ import annotations
from .devices import *


class Wire:
    def __init__(self, resistor: number, battery: number = 0.0, ac_battery=(0.0, 1.0, 0.0),
                 capacitor: Tuple[number, float] = (0.0, 0.0), inductor: number = 0.0,
                 galvanometer: Tuple[number, float] = (0.0, 0.0),
                 ammeter: Tuple[number, float] = (0.0, 0.0),
                 voltmeter: Tuple[number, float] = (0.0, 0.0)):

        self.switch = Switch()

        self.resistor = Resistor(resistor)
        self.battery = Battery(battery)
        self.ac_battery = ACBattery(*ac_battery)

        self.capacitor = Capacitor(*capacitor)
        self.inductor = Inductor(inductor)

        self.galvanometer = Galvanometer(*galvanometer)
        self.ammeter = Ammeter(*ammeter)
        self.voltmeter = Voltmeter(*voltmeter)

        self.starting_junction = Junction()
        self.ending_junction = Junction()

        self.q = []
        self.dq_dt = []
        self.d2q_dt2 = []

        self.potential_drop = []
        self.eq_resistance = []
        self.eq_capacitance = []
        self.eq_inductance = []

        self.power = []
        self.heat = []

        self.inductor_energy = []
        self.capacitor_energy = []

        self.galvanometer_reading = []
        self.ammeter_reading = []
        self.voltmeter_reading = []

    def net_resistance(self, t: float) -> float:
        return self.resistor.resistance(t) \
               + self.galvanometer.resistance(t) \
               + self.ammeter.resistance(t) \
               + self.voltmeter.resistance(t)

    def connect(self, junction: Junction):
        if self.starting_junction.is_null:
            self.starting_junction = junction
        elif self.ending_junction.is_null:
            self.ending_junction = junction
        else:
            raise CircuitError(f"Cannot 3rd connect wire to Wire: {self}")

    @property
    def is_null(self) -> bool:
        return self.switch.is_open or self.starting_junction.is_null or self.ending_junction.is_null\
               or self.starting_junction.is_singular or self.ending_junction.is_singular

    def swap_junctions(self):
        self.starting_junction, self.ending_junction = self.ending_junction, self.starting_junction

    def other_junction(self, junction: Junction):
        if junction is self.starting_junction:
            return self.ending_junction
        elif junction is self.ending_junction:
            return self.starting_junction
        else:
            raise CircuitError(f"given Junction: {junction} is not connected to wire")

    def __repr__(self) -> str:
        return f"({id(self)})"


class Junction:
    def __init__(self):
        self.__wires: Set[Wire] = set()
        self.potential = []

    def connect(self, wires: Set[Wire]):
        self.__wires.update(wires)

    @property
    def wires(self) -> Set[Wire]:
        return self.__wires

    @property
    def is_null(self) -> bool:
        return len(self.__wires) == 0

    @property
    def is_singular(self) -> bool:
        return len(self.__wires) == 1

    def other_wires(self, wire: Wire) -> Set[Wire]:
        return self.__wires.difference({wire})


class Earth(Junction):
    def __init__(self):
        super().__init__()

    def get_wire(self, *args, **kwargs):
        earth_wire = Wire(*args, **kwargs)
        earth_wire.connect(self)
        self.connect({earth_wire})
        return earth_wire
