from __future__ import annotations
from enum import IntEnum

from .devices import CircuitError
from .wire import Wire


class Sign(IntEnum):
    MINUS = -1
    ZERO = 0
    PLUS = 1


class SignedWires:
    def __init__(self):
        self.wires = []
        self.signs = []

    @staticmethod
    def shallow_copy(other: SignedWires) -> SignedWires:
        new_obj = SignedWires()
        new_obj.wires = [wire for wire in other.wires]
        new_obj.signs = [sign for sign in other.signs]
        return new_obj

    def append(self, wire: Wire, sign: Sign):
        if wire in self.wires:
            raise CircuitError(f"Wire: {wire} already exists in SignedWires: {self}")
        self.wires.append(wire)
        self.signs.append(sign)

    def get_sign(self, wire: Wire):
        return self.signs[self.wires.index(wire)]

    def __repr__(self) -> str:
        return f"({self.wires.__repr__()}, {self.signs.__repr__()})"
