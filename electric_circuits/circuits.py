from scipy import linalg, integrate
from random import sample
from .signed_wires import *
from .wire import *


class Circuit:
    def __init__(self):
        self.earth = Earth()
        self.wires: Set[Wire] = set()
        self.junctions: Set[Junction] = {self.earth}
        self.events = {}
        self.t = []

    def connect(self, junction: Junction, wires: Set[Wire]):
        """
        function that connects junction to wire and each wire to junction;
        adds the wires to wires field of Circuit object;
        """
        self.junctions.add(junction)
        self.wires.update(wires)
        junction.connect(wires)
        for wire in wires:
            wire.connect(junction)

    def __calculate_effective_wires(self):
        """
        sorts the wires into effective wires and null wires (which do not contribute to calculation);
        """
        self.effective_wires = set()
        self.null_wires = set()
        for wire in self.wires:
            if wire.is_null or len(Circuit.__get_loops(wire)) == 0:
                self.null_wires.add(wire)
            else:
                self.effective_wires.add(wire)

    @staticmethod
    def __get_loops(main_wire: Wire) -> List[SignedWires]:
        """
        returns all the closed loops with the main wire using __loops;
        main wire should not be null;
        """
        if main_wire.is_null:
            raise CircuitError(f"main_wire ({main_wire}) is null")
        loops = []
        loop_path = SignedWires()
        loop_path.append(main_wire, Sign.PLUS)
        Circuit.__loops(loops=loops, endpoint=main_wire.starting_junction,
                        current_junction=main_wire.ending_junction, current_wire=main_wire,
                        loop_path=loop_path, parsed_junctions=[main_wire.ending_junction])

        return loops

    @staticmethod
    def __loops(loops: list, endpoint: Junction,
                current_junction: Junction, current_wire: Wire,
                loop_path: SignedWires, parsed_junctions: list):
        """
        recursive function that calculates closed loops;
        """
        for new_wire in current_junction.other_wires(current_wire):
            if new_wire.is_null:
                continue
            elif current_junction is new_wire.starting_junction:
                sign = Sign.PLUS
            elif current_junction is new_wire.ending_junction:
                sign = Sign.MINUS
            else:
                raise CircuitError("sign could not be deduced")

            # new_loop_path must be copied and not referenced
            new_loop_path = SignedWires.shallow_copy(loop_path)
            new_loop_path.append(new_wire, sign)

            next_junction = new_wire.other_junction(current_junction)

            if next_junction is endpoint:
                loops.append(new_loop_path)
                continue
            elif next_junction in parsed_junctions:
                continue

            Circuit.__loops(loops=loops, endpoint=endpoint,
                            current_junction=next_junction, current_wire=new_wire,
                            loop_path=new_loop_path, parsed_junctions=[*parsed_junctions, next_junction])

    def __first_law_wires(self):
        """
        return the wires to be used in implementation of Kirchhoff's First Law;
        """
        equations = []
        possible_wire_combos = set()

        for junction in self.junctions:
            # remove wires that are null
            junction_wires = junction.wires - self.null_wires
            if len(junction_wires) in (0, 1):
                continue

            # if some wires are RC and others LCR, the LCR wires are neglected in calculating possible combos
            #      (i will be passed as argument in __derivatives), i.e. i is not a variable in LCR wires
            # if all are LCR effect_junction_wires are the same as junction_wires
            effective_junction_wires = set([wire for wire in junction_wires])
            if not all([wire.inductor.is_not_null() for wire in junction_wires]):
                for wire in junction_wires:
                    if wire.inductor.is_not_null():
                        effective_junction_wires.remove(wire)

            # check if these wire can be derived from existing wire loops
            if effective_junction_wires in possible_wire_combos:
                continue

            possible_wire_combos.add(frozenset(effective_junction_wires))
            possible_wire_combos.update(
                set(wire_combo.symmetric_difference(effective_junction_wires) for wire_combo in possible_wire_combos))

            loop = SignedWires()
            for wire in junction_wires:
                loop.append(wire, Sign.PLUS if junction is wire.starting_junction else Sign.MINUS)
            equations.append(loop)

        for wires in equations:
            for wire in wires.wires:
                assert wire in self.effective_wires

        return equations

    def __second_law_wires(self):
        """
        return the wires to be used in implementation of Kirchhoff's Second Law;
        """
        equations = []
        parsed_wires = set()

        for wire in self.effective_wires:
            for loop in Circuit.__get_loops(wire):
                if set(loop.wires).issubset(parsed_wires):
                    continue
                parsed_wires.update(set(loop.wires))
                equations.append(loop)

        for wires in equations:
            for wire in wires.wires:
                assert wire in self.effective_wires

        return equations

    def __set_earth_potential(self, time: List[int]):
        self.earth.t = [t for t in time]
        self.earth.potential = [0.0] * len(time)

    def __calculate_potentials(self, current_junction: Junction, calculated: List[Junction]):
        """
        calculate potentials of junctions based on initial earthed junction;
        """
        # wires that have not already been calculated
        not_calc_wires = []
        for wire in current_junction.wires:
            next_junction = wire.other_junction(current_junction)
            if next_junction in calculated:
                continue

            not_calc_wires.append(wire)
            calculated.append(next_junction)

            for i in range(len(self.t)):
                time = self.t[i]
                current_junction_potential = current_junction.potential[i]
                potential_drop = (wire.battery.emf(time) + wire.ac_battery.emf(time)
                                  - wire.dq_dt[i] * wire.net_resistance(time)
                                  - wire.capacitor.potential_drop(time, wire.q[i]))
                # for inductor wire
                if wire.inductor.is_not_null():
                    potential_drop -= wire.inductor.potential_drop(time, wire.d2q_dt2[i])

                sign = Sign.PLUS if current_junction is wire.starting_junction else Sign.MINUS
                next_junction.potential.append(current_junction_potential + sign * potential_drop)

        for wire in not_calc_wires:
            self.__calculate_potentials(wire.other_junction(current_junction), calculated)

    def add_event(self, time: float, event: Callable[[Any], Any], args=None, kwargs=None):
        """
        when solving, the event is called using args and kwargs at specified time;
        multiple events are allowed at same time, but add_event must be called multiple times;
        """
        if args is None:
            args = tuple()
        if kwargs is None:
            kwargs = dict()

        if time in self.events.keys():
            self.events[time].append((event, args, kwargs))
        else:
            self.events[time] = [(event, args, kwargs)]

    def __split_events(self, end: float, resolution: int):
        """
        event = {t0 = [events, ...], t1 = [events, ...], ...};
        splits time interval into smaller intervals separated by events (function calls);
        """
        if len(self.events) == 0:
            return [((0.0, end),
                     [(self.__set_earth_potential, ([end * i / (resolution - 1) for i in range(resolution)],), {}),
                      (self.__calculate_potentials, (self.earth, []), {})])]

        times = sorted([i for i in list(self.events.keys()) if i <= end])
        split_times = []
        old = 0.0
        for i in range(len(times)):
            split_times.append(((old, times[i]), self.events[times[i]]))
            old = times[i]
        split_times.append(((times[-1], end),
                            [(self.__set_earth_potential,
                              ([end * i / (resolution - 1) for i in range(resolution)],), {}),
                             (self.__calculate_potentials, (self.earth, []), {})]))

        return split_times

    @staticmethod
    def __derivatives(charges: List[float], time: float, wires: List[Wire],
                      first_law: List[SignedWires], second_law: List[SignedWires]):
        """
        charges: [q0, i0, q1, q2, q3, i3, ...];
        qk, ik is for wires with inductor, qk is for no inductor;
        p = [<list_of_wires>, <first_law_wires>, <second_law_wires>];
        """
        derivatives = [0.0] * len(charges)

        # _charges groups the inputs: [(q0, i0), (q1,), (q2,), (q3, i3), ...]
        _charges = []
        index = 0
        for wire in wires:
            # if LCR wire
            if wire.inductor.is_not_null():
                _charges.append((charges[index], charges[index + 1]))
                index += 1
            # if RC wire
            else:
                _charges.append((charges[index],))
            index += 1

        # di/dt for LCR and i for RC is calculated using matrices
        # R i = V
        R = []
        V = []

        # fill R, V using first law
        for loop in first_law:
            R_eq = [0.0] * len(wires)
            V_eq = 0.0

            # if all wires are LCR then use di/dt instead of i
            if all([wire.inductor.is_not_null() for wire in loop.wires]):
                for wire in loop.wires:
                    R_eq[wires.index(wire)] = float(loop.get_sign(wire))
            else:
                for wire in loop.wires:
                    # if LCR
                    if wire.inductor.is_not_null():
                        V_eq -= loop.get_sign(wire) * _charges[wires.index(wire)][1]
                    # if RC
                    else:
                        R_eq[wires.index(wire)] = float(loop.get_sign(wire))

            assert R_eq not in R and [-i for i in R_eq] not in R

            R.append(R_eq)
            V.append(V_eq)

        # fill R, V using second law
        for loop in second_law:
            R_eq = [0.0] * len(wires)
            V_eq = 0.0

            for wire in loop.wires:
                # if LCR
                if wire.inductor.is_not_null():
                    R_eq[wires.index(wire)] = loop.get_sign(wire) * wire.inductor.inductance(time)
                    V_eq += loop.get_sign(wire) * (
                            wire.battery.emf(time) + wire.ac_battery.emf(time)
                            - wire.capacitor.potential_drop(time, _charges[wires.index(wire)][0])
                            - wire.net_resistance(time) * _charges[wires.index(wire)][1]
                    )
                # if RC
                else:
                    R_eq[wires.index(wire)] = loop.get_sign(wire) * wire.net_resistance(time)
                    V_eq += loop.get_sign(wire) * (
                            wire.battery.emf(time) + wire.ac_battery.emf(time)
                            - wire.capacitor.potential_drop(time, _charges[wires.index(wire)][0])
                    )

            R.append(R_eq)
            V.append(V_eq)

        i = linalg.solve(R, V)

        # fill rest of derivatives
        i1 = i2 = 0
        for wire in wires:
            # if LCR
            if wire.inductor.is_not_null():
                # dq/dt is passed in argument
                derivatives[i1] = _charges[i2][1]
                derivatives[i1 + 1] = i[i2]
                i1 += 1
            # if RC
            else:
                derivatives[i1] = i[i2]
            i1 += 1
            i2 += 1

        return derivatives

    def solve(self, end: float, resolution: int):
        dt = end / (resolution - 1)
        self.__calculate_effective_wires()

        # using sample is not required but is useful for detecting bugs that occur only in some iter
        wires = sample(list(self.effective_wires), len(self.effective_wires))

        for event in self.__split_events(end, resolution):
            # calculate t
            _start, _end = event[0]
            steps = int(resolution * (_end - _start) / end)
            t = [_start + (_end - _start) * i / (steps - 1) for i in range(steps)]
            self.t.extend(t)

            # recalculate first_law and second_law
            first_law = self.__first_law_wires()
            second_law = self.__second_law_wires()
            args = (wires, first_law, second_law)

            # initial conditions are last values of charge and current
            init_conditions = []
            for wire in wires:
                init_conditions.append(wire.q[-1] if len(wire.q) != 0 else wire.capacitor.init_charge)
                # if LCR
                if wire.inductor.is_not_null():
                    init_conditions.append(wire.dq_dt[-1] if len(wire.dq_dt) != 0 else 0.0)

            if len(first_law) != 0 and len(second_law) != 0:
                # calulate solution and derivative
                _sol = integrate.odeint(Circuit.__derivatives, init_conditions, t, args=args)
                sol = list(zip(*_sol))
                sol_derivatives = list(zip(*[Circuit.__derivatives(list(i), t, *args)
                                             for t, i in list(zip(t, _sol))]))

                # save solutions in wires
                index = 0
                for wire in wires:
                    wire.q.extend(sol[index])
                    wire.dq_dt.extend(sol_derivatives[index])
                    # if LCR
                    if wire.inductor.is_not_null():
                        wire.d2q_dt2.extend(sol_derivatives[index + 1])
                        index += 1
                    index += 1

            index = 0
            for wire in self.null_wires:
                init_charge = wire.q[-1] if len(wire.q) != 0 else wire.capacitor.init_charge
                wire.q.extend([init_charge] * len(t))
                wire.dq_dt.extend([0.0] * len(t))  # no current through null wire
                index += 1
                # if LCR
                if wire.inductor.is_not_null():
                    wire.d2q_dt2.extend([0.0] * len(t))
                    index += 1

            # call event functions
            for func in event[1]:
                func[0](*func[1], **func[2])

            # recalculate effective wires
            self.__calculate_effective_wires()
            wires = sample(list(self.effective_wires), len(self.effective_wires))

        count = len(self.t)

        # calculate other quantities
        for wire in self.wires:
            for i in range(count):
                time = self.t[i]
                potential_drop = (wire.battery.emf(time) + wire.ac_battery.emf(time)
                                  - wire.dq_dt[i] * wire.net_resistance(time)
                                  - wire.capacitor.potential_drop(time, wire.q[i]))
                # for inductor wire
                if wire.inductor.is_not_null():
                    potential_drop -= wire.inductor.potential_drop(time, wire.d2q_dt2[i])

                wire.potential_drop.append(potential_drop)

            wire.eq_resistance = [wire.potential_drop[i] / wire.dq_dt[i] if wire.dq_dt[i] != 0.0 else 0.0
                                  for i in range(count)]
            wire.eq_capacitance = [wire.q[i] / wire.potential_drop[i] if wire.potential_drop[i] != 0.0 else 0.0
                                   for i in range(count)]
            if wire.inductor.is_not_null():
                wire.eq_inductance = [wire.potential_drop[i] / wire.d2q_dt2[i] if wire.d2q_dt2[i] != 0.0 else 0.0
                                      for i in range(count)]

            wire.power = [wire.dq_dt[i] * wire.dq_dt[i] * wire.net_resistance(self.t[i])
                          for i in range(count)]
            wire.heat = [sum([p * dt for p in wire.power[:i]]) for i in range(count)]

            if wire.inductor.is_not_null():
                wire.inductor_energy = [0.5 * wire.inductor.inductance(self.t[i]) * wire.dq_dt[i] * wire.dq_dt[i]
                                        for i in range(count)]
            wire.capacitor_energy = [0.5 * wire.q[i] * wire.capacitor.potential_drop(self.t[i], wire.q[i])
                                     for i in range(count)]

            wire.galvanometer_reading = [wire.galvanometer.get_reading(i) for i in wire.dq_dt]
            wire.ammeter_reading = [wire.ammeter.get_reading(i) for i in wire.dq_dt]
            wire.voltmeter_reading = [wire.voltmeter.get_reading(self.t[i], wire.dq_dt[i])
                                      for i in range(count)]
