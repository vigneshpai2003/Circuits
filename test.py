from electric_circuits import *
from matplotlib import pyplot


def main1():
    circuit = Circuit()
    wires = [Wire(lambda t: 1.5 + sin(t), battery=2.0, inductor=0.2),
             Wire(0.1, ac_battery=(3.0, 5.0, 0.5)),
             Wire(1.1, capacitor=(lambda t: 1.5 + sin(t), -10.0), inductor=1.0),
             Wire(1.0, inductor=lambda t: 1.5 + sin(t), ammeter=(0.01, 2))]
    junctions = [Junction(), Junction(), Junction()]

    circuit.connect(junctions[0], {wires[0], wires[2], wires[3]})
    circuit.connect(junctions[1], {wires[0], wires[1]})
    circuit.connect(junctions[2], {wires[1], wires[2], wires[3], circuit.earth.get_wire(0.0)})

    wires[0].swap_junctions()

    def toggle(_wire): return wires[_wire].switch.toggle()

    circuit.add_event(3.0, event=toggle, args=(3, ))
    circuit.add_event(7.0, event=toggle, args=(2, ))
    circuit.add_event(7.0, event=toggle, args=(3, ))
    circuit.add_event(8.0, event=toggle, args=(2, ))

    end = 20
    circuit.solve(end, 1000)

    pyplot.figure(num="Electric Circuits")
    pyplot.title("Complex LCR Circuit")
    pyplot.xlabel("time (s)")
    pyplot.xticks(range(end + 1))

    pyplot.plot(circuit.t, wires[2].q)
    pyplot.plot(circuit.t, wires[2].dq_dt)
    pyplot.plot(circuit.t, wires[2].d2q_dt2)
    pyplot.plot(circuit.t, wires[2].power)
    # pyplot.plot(circuit.t, junctions[1].potential)

    pyplot.legend(['Q', 'dQ/dt', 'd2Q/dt2', 'power'])
    pyplot.grid(True)
    pyplot.show()


if __name__ == '__main__':
    main1()
