from urllib import response
from matplotlib.pyplot import thetagrids
from numpy import identity
from qutip.qip.circuit import QubitCircuit, Gate, CircuitSimulator
from qutip.qip.qasm import read_qasm
from qutip import tensor, basis
from qutip.states import ket2dm
from qutip.measurement import measure, measurement_statistics
from qutip.operators import identity, sigmaz

import math
import cmath

alpha = 0
beta = 0
theta1 = 0
phi1 = 0
theta2 = 0
phi2 = 0

qc = QubitCircuit(N=2)
qc.add_gate("SNOT", targets=[0])
qc.add_gate("SNOT", targets=[1])
qc.add_gate("CSIGN", targets=[1], controls=[0])
qc.add_gate("RZ", targets=[0], arg_value=alpha)
qc.add_gate("RZ", targets=[1], arg_value=beta)
qc.add_gate("SNOT", targets=[0])
qc.add_gate("SNOT", targets=[1])

init_state = tensor(basis(2, 0), basis(2, 0))

print("### INIT ###")
print(init_state)

result = qc.run(state=init_state)
print("### RES ###")
print(result)

state0 = basis(2, 0)
state1 = basis(2, 1)

print("### states ###")
print(state0)
print(state1)

z1 = complex(0, phi1)
print("### Z ###")
print(z1)

E1 = (
    math.cos(theta1 / 2) * state0 + cmath.exp(z1) * math.sin(theta1 / 2) * state1
).unit()
E1_orth = (
    -math.cos(theta1 / 2) * state1 + cmath.exp(-z1) * math.sin(theta1 / 2) * state0
).unit()

print("### E1 ###")
print(E1)
print(E1_orth)

z2 = complex(0, phi2)
E2 = (
    math.cos(theta2 / 2) * state0 + cmath.exp(z2) * math.sin(theta2 / 2) * state1
).unit()
E2_orth = (
    -math.cos(theta2 / 2) * state1 + cmath.exp(-z2) * math.sin(theta2 / 2) * state0
).unit()

matrix = tensor((ket2dm(E1) - ket2dm(E1_orth)), (ket2dm(E2) - ket2dm(E2_orth)))

print("### DM ###")
print(ket2dm(E1))
print(ket2dm(E2))

print("### Matrix ###")
print(matrix)

results = {"00": 0, "01": 0, "10": 0, "11": 0}

for _ in range(100):
    value, new_state = measure(result, matrix)
    print("Value + New State")
    print(value)
    print(new_state)
    if value == 1:
        mx = tensor(sigmaz(), identity(2))
        print("mx")
        print(mx)
        val, nstate = measure(new_state, mx)
        print(val)
        print(nstate)
        if val == 1:
            results["00"] += 1
        elif val == -1:
            results["11"] += 1
        else:
            print("Unexpected error.")
    elif value == -1:
        mx = tensor(sigmaz(), identity(2))
        print("mx")
        print(mx)
        val, nstate = measure(new_state, mx)
        print(val)
        print(nstate)
        if val == 1:
            results["01"] += 1
        elif val == -1:
            results["10"] += 1
        else:
            print("Unexpected error.")

print(results)
