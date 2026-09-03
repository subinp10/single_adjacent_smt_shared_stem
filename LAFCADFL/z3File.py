from z3 import *

s = Optimize()

# Variable declaration
A_6_4_R3 = Int("A_6_4_R3")
A_6_4_R4 = Int("A_6_4_R4")
A_6_5_R3 = Int("A_6_5_R3")
A_6_5_R4 = Int("A_6_5_R4")
A_7_4_R3 = Int("A_7_4_R3")
A_7_4_R4 = Int("A_7_4_R4")
A_7_5_R3 = Int("A_7_5_R3")
A_7_5_R4 = Int("A_7_5_R4")

V_R3 = Int("V_R3")
V_R4 = Int("V_R4")

d_0_R3  = Int("d_0_R3")
d_0_R4  = Int("d_0_R4")
d_1_R3  = Int("d_1_R3")
d_1_R4  = Int("d_1_R4")
d_2_R3  = Int("d_2_R3")
d_2_R4  = Int("d_2_R4")

d_6_4_0_R3 = Int("d_6_4_0_R3")
d_6_4_1_R3 = Int("d_6_4_1_R3")
d_6_4_2_R3 = Int("d_6_4_2_R3")
d_6_4_0_R4 = Int("d_6_4_0_R4")
d_6_4_1_R4 = Int("d_6_4_1_R4")
d_6_4_2_R4 = Int("d_6_4_2_R4")
d_6_5_0_R3 = Int("d_6_5_0_R3")
d_6_5_1_R3 = Int("d_6_5_1_R3")
d_6_5_2_R3 = Int("d_6_5_2_R3")
d_6_5_0_R4 = Int("d_6_5_0_R4")
d_6_5_1_R4 = Int("d_6_5_1_R4")
d_6_5_2_R4 = Int("d_6_5_2_R4")
d_7_4_0_R3 = Int("d_7_4_0_R3")
d_7_4_1_R3 = Int("d_7_4_1_R3")
d_7_4_2_R3 = Int("d_7_4_2_R3")
d_7_4_0_R4 = Int("d_7_4_0_R4")
d_7_4_1_R4 = Int("d_7_4_1_R4")
d_7_4_2_R4 = Int("d_7_4_2_R4")
d_7_5_0_R3 = Int("d_7_5_0_R3")
d_7_5_1_R3 = Int("d_7_5_1_R3")
d_7_5_2_R3 = Int("d_7_5_2_R3")
d_7_5_0_R4 = Int("d_7_5_0_R4")
d_7_5_1_R4 = Int("d_7_5_1_R4")
d_7_5_2_R4 = Int("d_7_5_2_R4")

# Number of cells in the mixture filled with Rt reagent
s.add(And(A_6_4_R3>=0, A_6_4_R3<=1, A_6_4_R4>=0, A_6_4_R4<=1))
s.add(And(A_6_5_R3>=0, A_6_5_R3<=1, A_6_5_R4>=0, A_6_5_R4<=1))
s.add(And(A_7_4_R3>=0, A_7_4_R3<=1, A_7_4_R4>=0, A_7_4_R4<=1))
s.add(And(A_7_5_R3>=0, A_7_5_R3<=1, A_7_5_R4>=0, A_7_5_R4<=1))

s.add(A_6_4_R3 + A_6_5_R3 + A_7_4_R3 + A_7_5_R3 == V_R3)
s.add(A_6_4_R4 + A_6_5_R4 + A_7_4_R4 + A_7_5_R4 == V_R4)

# If a cell is filled with reagent Rt then no reagent Rk where k != t can be filled in that cell.
s.add(A_6_4_R3 + A_6_4_R4 == 1)
s.add(A_6_5_R3 + A_6_5_R4 == 1)
s.add(A_7_4_R3 + A_7_4_R4 == 1)
s.add(A_7_5_R3 + A_7_5_R4 == 1)

# To get traceability and connectivity.
s.add(If(And(A_6_4_R3 == 1, A_6_5_R3 + A_7_4_R3 == 0), (d_6_4_0_R3 == 1), (d_6_4_0_R3 == 0)))
s.add(If(And(A_6_5_R3 == 1, A_7_5_R3 + A_6_4_R3 == 0), (d_6_5_0_R3 == 1), (d_6_5_0_R3 == 0)))
s.add(If(And(A_7_4_R3 == 1, A_6_4_R3 + A_7_5_R3 == 0), (d_7_4_0_R3 == 1), (d_7_4_0_R3 == 0)))
s.add(If(And(A_7_5_R3 == 1, A_6_5_R3 + A_7_4_R3 == 0), (d_7_5_0_R3 == 1), (d_7_5_0_R3 == 0)))
s.add(If(And(A_6_4_R4 == 1, A_6_5_R4 + A_7_4_R4 == 0), (d_6_4_0_R4 == 1), (d_6_4_0_R4 == 0)))
s.add(If(And(A_6_5_R4 == 1, A_7_5_R4 + A_6_4_R4 == 0), (d_6_5_0_R4 == 1), (d_6_5_0_R4 == 0)))
s.add(If(And(A_7_4_R4 == 1, A_6_4_R4 + A_7_5_R4 == 0), (d_7_4_0_R4 == 1), (d_7_4_0_R4 == 0)))
s.add(If(And(A_7_5_R4 == 1, A_6_5_R4 + A_7_4_R4 == 0), (d_7_5_0_R4 == 1), (d_7_5_0_R4 == 0)))
s.add(If(And(A_6_4_R3 == 1, A_6_5_R3 + A_7_4_R3 == 1), (d_6_4_1_R3 == 1), (d_6_4_1_R3 == 0)))
s.add(If(And(A_6_5_R3 == 1, A_7_5_R3 + A_6_4_R3 == 1), (d_6_5_1_R3 == 1), (d_6_5_1_R3 == 0)))
s.add(If(And(A_7_4_R3 == 1, A_6_4_R3 + A_7_5_R3 == 1), (d_7_4_1_R3 == 1), (d_7_4_1_R3 == 0)))
s.add(If(And(A_7_5_R3 == 1, A_6_5_R3 + A_7_4_R3 == 1), (d_7_5_1_R3 == 1), (d_7_5_1_R3 == 0)))
s.add(If(And(A_6_4_R4 == 1, A_6_5_R4 + A_7_4_R4 == 1), (d_6_4_1_R4 == 1), (d_6_4_1_R4 == 0)))
s.add(If(And(A_6_5_R4 == 1, A_7_5_R4 + A_6_4_R4 == 1), (d_6_5_1_R4 == 1), (d_6_5_1_R4 == 0)))
s.add(If(And(A_7_4_R4 == 1, A_6_4_R4 + A_7_5_R4 == 1), (d_7_4_1_R4 == 1), (d_7_4_1_R4 == 0)))
s.add(If(And(A_7_5_R4 == 1, A_6_5_R4 + A_7_4_R4 == 1), (d_7_5_1_R4 == 1), (d_7_5_1_R4 == 0)))
s.add(If(And(A_6_4_R3 == 1, A_6_5_R3 + A_7_4_R3 == 2), (d_6_4_2_R3 == 1), (d_6_4_2_R3 == 0)))
s.add(If(And(A_6_5_R3 == 1, A_7_5_R3 + A_6_4_R3 == 2), (d_6_5_2_R3 == 1), (d_6_5_2_R3 == 0)))
s.add(If(And(A_7_4_R3 == 1, A_6_4_R3 + A_7_5_R3 == 2), (d_7_4_2_R3 == 1), (d_7_4_2_R3 == 0)))
s.add(If(And(A_7_5_R3 == 1, A_6_5_R3 + A_7_4_R3 == 2), (d_7_5_2_R3 == 1), (d_7_5_2_R3 == 0)))
s.add(If(And(A_6_4_R4 == 1, A_6_5_R4 + A_7_4_R4 == 2), (d_6_4_2_R4 == 1), (d_6_4_2_R4 == 0)))
s.add(If(And(A_6_5_R4 == 1, A_7_5_R4 + A_6_4_R4 == 2), (d_6_5_2_R4 == 1), (d_6_5_2_R4 == 0)))
s.add(If(And(A_7_4_R4 == 1, A_6_4_R4 + A_7_5_R4 == 2), (d_7_4_2_R4 == 1), (d_7_4_2_R4 == 0)))
s.add(If(And(A_7_5_R4 == 1, A_6_5_R4 + A_7_4_R4 == 2), (d_7_5_2_R4 == 1), (d_7_5_2_R4 == 0)))

s.add(d_6_4_0_R3 + d_6_5_0_R3 + d_7_4_0_R3 + d_7_5_0_R3 == d_0_R3)
s.add(d_6_4_1_R3 + d_6_5_1_R3 + d_7_4_1_R3 + d_7_5_1_R3 == d_1_R3)
s.add(d_6_4_2_R3 + d_6_5_2_R3 + d_7_4_2_R3 + d_7_5_2_R3 == d_2_R3)
s.add(d_6_4_0_R4 + d_6_5_0_R4 + d_7_4_0_R4 + d_7_5_0_R4 == d_0_R4)
s.add(d_6_4_1_R4 + d_6_5_1_R4 + d_7_4_1_R4 + d_7_5_1_R4 == d_1_R4)
s.add(d_6_4_2_R4 + d_6_5_2_R4 + d_7_4_2_R4 + d_7_5_2_R4 == d_2_R4)

s.add(Implies(V_R3 == 1, And(d_0_R3 == 1, d_1_R3 == 0, d_2_R3 == 0)))
s.add(Implies(V_R3 == 2, And(d_0_R3 == 0, d_1_R3 == 2, d_2_R3 == 0)))
s.add(Implies(V_R3 == 3, And(d_0_R3 == 0, d_1_R3 == 2, d_2_R3 == 1)))
s.add(Implies(V_R3 == 4, And(d_0_R3 == 0, d_1_R3 == 0, d_2_R3 == 4)))
s.add(Implies(V_R4 == 1, And(d_0_R4 == 1, d_1_R4 == 0, d_2_R4 == 0)))
s.add(Implies(V_R4 == 2, And(d_0_R4 == 0, d_1_R4 == 2, d_2_R4 == 0)))
s.add(Implies(V_R4 == 3, And(d_0_R4 == 0, d_1_R4 == 2, d_2_R4 == 1)))
s.add(Implies(V_R4 == 4, And(d_0_R4 == 0, d_1_R4 == 0, d_2_R4 == 4)))

s.add(And(V_R3 == 3, V_R4 == 1))

if s.check() == unsat:
	print("Not possible to create traceable graph for all reagents")
else:
	fp = open('output.txt','w')
	values = s.model()
	if values[A_6_4_R3] == 1:
		fp.write("6,4,R3\n")
	if values[A_6_4_R4] == 1:
		fp.write("6,4,R4\n")
	if values[A_6_5_R3] == 1:
		fp.write("6,5,R3\n")
	if values[A_6_5_R4] == 1:
		fp.write("6,5,R4\n")
	if values[A_7_4_R3] == 1:
		fp.write("7,4,R3\n")
	if values[A_7_4_R4] == 1:
		fp.write("7,4,R4\n")
	if values[A_7_5_R3] == 1:
		fp.write("7,5,R3\n")
	if values[A_7_5_R4] == 1:
		fp.write("7,5,R4\n")
