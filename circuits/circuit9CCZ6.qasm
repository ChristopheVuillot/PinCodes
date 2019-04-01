OPENQASM 2.0;
include "qelib1.inc";
qreg q[9];
ccz q[0], q[4], q[8];
ccz q[2], q[3], q[7];
ccz q[1], q[5], q[6];
ccz q[0], q[1], q[2];
ccz q[3], q[4], q[5];
ccz q[6], q[7], q[8];
