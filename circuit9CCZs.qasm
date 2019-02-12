OPENQASM 2.0;
include "qelib1.inc";
qreg q[9];
gate ccz2t a,b,c
{
 h c;
 s c;
 z c;
 h c;
 x c;
 cx c, b;
 cx c, a;
 x c;
 h c;
 cz c, b;
 cz c, a;
 h c;
 t c;
 s c;
 z c;
 h c;
 cz c, b;
 cz c, a;
 h c;
}
ccz q[0], q[4], q[8];
ccz q[2], q[3], q[7];
ccz q[1], q[5], q[6];

ccz2t q[0], q[1], q[2];
ccz2t q[3], q[4], q[5];
ccz2t q[6], q[7], q[8];
