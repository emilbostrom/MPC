load distilldata

C = eye(8);
D = zeros(8,2);

Ts = 1;

sys = ss(A,B,C,D);
sys_d = c2d(sys,Ts,'tustin');

F = sys_d.A;
G = sys_d.B;

N = 30;

Q1 = [1 0;
      0 15];

Q2 = 0.001*[1 0;
            0 1];

ubounds = [0 5;
           0 5];

global uold
uold = [0 0]';

xcross = [0.98 0.90 0.76 0.53 0.37 0.20 0.08 0.02];
%x0 = xcross;
x0 = xstar;

sim('distillmpcbasic.mdl')