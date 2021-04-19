load distilldata

A = [A zeros(8,1);
     zeros(1,9)];
 
B =[B zeros(8,1); zeros(1,2) -1/3];

C = eye(9);
 
D = zeros(9,3);

Ts = 1;

sys = ss(A,B,C,D);
sys_d = c2d(sys,Ts,'tustin');

F = sys_d.A;
G = sys_d.B;

N = 30;

Q1 = [1 0  0;
      0 15 0;
      0 0  0.001];

Q2 = 0.001*[1 0 0;
            0 1 0;
            0 0 1000];
             
ubounds = [0 5;
           0 5;
           0 5];

Mz = [Mz zeros(2,1);
      zeros(1,8) 1];
       
global uold
uold = [0 0 0]';

xcross = [0.98 0.90 0.76 0.53 0.37 0.20 0.08 0.02];
%x0 = xcross;
x0 = xstar;

sim('distillmpcbuffer.mdl')