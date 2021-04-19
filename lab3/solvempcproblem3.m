function un = solvempcproblem3(x,F,G,M,N,Q1,Q2,ubounds,Ts,t)

global uold

[H,S] = createpredictors(F,G,N);
Q1b = blockrepeat(Q1,N);
Q2b = blockrepeat(Q2,N);
Mb = blockrepeat(M,N);
[n,m]=size(G);
simtim = 100;

Au = [eye(N*m);
      -eye(N*m)];

bu = [repmat(ubounds(:,2),N,1);repmat(-ubounds(:,1),N,1)];

n1 = round(20/Ts);
n2 = round(60/Ts)-n1;
n3 = round(simtim/Ts)-n1-n2;

r1 = [0*ones(n1,1);0*ones(n2,1);0*ones(n3+N,1)];
r2 = [0*ones(n1,1);0*ones(n2,1);0*ones(n3+N,1)];
r3 = [0*ones(n1,1);0*ones(n2,1);0*ones(n3+N,1)];
Rlong = zeros(round(simtim/Ts)*m+N*m,1);
Rlong(1:m:end-2) = r1;
Rlong(2:m:end-1) = r2;
Rlong(3:m:end) = r3;

k=round(t/Ts);
R = Rlong(k*m+1:k*m+N*m);

omega = eye(N*m)-[zeros(m,N*m);eye((N-1)*m) zeros((N-1)*m,m)];


delta = [uold;
         zeros((N-1)*m,1)];

    
U = quadprog(S'*Mb'*Q1b*Mb*S+omega'*Q2b*omega,S'*Mb'*Q1b*(Mb*H*x-R)-omega'*Q2b*delta,Au,bu);
%U = quadprog(S'*Mb'*Q1b*Mb*S+Q2b,S'*Mb'*Q1b*Mb*H*x,Au,bu);

un = U(1:m);
%un = [0 0];

uold = un;

end