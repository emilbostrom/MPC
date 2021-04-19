function [H,S]=createpredictors(F,G,N);
% Skapa H och S i prediktionen X=Hx(k)+SU
% Systemdimension
[n,m] = size(G);
% H
Fprod=eye(n);
H = zeros(n*N,n);
for j=1:N,
    H(1+n*(j-1):n*j,:)=F^(j-1);
end;
% S
S = zeros(N*n,N*m);
for j=2:N,
    for k=1:j-1,
        S(1+n*(j-1):n*j,1+m*(k-1):k*m)=F^(j-k-1)*G;
    end;
end;
