function X=blockrepeat(A,N)
% Skapar en blockdiagonal matris med samma block 
% repeterat N gånger 

X = kron(eye(N),A);