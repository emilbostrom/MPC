function X=blockrepeat(A,N)
% Skapar en blockdiagonal matris med samma block 
% repeterat N g�nger 

X = kron(eye(N),A);