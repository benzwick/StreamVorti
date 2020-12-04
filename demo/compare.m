a=dlmread('square10x10.dx.dat');
A=sparse(a(:,1),a(:,2),a(:,3));
b=dlmread('mfem_square10x10.dx.dat');
B=sparse(b(:,1),b(:,2),b(:,3));
A-B
spy(A)
spy(B)
