function T=chebpolys(deg,x)

%--------------------------------------------------------------------------
% Object:
% This routine computes the Chebyshev-Vandermonde matrix on the real line
% by recurrence.
%--------------------------------------------------------------------------
% Input:
% deg: maximum polynomial degree
% x: 1-column array of abscissas
%--------------------------------------------------------------------------
% Output:
% T: Chebyshev-Vandermonde matrix at x, T(i,j+1)=T_j(x_i), j=0,...,deg.
%--------------------------------------------------------------------------
% Authors:
% Alvise Sommariva and Marco Vianello
% University of Padova, December 15, 2017
%--------------------------------------------------------------------------

T=zeros(length(x),deg+1);
t0=ones(length(x),1); T(:,1)=t0;
t1=x; T(:,2)=t1;

for j=2:deg
    t2=2*x.*t1-t0; 
    T(:,j+1)=t2;
    t0=t1; 
    t1=t2;
end