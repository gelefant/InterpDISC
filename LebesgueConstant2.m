function LebCon = LebesgueConstant2(deg,centers,radii)
% -------------------------------------------------------------------------
% Function to approximate the Lebesgue constant of the interpolant constructed
% via integrals on discs, it requires the degree of the interpolant, the
% centers and the radii.
% It considers M number of random points in the discs where compute the
% argument of the supremum on the definition of the Lebesgue constant. 
% 
% INPUT:
% deg     - degree of the interpolant
% centers - a matrix Nx2 of the coordinates of the centers in the unitarian
%           disc
% radii   - a column vector di dimension N of the radii of the discs
% OUTPUT
% LebCon  - an approximation of the value of the Lebesgue constant
% -------------------------------------------------------------------------
% Dates
%--------------------------------------------------------------------------
% First version: November 15, 2023;
% Checked: October 08, 2024.
%--------------------------------------------------------------------------
% Authors
%--------------------------------------------------------------------------
% L. Bruni Bruno and G. Elefante
%--------------------------------------------------------------------------
% Paper
%--------------------------------------------------------------------------
% "Uniform approximation of diffused data"
% L. Bruni Bruno and G. Elefante
%--------------------------------------------------------------------------

% Number of repetition due to the randomness
rep = 50;

% Number of random disk where to evaluate 
M = 10001;

% Number of random disk where to evaluate 
M = 10001;

S = diag(radii.^2)*pi;
Vcheb = chebVandInt(deg,centers,radii);
U = Vcheb\S;

% Centers on halton points moved to the disk
for h = 1:rep
    EvCent = rand(M,2);
    EvCent = EvCent*2-1;
    EvCent = [EvCent(:,1).*sqrt(1-EvCent(:,2).^2/2),EvCent(:,2).*sqrt(1-EvCent(:,1).^2/2)];
    W = chebVand(deg,EvCent);
    Leb(h) = norm(W*U,'inf');
end
LebCon = max(Leb);