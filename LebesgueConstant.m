function LebCon = LebesgueConstant(deg,centers,radii)
% -------------------------------------------------------------------------
% Function to approximate the Lebesgue constant of the interpolant constructed
% via integrals on discs, it requires the degree of the interpolant, the
% centers and the radii.
% It considers M number of points as centers of the discs where compute the
% argument of the supremum on the definition of the Lebesgue constant. It 
% is computed by using a random radius, hence it is repeted rep times and, 
% then it is considered the maximum over these values 
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
% Checked: December 07, 2023.
%--------------------------------------------------------------------------
% Authors
%--------------------------------------------------------------------------
% L. Bruni Bruno and G. Elefante
%--------------------------------------------------------------------------
% Paper
%--------------------------------------------------------------------------
% "Interpolation by integrals on discs"
% L. Bruni Bruno and G. Elefante
%--------------------------------------------------------------------------

% Number of repetition due to the randomness
rep = 50;

% Number of random disk where to evaluate 
M = 5001;

% Centers on halton points moved to the disk
p = sobolset(2);
EvCent = net(p,M);
EvCent = EvCent*2-1; EvCent(1,:) = [];
EvCent = [EvCent(:,1).*sqrt(1-EvCent(:,2).^2/2),EvCent(:,2).*sqrt(1-EvCent(:,1).^2/2)];
% Evaluation of the distance from the center
NormCent = vecnorm(EvCent')';
% Vector containing the maximum value of the radius for each disk
RadMax = (1-NormCent)/2;
% Matrix of the 
S = diag(radii.^2);
Vcheb = chebVandInt(deg,centers,radii);
U = Vcheb\S;

for h = 1:rep
    REv = RadMax.*rand(size(NormCent));
    C = diag(1./(REv).^2);
    W = chebVandInt(deg,EvCent,REv);
    Leb(h) = norm(C*W*U,'inf');
end
LebCon = max(Leb);