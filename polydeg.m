function PolDeg = polydeg(d)
% -------------------------------------------------------------------------
%
% INPUT:
% d       - degree
% OUTPUT
% PolDeg  - matrix Nx2 with all the multiindices of the degrees of the 
%           polynomials of degree <= d
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
i = 1;
PolDeg(i,:) = [0,0];

while PolDeg(i,1) ~= d
    x = PolDeg(i,:);
    PolDeg(i+1,:) = mono_next_grlex(2,x);
    i = i+1;
end

end