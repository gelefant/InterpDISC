function demo_vandermonde_orb
% -------------------------------------------------------------------------
% In this demo we compute the conditioning of the Vandermonde matrices 
% construct via monomial basis and Chebyshev polynomials and the Lebesgue 
% constant of the interpolant construct by integrals on disks with 
% centers lying on orbits 
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


%--------------------------------------------------------------------------
% Parameters
degSet = [1:20];
do_plot = 1;
intersection = 0; % 0 - no intersection, 1 - fixed radmax, 2 - caos 
typeset_rho = 4; % 1 - rand, 2 - ordered rand, 3 - equispaced, 4 - CL no extrema, 5 - Cheby
RadMaxDiscs = 0; % 0 - fixed radius, 1 - random radius
%--------------------------------------------------------------------------


h = 1;
for deg = degSet

fprintf('\n')
fprintf('...................................................')
fprintf('\n \t          degree   : %5d',deg)

dimP = nchoosek(deg+2,2);

NOrb = ceil((deg+1)/2); %floor(deg/2)+1;

dimPdisk = 2*[deg:-2:0] + 1;

centers = [];

switch typeset_rho
    case 1
        rho = rand(size(dimPdisk));
    case 2
        rho = sort(rand(size(dimPdisk)),'descend');
    case 3
        if dimPdisk(end)==1
            rho = linspace(0,1,size(dimPdisk,2)+1); 
            rho(end) = []; 
            rho = rho(end:-1:1);
        else
            rho = linspace(0,1,size(dimPdisk,2)+2); 
            rho(end) = []; rho(1) = [];
            rho = rho(end:-1:1);
        end
    case 4
        if dimPdisk(end)==1
            rho = cos([0:(2*NOrb)]*pi/(2*NOrb)); 
            rho = rho(2:NOrb+1);
        else
            rho = cos([0:(2*NOrb+1)]*pi/(2*NOrb+1));  
            rho = rho(2:NOrb+1);
        end
    case 5
        rho = cos((2*[1:(2*NOrb)]-1)*pi/(4*NOrb));
        rho = rho(1:NOrb);
end

for j = 1:length(dimPdisk)
    theta = linspace(0,2*pi,dimPdisk(j)+1); theta(end) = [];
    theta = mod(pi*j/3+theta,2*pi);
    centers = [centers, rho(j)*exp(1i*theta)];
end
centers = [real(centers'),imag(centers')];

C = exp(1i*linspace(0,2*pi,150));

if do_plot
    figure(1)
    plot(centers(:,1),centers(:,2),'.','Color',[245,119,34]/255); % 'Color',[147,147,147]/255
    hold on;
    plot(real(C),imag(C),'Color',[57,57,57]/255);
end

if intersection == 2
    Rmax = 1-rho';
elseif intersection == 1
    Rmax = 1-max(rho);
else
    if dimPdisk(end)==1
        R = [1,rho];
        Rmax = min(abs(R(1:end-1)-R(2:end)))/2;
    else
         R = [1,rho,0];
        Rmax = min(abs(R(1:end-1)-R(2:end)))/2;
    end
    Rmax = min(Rmax,2*pi/(2*dimPdisk(1)+2));
end

if RadMaxDiscs == 0
    radii_orb = Rmax.*ones(NOrb,1);
else
    radii_orb = Rmax*rand(NOrb,1);
end

radii = [];
for i = 1:NOrb
radii = [radii;repmat(radii_orb(i),dimPdisk(i),1)];
end

if do_plot
    for j=1:size(centers,1)
        CC = C*radii(j)+centers(j,1)+1i*centers(j,2);
        plot(real(CC),imag(CC),'Color',[147,147,147]/255) % 'Color',[245,119,34]/255
    end
end



Vmon = VandInt(deg,centers,radii);

Vcheb = chebVandInt(deg,centers,radii);

condVmon(h) = cond(Vmon);
condVcheb(h) = cond(Vcheb);

fprintf('\n...................................................')
fprintf('\n \t Cond Vand         : %2.3e',condVmon(h))
fprintf('\n \t Cond Vand cheb    : %2.3e',condVcheb(h))
fprintf('\n \t ratio cond        : %2.3e',condVmon(h)/condVcheb(h))
fprintf('\n...................................................')

% Lebesgue constant computation
LebCon(h) = LebesgueConstant2(deg,centers,radii);


fprintf('\n \t Leb constant      : %2.3e',LebCon(h))
fprintf('\n...................................................')
fprintf('\n')

h = h + 1;
end

if do_plot
    figure(2)
    semilogy(degSet,condVmon,'-o','Color',[0,142,142]/255,'LineWidth',1.5)
    hold on;
    semilogy(degSet,condVcheb,'-o','Color',[142,0,0]/255,'LineWidth',1.5)

    figure(3)
    semilogy(degSet,LebCon,'-o','Color',[142,0,0]/255,'LineWidth',1.5);
    hold on
    if typeset_rho == 4
        semilogy(degSet,LebCon(1)*degSet.^(1/2),'-o','Color',[142,142,142]/255,'LineWidth',1.5)
        semilogy(degSet,LebCon(1)*degSet,'-o','Color',[142,142,142]/255,'LineWidth',1.5)
    end
end


