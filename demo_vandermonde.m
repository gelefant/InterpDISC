function demo_vandermonde
% -------------------------------------------------------------------------
% In this demo we compute the conditioning of the Vandermonde matrices 
% construct via monomial basis and Chebyshev polynomials and the Lebesgue 
% constant of the interpolant construct by integrals on disks with 
% scattered centers 
% 
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
type_centers = 1; % 0 - halton pts, 1 - Meurant Sommariva pts, 2 - Carnicer Godes pts,  3 - mapped grid, 4 - cheby mesh
RadMaxDiscs = 0; % 0 - fixed radius, 1 - random radius
%--------------------------------------------------------------------------

p = haltonset(2);
h = 1;
fprintf('\n\n')

for deg = degSet
    
fprintf('...................................................')
fprintf('\n \t degree          : %5d \n',deg)

dimP = nchoosek(deg+2,2);

switch type_centers
    case 0
        centers = net(p,dimP+1); centers(1,:) = [];
        centers = centers*2 - 1;
        centers = [centers(:,1).*sqrt(1-centers(:,2).^2/2),centers(:,2).*sqrt(1-centers(:,1).^2/2)];
    case 1
        epsilon = 1e-4;
        [centers,stats_matrix]=set_disk_lebesgue(deg);
        centers = (1-epsilon)*centers;
    case 2
        epsilon = 1e-4;
        [centers,stats_matrix]=set_disk_cg(deg);
        centers = (1-epsilon)*centers;
    case 3
        if mod(deg+1,2)==0
            x1 = linspace(-1,1,(deg+1)/2+2); x1(1) = []; x1(end) = [];
            x2 = linspace(-1,1,deg+4); x2(1) = []; x2(end) = [];
        else
            x1 = linspace(-1,1,deg+3); x1(1) = []; x1(end) = [];
            x2 = linspace(-1,1,(deg+2)/2+2); x2(1) = []; x2(end) = [];
        end
        [xx1,xx2] = meshgrid(x1,x2); xx1 = xx1(:); xx2 = xx2(:);
        centers = [xx1.*sqrt(1-xx2.^2/2),xx2.*sqrt(1-xx1.^2/2)];
    case 4
        if mod(deg+1,2)==0
            x1 = cos((2*[1:(deg+1)/2]-1)*pi/(deg+1));
            x2 = cos((2*[1:(deg+2)]-1)*pi/(2*(deg+2)));
            x1 = x1/2+1/2;
            x2 = 2*asin(sin(pi/2)*x2)+pi;
            [xx1,xx2] = meshgrid(x1,x2);
            y1 = xx1.*cos(xx2);
            y2 = xx1.*sin(xx2);
            centers = [y1(:),y2(:)];
        else
            x1 = cos((2*[1:(deg+1)]-1)*pi/(2*(deg+1)));
            x2 = cos((2*[1:(deg+2)/2]-1)*pi/(deg+2));
            x1 = x1/2+1/2;
            x2 = 2*asin(sin(pi/2)*x2)+pi;
            [xx1,xx2] = meshgrid(x1,x2);
            y1 = xx1.*cos(xx2);
            y2 = xx1.*sin(xx2);
            centers = [y1(:),y2(:)];
        end
end

rho = vecnorm(centers,2,2);

[rho,idx] = sort(rho,'descend');

centers = centers(idx,:);

C = exp(1i*linspace(0,2*pi,150));

if do_plot
    figure(1)
    plot(centers(:,1),centers(:,2),'.','Color',[147,147,147]/255);
    axis equal
    hold on;
    plot(real(C),imag(C),'Color',[57,57,57]/255);
end

if intersection == 2
    Rmax = 1-rho;
    if RadMaxDiscs==0 
        radii = Rmax.*ones(size(rho));
    else 
        radii = Rmax.*rand(size(rho));
    end
elseif intersection == 1
    Rmax = 1-max(rho);
    if RadMaxDiscs==0 
        radii = Rmax.*ones(size(rho));
    else 
        radii = Rmax.*rand(size(rho));
    end
else
    R = DistanceMatrix(centers,centers) + 2*diag(1-rho);
    Rmax = min(R,[],2)/2;
    if RadMaxDiscs==0 
        radii = Rmax.*ones(size(rho));
    else 
        radii = Rmax.*rand(size(rho));
    end
end


if do_plot
    for j=1:size(centers,1)
        CC = C*radii(j)+centers(j,1)+1i*centers(j,2);
        plot(real(CC),imag(CC),'Color',[245,119,34]/255)
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
    hold on;
    if type_centers == 1 || type_centers == 2
        semilogy(degSet,LebCon(1)*degSet.^(1/2),'--','Color',[142,142,142]/255,'LineWidth',1.5)
        semilogy(degSet,LebCon(1)*degSet,'--','Color',[142,142,142]/255,'LineWidth',1.5)
    end
end
