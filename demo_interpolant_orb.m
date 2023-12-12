function demo_interpolant_orb
% -------------------------------------------------------------------------
% In this demo we compute the error with two test functions of the 
% interpolant construct by integrals on disks with centers lying on orbits 
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
typeset_rho = 4; % 1 - rand, 2 - ordered rand, 3 - equispaced, 4 - cheby
type_fun = 2;
%--------------------------------------------------------------------------


h = 1;
fprintf('\n\n')
for deg = degSet
    
fprintf('...................................................')
fprintf('\n \t degree          : %5d \n',deg)

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
    plot(centers(:,1),centers(:,2),'.','Color',[147,147,147]/255);
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

% radii_orb = Rmax*rand(NOrb,1);
radii_orb = Rmax.*ones(NOrb,1);

radii = [];
for i = 1:NOrb
radii = [radii;repmat(radii_orb(i),dimPdisk(i),1)];
end

if do_plot
    for j=1:size(centers,1)
        CC = C*radii(j)+centers(j,1)+1i*centers(j,2);
        plot(real(CC),imag(CC),'Color',[245,119,34]/255)
    end
end

Vcheb = chebVandInt(deg,centers,radii);

switch type_fun 
    case 1
        f = @(x,y) exp(x).*sin(x+y);
    case 2
        f = @(x,y) 1./(25*(x.^2+y.^2)+1);
end

[XYW,~] = cub_unit_disk_sets(77);

for i = 1:length(radii)
    F(i,1) = f(XYW(:,1)*radii(i)+centers(i,1),XYW(:,2)*radii(i)+centers(i,2))'*XYW(:,3)*radii(i)^2;
end
% F = F';

a = Vcheb\F;

x = linspace(-1,1,100);
% x(1) = []; x(end) = [];
[X_old,Y_old] = meshgrid(x);
X = X_old(:).*sqrt(1-Y_old(:).^2/2); Y = Y_old(:).*sqrt(1-X_old(:).^2/2);

PolDeg = polydeg(deg); 

VxCheb = chebpolys(deg,X);
VyCheb = chebpolys(deg,Y);
for j = 1:dimP
        V(:,j) = VxCheb(:,PolDeg(j,1)+1).*VyCheb(:,PolDeg(j,2)+1);
end

P = V*a;

P = reshape(P,size(X_old));
X = reshape(X,size(X_old)); Y = reshape(Y,size(X_old));

if do_plot
    figure(4)
    mesh(X,Y,P,'FaceAlpha','0.8','edgecolor','none','facecolor',[174,208,232]/255)
    hold on;
    mesh(X,Y,f(X,Y),'FaceAlpha','0.8','edgecolor','none','facecolor',[232,198,174]/255) 
end


err(h) = norm(P(:)-f(X(:),Y(:)),'inf');

fprintf('\n \t error           : %1.4e \n',err(h))
fprintf('...................................................\n')
h = h + 1;
end

if do_plot
    figure(5)
%     semilogy(degSet,err,'-o','Color',[0,142,142]/255,'LineWidth',1.5) % blue
    semilogy(degSet,err,'-o','Color',[142,0,0]/255,'LineWidth',1.5) % red
%     semilogy(degSet,err,'-o','Color',[142,142,0]/255,'LineWidth',1.5) % mustard
    hold on;
end

end