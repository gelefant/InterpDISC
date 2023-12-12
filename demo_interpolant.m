function demo_interpolant
% -------------------------------------------------------------------------
% In this demo we compute the error with two test functions of the 
% interpolant construct by integrals on disks with scattered centers 
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
type_centers = 2; % 0 - halton pts, 1 - Meurant Sommariva pts, 2 - Carnicer Godes pts,  3 - mapped grid, 4 - cheby mesh
type_fun = 2;
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
    R = DistanceMatrix(centers,centers)/2 + diag(1-rho);
    Rmax = min(R,[],2);
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
%     semilogy(degSet,err,'-o','Color',[142,0,0]/255,'LineWidth',1.5) % red
    semilogy(degSet,err,'-o','Color',[142,142,0]/255,'LineWidth',1.5) % mustard
    hold on;
end