clear all
close all

load('ENERGIES_3D.mat')

cs1 = 1.4:0.1:2.0;;
cs2 = 1.4:0.1:2.0;;
ang = 140:5:180;

[na, nc1, nc2, nstate] = size(energies);

energiesGS = min(energies(:));
energies = energies(:,:,:,1);
energies = (energies-energiesGS)*27.2114;

[CS1, ANG, CS2] = meshgrid(cs1, ang, cs2);

figure
scatter3(CS1(:), ANG(:), CS2(:), [], energies(:))


%%

cs1_dense = linspace(min(cs1), max(cs1), 20);
cs2_dense = linspace(min(cs2), max(cs2), 20);
ang_dense = linspace(min(ang), max(ang), 20);

[CS1d, ANGd, CS2d] = meshgrid(cs1_dense, ang_dense, cs2_dense);

E1 = fillmissing(energies, 'spline');

E1F = interp3(cs1, ang, cs2, E1, CS1d, ANGd, CS2d, 'spline');
%E1F = interp2(cs1, ang, E1, CS1, ANG);
%surf(ANG, CS1, E1F)

figure
%slice(CS1d,ANGd,CS2d,E1F,[], 180, [])
slice(CS1d,ANGd,CS2d,E1F,[], [], [1.5, 1.7, 2.0])
%slice(CS1d,ANGd,CS2d,E1F,cs1, ang, cs2)
shading interp
colormap(jet)
caxis([0 10])
ylabel('$\theta_{\textrm{scs}}$ (Degrees)', 'interpreter','latex')
xlabel('$\mathit{r}_{\textrm{a}}$ (\r{A})', 'interpreter','latex')
zlabel('$\mathit{r}_{\textrm{b}}$ (\r{A})', 'interpreter','latex')
colorbar

figure
slice(CS1d,ANGd,CS2d,E1F,[], [180, 160, 140], [])
shading interp
colormap(jet)
caxis([0 10])
ylabel('$\theta_{\textrm{scs}}$ (Degrees)', 'interpreter','latex')
xlabel('$\mathit{r}_{\textrm{a}}$ (\r{A})', 'interpreter','latex')
zlabel('$\mathit{r}_{\textrm{b}}$ (\r{A})', 'interpreter','latex')
colorbar

lvls = 0:0.05:5;
figure
slice(CS1d,ANGd,CS2d,E1F,[1.5, 1.7, 2.0], [], [])
shading interp
colormap(jet)
caxis([0 10])
ylabel('$\theta_{\textrm{scs}}$ (Degrees)', 'interpreter','latex')
xlabel('$\mathit{r}_{\textrm{a}}$ (\r{A})', 'interpreter','latex')
zlabel('$\mathit{r}_{\textrm{b}}$ (\r{A})', 'interpreter','latex')
colorbar

figure
[FX, FY, FZ] = gradient(E1F);
slice(CS1d,ANGd,CS2d,FX,[], [], 1.5)

figure
slice(CS1d,ANGd,CS2d,FY,[], [], 1.5)

figure
slice(CS1d,ANGd,CS2d,FZ,[], [], 1.5)
