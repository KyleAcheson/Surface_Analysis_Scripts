clear all
close all

energyVar = 'energies'; % load in variables
nacmeVar = 'nacmes';
socVar = 'socs';
bondVar = 'stretch';
angVar = 'ang';
E = load('ENERGIES_ST2DF.mat', energyVar);
N = load('NACMES_ST2DF.mat', nacmeVar);
S = load('SOCS_ST2DF.mat', socVar);
B = load('ProposedGrid2.mat', bondVar);
A = load('ProposedGrid2.mat', angVar);
energies = E.(energyVar);
nacmes = N.(nacmeVar);
soc = S.(socVar);
bondlen = B.(bondVar);
ang = A.(angVar);
clear E N B A 
%% Define grid and reshaping arrays

FLAGsave = 1;
ns = 9;
nt = 10;
nsinglet = 9; % for soc array
ntriplet = 10;
nsingletE = 9;
ntripletE = 19; % for energy array
stateTitle_S = 'Singlet 4A"'; % surface plot title
stateTitle_T = 'Triplet 5A"';
file1 = 'Surface_S4App_T5App.png'; % surface plot file name
file2 = 'SOC_S4App_T5App.png'; % soc plot file name
file3 = 'idE_S4App_T5App.png'; % 1/dE plot file name

minAng = min(ang); % defining limits and constants
maxAng = max(ang);
minBL = min(bondlen);
maxBL = max(bondlen);
energyGS = min(energies(:));

energies=(energies-energyGS)*27.2114; % Convert to eV
BondLenDense = linspace(minBL,maxBL,100); % Set up query points
AngDense = linspace(minAng,maxAng,100);
[aq, bq] = meshgrid(ang, bondlen); 
[Aq,Bq] = meshgrid(AngDense, BondLenDense); % Grid of query points

[istate, jstate, rim, na, nb] = size(soc);
[na, nb, nst] = size(energies); % (no. angles, no. bond length, no. states)
%socr = squeeze(soc([1:ns], [ns+1:end], 1, :, :));

socr = squeeze(soc(:, :, 1, :, :));
soci = squeeze(soc(:, :, 2, :, :));
socrA = socr([1:ns], [ns+1:(nt*2)-1], :, :);
socrB = socr([1:ns], [nt*2:(nt*3)-1], :, :);
socrC = socr([1:ns], [nt*3:end], :, :);
socrAi = soci([1:ns], [ns+1:(nt*2)-1], :, :);
socrBi = soci([1:ns], [nt*2:(nt*3)-1], :, :);
socrCi = soci([1:ns], [nt*3:end], :, :);

SOCr = zeros(ns, nt, na, nb);
SOCi = zeros(ns, nt, na, nb);

for s=1:ns
    for t=1:nt
        for a=1:na
            for b=1:nb
                
                if socrA(s, t, a, b) ~= 0;
                    SOCr(s, t, a, b) = socrA(s, t, a, b);
                elseif socrB(s, t, a, b) ~= 0;
                    SOCr(s, t, a, b) = socrB(s, t, a, b);
                elseif socrC(s, t, na, nb) ~= 0;
                    SOCr(s, t, a, b) = socrC(s, t, a, b);
        
                end
            end
        end
    end
end
        
for s=1:ns
    for t=1:nt
        for a=1:na
            for b=1:nb
                
                if socrAi(s, t, a, b) ~= 0;
                    SOCi(s, t, a, b) = socrAi(s, t, a, b);
                elseif socrBi(s, t, a, b) ~= 0;
                    SOCi(s, t, a, b) = socrBi(s, t, a, b);
                elseif socrCi(s, t, na, nb) ~= 0;
                    SOCi(s, t, a, b) = socrCi(s, t, a, b);
        
                end
            end
        end
    end
end

%% Check crashed calculations

[aNan,bNan, ~] = ind2sub(size(energies),find(isnan(energies))); % get NaN index
energyNan = [aNan bNan];
energyNan = unique(energyNan, 'rows'); % Indexes of failed energy calculations
[~, ~, ~, dNan, eNan] = ind2sub(size(nacmes),find(isnan(nacmes)));
nacmeNan = [dNan, eNan];
nacmeNan = unique(nacmeNan, 'rows'); % Indexes of failed nacme calculations
energyNote = ['Failed Energy Calculations: ', num2str(size(energyNan, 1)), ' /', num2str(na*nb)];
nacmeNote = ['Failed NACME Calculations: ', num2str(size(nacmeNan, 1)), ' /', num2str(na*nb)];
disp(energyNote) % number of failed calculations out of total - these are splined over
disp(nacmeNote)

%% Interpolation

E1 = squeeze(energies(:,:,nsingletE));
E2 = squeeze(energies(:,:,ntripletE));
E1 = fillmissing(E1, 'spline'); % Spline missing energy points
E2 = fillmissing(E2, 'spline');

SOCij_r = squeeze(SOCr(nsinglet, ntriplet, :, :));
SOCij_i = squeeze(SOCi(nsinglet, ntriplet, :, :));

SOCij_r = fillmissing(SOCij_r, 'makima');
SOCij_i = fillmissing(SOCij_i, 'makima');

E1F = griddata(ang,bondlen,E1',Aq,Bq,'cubic'); % fitted energies
E2F = griddata(ang,bondlen,E2',Aq,Bq,'cubic');
SOCF_r = interp2(ang,bondlen,SOCij_r',Aq,Bq,'makima');
SOCF_i = interp2(ang,bondlen,SOCij_i',Aq,Bq,'makima');

dE = E2F-E1F; % Energy difference of two states
idE = 1./dE;  % 1/Energy difference - brings out singularities

%% Plotting 2 Surfaces

Fig1 = figure
subplot(1,2,1)
surf(Aq,Bq,E1F);
title(stateTitle_S)           %%%%%%%%%% EDIT
shading interp
view(120,30);
xlim([90 180])
ylim([1.2 4.0])
xlabel('$\theta_{\textrm{scs}}$ (Degrees)', 'interpreter','latex')
ylabel('$\mathit{r}_{\textrm{cs}}$ (\r{A})', 'interpreter','latex')
zlabel('Energy (eV)')
axs = gca;
axs.XAxis.MinorTick = 'on';
axs.YAxis.MinorTick = 'on';
axs.XAxis.Label.Position = [120 maxBL min(E1F(:))-5];
axs.XAxis.Label.HorizontalAlignment = 'left';
axs.YAxis.Label.Position = [maxAng 2.0 min(E1F(:))-5 ];
axs.XRuler.MinorTickValues = ang;
axs.YRuler.MinorTickValues = bondlen;
colormap(jet);
caxis([0 10])

subplot(1,2,2)
surf(Aq,Bq,E2F);
title(stateTitle_T)           %%%%%%%%%% EDIT
xlim([90 180])
ylim([1.2 4.0])
xlabel('$\theta_{\textrm{scs}}$ (Degrees)', 'interpreter','latex')
ylabel('$\mathit{r}_{\textrm{cs}}$ (\r{A})', 'interpreter','latex')
zlabel('Energy (eV)')
view(120,30);
axs = gca;
axs.XAxis.MinorTick = 'on';
axs.YAxis.MinorTick = 'on';
axs.XAxis.Label.Position = [120 maxBL min(E2F(:))-5];
axs.XAxis.Label.HorizontalAlignment = 'left';
axs.YAxis.Label.Position = [maxAng 2.0 min(E2F(:))-5];
axs.XRuler.MinorTickValues = ang;
axs.YRuler.MinorTickValues = bondlen;
colormap(jet);
caxis([0 10])
shading interp
cbp = get(subplot(1,2,2),'Position');
colorbar('Position', [cbp(1)+cbp(2)*3.5  cbp(2)*2  0.015  cbp(2)+cbp(3)+0.15])

if FLAGsave == 1
    Fig1.WindowState = 'fullscreen';
    saveas(Fig1, file1)
end


%% Plotting SOC

Fig2 = figure
subplot(1,2,1)
surf(Aq,Bq,SOCF_r);
title('SOC - Real', 'interpreter', 'latex')
hold on
view(89.999,90);
ylim([1.2 4.0]);
xlabel('$\theta_{\textrm{scs}}$ (Degrees)', 'interpreter','latex')
ylabel('$\mathit{r}_{\textrm{cs}}$ (\r{A})', 'interpreter','latex')
axs = gca;
axs.XAxis.MinorTick = 'on';
axs.YAxis.MinorTick = 'on';
axs.XRuler.MinorTickValues = ang;
axs.YRuler.MinorTickValues = bondlen;
colormap(jet);
caxis([-200 200])
shading interp

subplot(1,2,2)
surf(Aq,Bq,SOCF_i);
title('SOC - Imag', 'interpreter', 'latex')
view(89.999,90);
ylim([1.2 4.0]);
xlabel('$\theta_{\textrm{scs}}$ (Degrees)', 'interpreter','latex')
ylabel('$\mathit{r}_{\textrm{cs}}$ (\r{A})', 'interpreter','latex')
axs = gca;
axs.XAxis.MinorTick = 'on';
axs.YAxis.MinorTick = 'on';
axs.XRuler.MinorTickValues = ang;
axs.YRuler.MinorTickValues = bondlen;
colormap(jet);
caxis([-200 200])
shading interp
cbp = get(subplot(1,2,2),'Position');
colorbar('Position', [cbp(1)+cbp(2)*3.5  cbp(2)*2  0.015  cbp(2)+cbp(3)+0.15])

if FLAGsave == 1
    Fig2.WindowState = 'fullscreen';
    saveas(Fig2, file2)
end


Fig3 = figure
surf(Aq,Bq,idE);
title('1/dE', 'interpreter', 'latex')
hold on
view(89.999,90);
ylim([1.2 4.0]);
xlabel('$\theta_{\textrm{scs}}$ (Degrees)', 'interpreter','latex')
ylabel('$\mathit{r}_{\textrm{cs}}$ (\r{A})', 'interpreter','latex')
axs = gca;
axs.XAxis.MinorTick = 'on';
axs.YAxis.MinorTick = 'on';
axs.XRuler.MinorTickValues = ang;
axs.YRuler.MinorTickValues = bondlen;
colormap(jet);
caxis([0 20])
shading interp
colorbar

if FLAGsave == 1
    Fig3.WindowState = 'fullscreen';
    saveas(Fig3, file3)
end

