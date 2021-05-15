%% Setting up data 
clear all
close all

energyVar = 'energies'; % load in variables
nacmeVar = 'nacmes';
bondVar = 'stretch';
angVar = 'ang';
E = load('energies_ST2DF.mat', energyVar);
N = load('nacmes_ST2DF.mat', nacmeVar);
B = load('ProposedGrid2.mat', bondVar);
A = load('ProposedGrid2.mat', angVar);
energies = E.(energyVar);
nacmes = N.(nacmeVar);
bondlen = B.(bondVar);
ang = A.(angVar);
clear E N B A 
%energies = permute(energies, [2 3 1]);
%energies(energies==0) = NaN;

FLAGsave = 0; % 1 if saving plots
FLAGsymm = 0; % 1 for A" 
filename1 = '_Aps_PES.png';
filename2 = '_Aps_Coupling2D.png';
filename3 = '_Apps_VectorField.png';
c = 6; % coupling combination
s = 3; % state 1
ds = 1; % increase in state

[na, nb, ns] = size(energies); % (no. angles, no. bond length, no. states)
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

%% Interpolation - Energy, NACME magnitudes, 1/dE etc.

nacme2S = squeeze(nacmes(:,:,c,:,:)); % nacme between two states
nacmeMagnitude = zeros(na,nb); 
for i=1:na                   % Sum over atoms of vector magnitudes
    for j=1:nb
        magnitude = 0;
        for a=1:3
            magnitude = magnitude + norm(nacme2S(a,:,i,j));
        end
        nacmeMagnitude(i,j) = magnitude;
    end
end

E1 = squeeze(energies(:,:,s));
E2 = squeeze(energies(:,:,s+ds));
E1 = fillmissing(E1, 'spline'); % Spline missing energy points
E2 = fillmissing(E2, 'spline');
nacmeMagnitude = fillmissing(nacmeMagnitude, 'makima'); % Akima pchi missing nacme magnitudes
NAC = interp2(ang,bondlen,nacmeMagnitude',Aq,Bq,'makima'); % fitted NACME for 2 states
E1F = griddata(ang,bondlen,E1',Aq,Bq,'cubic'); % fitted energies
E2F = griddata(ang,bondlen,E2',Aq,Bq,'cubic');

dE = E2F-E1F; % Energy difference of two states
idE = 1./dE;  % 1/Energy difference - brings out singularities

%% plotting 2 state surfaces

Fig1 = figure
subplot(1,2,1)
surf(Aq,Bq,E1F);
title('2A''')           %%%%%%%%%% EDIT
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
title('4A''')           %%%%%%%%%% EDIT
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
    if FLAGsymm == 1
        s = s-5;
    end
    Fig1.WindowState = 'fullscreen';
    file1 = [num2str(s), num2str(s+ds), filename1];
    saveas(Fig1, file1)
end

%% plotting 1/dE and NACME magnitudes
Fig2 = figure
subplot(1,2,1)
surf(Aq,Bq,idE);
title('$\frac{1}{\Delta E}$', 'interpreter', 'latex')
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

subplot(1,2,2)
surf(Aq,Bq,NAC);
title('$\sum_{n}^{n=3}\|\mathbf{d}_{ij}\|$', 'interpreter', 'latex')
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
caxis([0 10000])
shading interp
cbp = get(subplot(1,2,2),'Position');
colorbar('Position', [cbp(1)+cbp(2)*3.5  cbp(2)*2  0.015  cbp(2)+cbp(3)+0.15])

if FLAGsave == 1
    Fig2.WindowState = 'fullscreen';
    file2 = [num2str(s), num2str(s+ds), filename2];
    saveas(Fig2, file2)
end

%% Plotting NAC vector fields

colors = [1 0 0; 0 1 0; 1 0 1];
viewAng = [120 25; 30 25];
AxesLoc = [ 0.05 0.65 0.4 0.25; 0.5 0.65 0.4 0.25; 0.05 0.35 0.4 0.25; 0.5 0.35 0.4 0.25; 0.05 0.05 0.4 0.25; 0.5 0.05 0.4 0.25];
count = 0;
Fig3 = figure;

for i=1:3
    for j=1:2
        count = count + 1;
    ax(count) = axes('Position', AxesLoc(count,:));
    contourf(ax(count), Aq, Bq, NAC, 'LevelStep', 3.0);
    
xlim([90 180])
ylim([1.2 3.0])
set(gca,'ztick',[])
set(gca,'zticklabel',[])
xlabel('$\theta_{\textrm{scs}}$ (Degrees)', 'interpreter','latex')
ylabel('$\mathit{r}_{\textrm{cs}}$ (\r{A})', 'interpreter','latex')
axs = gca;
axs.XAxis.MinorTick = 'on';
axs.YAxis.MinorTick = 'on';
axs.XRuler.MinorTickValues = ang;
axs.YRuler.MinorTickValues = bondlen;
% axs.XAxis.Label.HorizontalAlignment = 'center';
% axs.XAxis.Label.Position = [138.8109913379758,3.206670279551101,-0.701754562425361];
axs.YAxis.Label.HorizontalAlignment = 'center';
axs.YAxis.Label.Position = [185.4129137218747,1.993493869611452,-0.707822651533152];
colormap(jet);
shading interp
hold on
caxis([0 20])
qp(count) = quiver3(ax(count), aq', bq', zeros(21,30), squeeze(nacme2S(i,1,:,:)), squeeze(nacme2S(i,2,:,:))...
    ,squeeze(nacme2S(i,3,:,:)));
qp(count).LineWidth = 1.0;
set(qp(count),'AutoScale','off', 'AutoScaleFactor', 0)
qp(count).Head.LineStyle = 'solid';
qp(count).MaxHeadSize = 0.005;
qp(count).Color = colors(i,:);
%view(viewAng(j,:))
    end
end
colorbar('Position', [0.95 0.1 0.01 0.7])
annotation(Fig3,'textbox',...
    [0.860027777777778 0.892815076560659 0.0774722222222218 0.0895170789163685],...
    'String',' NAC Vectors',...
    'Interpreter','latex',...
    'FontWeight','bold',...
    'FontSize',16,...
    'FitBoxToText','off');
annotation(Fig3,'textbox',...
    [0.864194444444445 0.937573616018845 0.058722222222222 0.0200235571260284],...
    'Color',[1 0 0],...
    'String','$\rightarrow$ Carbon',...
    'Interpreter','latex',...
    'FontSize',14,...
    'FitBoxToText','off',...
    'EdgeColor','none');
annotation(Fig3,'textbox',...
    [0.864194444444445 0.919905771495877 0.0663611111111109 0.0200235571260284],...
    'Color',[0 1 0],...
    'String','$\rightarrow$ Sulphur 1',...
    'Interpreter','latex',...
    'FontSize',14,...
    'FitBoxToText','off',...
    'EdgeColor','none');
annotation(Fig3,'textbox',...
    [0.864194444444445 0.902237926972908 0.0663611111111109 0.0200235571260284],...
    'Color',[1 0 1],...
    'String','$\rightarrow$ Sulphur 2',...
    'Interpreter','latex',...
    'FontSize',14,...
    'FitBoxToText','off',...
    'EdgeColor','none');

if FLAGsave == 1
    if FLAGsymm == 1
        s = s -5;
    end
    Fig3.WindowState = 'fullscreen';
    file3 = [num2str(s), num2str(s+ds), filename3];
    saveas(Fig3, file3)
end