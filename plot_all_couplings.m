%clear all
close all

%load('nacmes_linear30_dR0000001')
load('energies_20_180.mat')
% c = 16 % coupling combination
% s = 8; % min state
% ds = 1;

au2ev = 27.2114;
bondlen = 20:5:180;
%bondlen = stretch;
gs_energy = -833.0619;% min(energies(:,1));

energiesev = (energies-gs_energy)*au2ev;
step = linspace(20,180,1000);
[n,m] = size(energiesev);
[j,k]=size(step);
splined_energies=zeros(k,m);
for i=1:m
    splined_energies(:,i)=spline(bondlen,energiesev(:,i),step)';
end

figure
for s=1:9
hold on
if s < 6
p1 = plot(step, splined_energies(:,s), 'LineWidth', 2);
else
p1 = plot(step, splined_energies(:,s), 'LineWidth', 2, 'LineStyle', '--');
end
xlabel(['Bond Length (' char(197) ')' ])
ylabel('Energy (eV)')
ylim([0 20])
end

for c=1:16
nacme_single = nacmes(:,:,c,:);
nacme_mag = zeros(1,n);
for i=1:n
    magnitude = 0;
    for a=1:3
        magnitude = magnitude + norm(nacme_single(a,:,:,i));
    end
    nacme_mag(i) = magnitude;
end
splined_nacme=zeros(k,1);
makima_nacme=zeros(k,1);
splined_nacme(:,1)=spline(bondlen,nacme_mag(1,:),step)';
makima_nacme(:,1)=makima(bondlen,nacme_mag(1,:),step)';
p2 = plot(step, makima_nacme(:,1));
end