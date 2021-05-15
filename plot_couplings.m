clear all
close all

load('nacmes_AngularG13_dR0001')
load('energies_AngularG13_dR0001')
c = 16 % coupling combination
s = 8; % min state
ds = 1;

au2ev = 27.2114;
%bondlen = 1.4:0.025:3.0;
bondlen = 180:-5:120;
gs_energy = min(energies(:,1));

energiesev = (energies-gs_energy)*au2ev;
%step = linspace(1.4,3.0,1000);
step = linspace(180,120,1000);
[n,m] = size(energiesev);
[j,k]=size(step);
splined_energies=zeros(k,m);
for i=1:m
    splined_energies(:,i)=spline(bondlen,energiesev(:,i),step)';
end

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
splined_nacme = splined_nacme + min(splined_energies(:,s));
makima_nacme = makima_nacme + min(splined_energies(:,s));

Emin = min(energiesev(:,s));
figure
hold on
%[0.8500 0.3250 0.0980]
p0 = plot(step, splined_energies(:,1), 'LineWidth', 2, 'color', [0.4940 0.1840 0.5560]);
p1 = plot(step, splined_energies(:,s), 'LineWidth', 2, 'LineStyle', '--', 'color', [0.8500 0.3250 0.0980]);
p2 = plot(step, splined_energies(:,s+ds), 'LineWidth', 2, 'LineStyle', '--', 'color', [0.9290 0.6940 0.1250]);
ylim([0, 15])
%xlabel(['Bond Length (' char(197) ')' ])
xlabel(['Bond Angle (Degrees)'])
ylabel('Energy (eV)')
%p3 = plot(step, splined_nacme(:,1), 'color', [0 0.4470 0.7410]);
p4 = plot(step, makima_nacme(:,1), 'color', [0.6350 0.0780 0.1840]);
p5 = plot(bondlen, Emin + nacme_mag(1,:), 'color', [0.4660 0.6740 0.1880]);
plot(bondlen, Emin + nacme_mag(1,:),'o', 'color', [0.4660 0.6740 0.1880]);
p6 = plot(bondlen, linspace(Emin, Emin, n), '-k');
l = legend([p0 p1 p2 p4 p5 p6], {'$\bar{X}$','3A"', '4A"', '$\sum\|\mathbf{d_{ij}}\|$ Akima', '$\sum\|\mathbf{d_{ij}}\|$ Raw', 'Shifted Zero'});
set(l,'Interpreter','latex');
s = s-5;
file = [num2str(s), num2str(s+ds), '_Apps_AngularG13_dR0001.png'];
saveas(gcf, file);


        