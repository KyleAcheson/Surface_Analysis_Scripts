clear all
close all

%load('SOC_TEST.mat')
load('energiesAng.mat')
au2ev = 27.2114;
bondlen = 120:5:180;
ns = 9;
nt = 10;
nsinglet = 9;
ntriplet = 10;
%gsenergy = min(energies(:));
gsenergy =  -833.0620;

%[istate, jstate, rim, n] = size(soc);
[n, nstate] = size(energies);
%socr = squeeze(soc([1:ns], [ns+1:end], 1, :));

energiesev = (energies-gsenergy)*au2ev;
%energiesev = fillmissing(energiesev, 'spline');
step = linspace(120,180,100);
[n,m] = size(energiesev);
[j,k]=size(step);
splined_energies=zeros(k,m);
for i=1:m
    splined_energies(:,i)=spline(bondlen,energiesev(:,i),step)';
end

% socr(isnan(socr)) = 0;
% splined_soc = zeros(nsinglet,3*ntriplet,k);
% for i=1:nsinglet
%     for j=1:ntriplet
%         splined_soc(i,j,:) = makima(bondlen,socr(i,j,:),step)';
%     end
% end
        

figure
for s=1:nsinglet
    plot(step, splined_energies(:,s), '-r')
    hold on
end
for s=1:ntriplet
    s = ns + s;
    plot(step, splined_energies(:,s), '--b')
    hold on
end

% for i=1:nsinglet
%     for j=1:ntriplet
%         soc_ij_m0 = squeeze(splined_soc(i,j,:));
%         soc_ij_m1 = squeeze(splined_soc(i,j+10,:));
%         soc_ij_mm1 = squeeze(splined_soc(i,j+20,:));
%         if  any(soc_ij_m0(:));
%             plot(step, abs(soc_ij_m0))
%             hold on
%             plot(bondlen, squeeze(abs(socr(i,j,:))), 'o')
%         elseif any(soc_ij_m1(:));
%             plot(step, abs(soc_ij_m1))
%             hold on
%             plot(bondlen, squeeze(abs(socr(i,j+10,:))), 'o')
%         elseif any(soc_ij_m1(:));
%             plot(step, abs(soc_ij_mm1))
%             hold on
%             plot(bondlen, squeeze(abs(socr(i,j+20,:))), 'o')
%         else
%             continue
%         end
%     end
% end
xlabel(['Bond Length (' char(197) ')' ])
ylabel('Energy (eV)')
ylim([ 0 10])            
     
       




