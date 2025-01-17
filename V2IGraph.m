% V2I Graph
% Plots Secrecy Capacity (color) as a function of both Bob and Eves'
% positions. Plots image as an imagesc, but there is also a much-cooler
% surface that can be plotted by uncommenting the appropriate lines at the
% bottom. Not as useful, but way cooler looking. 
%
% Data is normalized by channel gain g, and the capacity is the 
% gaussian capacity. 


clear; 
close all;
addpath('Functions');

    % Plot Test Point A
% load('Data/GraphPwelchedData/test-point-A_graph-pwelched-data.mat');
load('Data/CorrectlyAveragedData/test-point-A_jagged-mid32.mat');
% load('Data/CorrectlyAveragedData/test-point-A_jagged.mat');

gSquared = pwelched;
[numRows, numCarriers] = size(gSquared);
SNR = 1000;

normalizedGSquared = normalizeChannelGain(gSquared);

    % With new data - find the average normalized g^2 value
sum = 0;
count = 0;
for loc = 1:numRows
    for carrier = 1:numCarriers
        sum = sum + normalizedGSquared(loc, carrier);
        count = count + 1;
    end
end

averageGSquared = sum/count;

[capacityPerLocation, capacityPerCarrier] = gaussian_capacity(normalizedGSquared, SNR);
[secrecyCapacityPerLocation, secrecyCapacityPerCarrier] = ...
    calculateSecrecyCapacity(capacityPerCarrier, capacityPerCarrier);

%% Plot Figures

[maxCapacity, indexOfMaxCapacity] = max(capacityPerLocation);

% figure()
% h = surf(secrecyCapacityPerLocation);
% set(h, 'LineStyle','none');     % Removes black lines so you can see the data
% ylabel("Bob's distance from Alice (m)");
% xlabel("Eve's distance from Alice (m)");
% zlabel("Secrecy Capacity");
% title("Secrecy Capacity");
% xticks([indexOfMaxCapacity-2000 indexOfMaxCapacity-1000 ...
%         indexOfMaxCapacity indexOfMaxCapacity+1000 ...
%         indexOfMaxCapacity+2000]);
% xticklabels([-200, -100, 0, 100, 200]); 
% yticks([indexOfMaxCapacity-2000 indexOfMaxCapacity-1000 ...
%         indexOfMaxCapacity indexOfMaxCapacity+1000 ...
%         indexOfMaxCapacity+2000]);
% yticklabels([-200, -100, 0, 100, 200]);
% xlim([indexOfMaxCapacity-2000 indexOfMaxCapacity+2000]);
% ylim([indexOfMaxCapacity-2000 indexOfMaxCapacity+2000]);
% view(45, 60);
% y = colorbar;
% ylabel(y, 'Equivocation at Eve (bits/channel use)');

figure()
imagesc(secrecyCapacityPerLocation');
ylabel("Eve's Position (m)");
xlabel("Bob's Position (m)");
title("Secrecy Capacity");
xticks([indexOfMaxCapacity-2000 indexOfMaxCapacity-1000 ...
        indexOfMaxCapacity indexOfMaxCapacity+1000 ...
        indexOfMaxCapacity+2000]);
xticklabels([-200, -100, 0, 100, 200]); 
yticks([indexOfMaxCapacity-2000 indexOfMaxCapacity-1000 ...
        indexOfMaxCapacity indexOfMaxCapacity+1000 ...
        indexOfMaxCapacity+2000]);
yticklabels([-200, -100, 0, 100, 200]);
xlim([indexOfMaxCapacity-2000 indexOfMaxCapacity+2000]);
ylim([indexOfMaxCapacity-2000 indexOfMaxCapacity+2000]);
colormap(flipud(parula))
y = colorbar;
ylabel(y, 'Secrecy Capacity (bits/channel use)');



