%% Plots the number of good carriers and unique good carriers for harrison
%% at each SNR level
%% Created by Dakota Flanary

clear;
close all;
load('Data/CorrectlyAveragedData/test-point-A_jagged-mid32.mat'); % loads in the signal data
alpha = sqrt(pwelched);
load('Data/CorrectlyAveragedData/test-point-F_average41-mid32.mat'); % loads in the signal data
foxtrot = pwelched;
load('Data/CorrectlyAveragedData/test-point-G_average41-mid32.mat'); % loads in the signal data
golf = pwelched;
load('Data/CorrectlyAveragedData/test-point-H_average41-mid32.mat'); % loads in the signal data
hotel = pwelched;
load('Data/CorrectlyAveragedData/test-point-I_average41-mid32.mat'); % loads in the signal data
india = pwelched;
load('Data/CorrectlyAveragedData/test-point-J_average41-mid32.mat'); % loads in the signal data
juliet = pwelched;
clear pwelched;

% change to true if you want to display figures
iWantFigures = true;

alpha_max = max(max(alpha));
alpha = (alpha/alpha_max).^2; %normalize the alpha (v2i) run

% find the max from each run
fox_max = max(max(foxtrot));
golf_max = max(max(golf));
hotel_max = max(max(hotel));
india_max = max(max(india));
juliet_max = max(max(juliet));

%normalize all of the v2v runs based on the overall max
normalize_max = sqrt(max([fox_max,golf_max,hotel_max,india_max,juliet_max]));
foxtrot = (sqrt(foxtrot)/normalize_max).^2;
golf = (sqrt(golf)/normalize_max).^2;
hotel = (sqrt(hotel)/normalize_max).^2;
india = (sqrt(india)/normalize_max).^2;
juliet = (sqrt(juliet)/normalize_max).^2; 

addpath(genpath('Functions')) % add the functions to this path

[v2i_num_locations,v2i_num_carriers] = size(alpha); % get the size of the v2i arrays
[v2v_num_locations,v2v_num_carriers] = size(foxtrot); % get the size of the v2v arrays


num_loops = 50; % number of snr values to test
carriers_alpha = zeros(v2i_num_locations,num_loops); % capacity from the alpha run
carriers_foxtrot = zeros(v2v_num_locations,num_loops); % capacity from the foxtrot run
carriers_golf = zeros(v2v_num_locations,num_loops); % capacity from the golf run
carriers_hotel = zeros(v2v_num_locations,num_loops); % capacity from the hotel run
carriers_india = zeros(v2v_num_locations,num_loops); % capacity from the india run
carriers_juliet = zeros(v2v_num_locations,num_loops); % capacity from the juliet run

snr = logspace(-0,5,num_loops); % snr values to test(in linear)
% snr = 6.158482110660267e+02;
snrdB = 10*log10(snr);
db_threshold = 0; % the threshold to test in db
threshold = 10^(db_threshold/10); % the threshold to test converted to linear
% threshold = 0.057634511544837;
% calculate the capacities at different SNR values for each v2v run
for index = 1:num_loops
    % Find the number of carriers(capacity) of each run at every location
    carriers_foxtrot(:,index) = find_num_carriers(foxtrot,snr(index),threshold);
    carriers_golf(:,index) = find_num_carriers(golf,snr(index),threshold);
    carriers_hotel(:,index) = find_num_carriers(hotel,snr(index),threshold);
    carriers_india(:,index) = find_num_carriers(india,snr(index),threshold);
    carriers_juliet(:,index) = find_num_carriers(juliet,snr(index),threshold);
    
    % code for displaying figures of the capacities
    if iWantFigures && index == 16
        figure()
        hold on
        plot(carriers_foxtrot(:,index));
        plot(carriers_golf(:,index));
        plot(carriers_hotel(:,index));
        plot(carriers_india(:,index));
        plot(carriers_juliet(:,index));
        hold off
    end
end
num_loops_carriers_per_location = 50;

nmax = 41;
for n = 1:nmax
   b_k(n) = 1/nmax; 
end

equiv_foxtrot = 16 - carriers_foxtrot(:,16)/2;
equiv_golf = 16 - carriers_golf(:,16)/2;
equiv_hotel = 16 - carriers_hotel(:,16)/2;
equiv_india = 16 - carriers_india(:,16)/2;
equiv_juliet = 16 - carriers_juliet(:,16)/2;

car_foxtrot = zeros(2330,1);
car_golf = zeros(2330,1);
carriers_hotel = zeros(2330,1);
carriers_india = zeros(2330,1);
carriers_juliet = zeros(2330,1);

figure();
hold on;
grid on;
plot(equiv_foxtrot);
plot(equiv_golf);
plot(equiv_hotel);
plot(equiv_india);
plot(equiv_juliet);
legend('bob', 'eve1', 'eve2', 'eve3', 'eve4');
xlabel("Alice's Location");
ylabel('Equivocation in bits');
hold off

v2i_snr = logspace(2,3,num_loops_carriers_per_location);

% calculate the capacity for the v2i run at different SNR values 
for index = 1:num_loops_carriers_per_location
    carriers_alpha(:,index) = find_num_carriers(alpha,v2i_snr(index),threshold);
    
    % code for displaying a figure for each SNR value
%     if (iWantFigures)
%         figure()
%         hold on
%         plot(cap_alpha(:,index));
%         hold off
%     end
end 

carriers_foxtrot = carriers_foxtrot/2;
carriers_golf = carriers_golf/2;
carriers_hotel = carriers_hotel/2;
carriers_india = carriers_india/2;
carriers_juliet = carriers_juliet/2;

snr = snr(16);

% figure()
% hold on
% grid on
% plot(carriers_foxtrot)
% plot(carriers_golf)
% plot(carriers_hotel)


save('Data/dataForUncodedCase.mat', ...
    'equiv_foxtrot', 'equiv_golf', 'equiv_hotel', ... 
    'equiv_india', 'equiv_juliet', 'snr');

%%



% % plot the capacities and secrecy capacities
% snr = 10*log10(snr);
% figure()
% hold on;
% plot(snr,cap_har,'DisplayName','Bob');
% plot(snr,eve,'DisplayName','Eve');
% plot(snr,sec_har, 'DisplayName', 'Secrecy Capacity');
% xlabel('SNR (dB)');
% ylabel('Bits per channel use');
% legend;
% hold off;
