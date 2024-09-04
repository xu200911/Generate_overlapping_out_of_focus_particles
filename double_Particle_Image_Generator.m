
clc; clear; close all
%  This method has been patented and is only used for academic research. If used, please use the paper
%  OLR concept can refer to Xu, R., Z. Huang, W. Gong, W. Zhou, and C. Tropea. 2024. Depth from defocus technique for high number densities and non-spherical particles. Measurement, accepted for publication .
%% Pictures and particle parameters
height = 480;
width = 480;
D0 = 80;
radius = D0/2;
%% First generate out-of-focus particles far from the lens
J1 = ones(height,width);  %Set the background brightness to 0
centerX = height / 2;
centerY = height / 2;
for i = 1:height
    for j = 1:height
        % Calculate the distance from the pixel to the center of the circle
        distance = sqrt((i - centerX)^2 + (j - centerY)^2);
        % If the distance is less than the radius, set the pixel value to 1
        if distance <= radius
            J1(i, j) = 0;
            
        end
    end
end

sigma_D1 = 0.2;%The blur level of particle 1
sigma1 = sigma_D1*D0;
Ij = imgaussfilt(J1, sigma1);
% figure;imshow(J1);
% figure;imshow(Ij);

%% Second generate out-of-focus particles close from the lens
J2 = zeros(height,width);   %Set the background brightness to 0
OLR = 0.7;% Given OLR
% Determine the position of the second particle
par_A_B_distance = ((D0/2) +(D0/2))* OLR;
theta = linspace(0, 2*pi, 100);
particleX_ori = centerX + par_A_B_distance * cos(theta);
particleY_ori = centerY + par_A_B_distance * sin(theta);
rand_location_B = randi([1 length(particleX_ori)],1);
particleX = particleX_ori(rand_location_B);
particleY = particleY_ori(rand_location_B);
% particleX = height / 2;
% particleY = (height / 2) + 50;

for i = 1:height
    for j = 1:height
        distance = sqrt((i - particleX)^2 + (j - particleY)^2);
        % Copy the defocused image of particle 1
        if distance <= radius
            J2(i, j) = Ij(i, j);          
        end
    end
end
sigma_D2 = 0.1;%The blur level of particle 2
sigma2 = sigma_D2*D0;
Ij3 = imgaussfilt(J2, sigma2);
I_final = Ij-Ij3;
figure;imshow(I_final);

