
clc; clear; close all
%  This method has been patented and is only used for academic research. If used, please use the paper
%  OLR concept can refer to Xu, R., Z. Huang, W. Gong, W. Zhou, and C. Tropea. 2024. Depth from defocus technique for high number densities and non-spherical particles. Measurement, accepted for publication .
%% ͼƬ����
height = 1080;%ͼƬ�ֱ���
width = 1080;
Mag = 1.6;%�Ŵ���
pixelsize = 3.45;%��Ԫ��С��m
belta = 0.0118;%����ϵͳ�궨
cv = 0.01/100;%���Ũ�� 0.001%�����������趨
% max_sigma = 80;%���ֱ����Ӧ��������sigma
z = 20000;%��������Ȧ�m  �����ƣ�����ϵͳ���Դ����߸��ݹ�Դǰ���浽��ͷǰ����ľ��룻�ܵ����ճߴ磬��������������ʲô��Χ�ڵĵĿ��������񣬺�����Ҫת��Ϊ�����ٵĦҡ�
picture_num = 10;%ͼƬ����
%% �����ֲ�
mu = 100 * Mag /  pixelsize;  % �ֲ��ľ�ֵ��m
sigma = 0.1;  % �ֲ��ı�׼��ֲ�����
x = (5 * Mag /  pixelsize)  : (500 * Mag /  pixelsize) ;%������Χ��m
x = x';
%�����ֲ����ͺͼ���ֲ��ĸ����ܶȺ���
y = DistributionGeneration(x, mu, sigma,'LogNormal');% 'LogNormal'  'Normal'  'RR'  'Equality' ���ֲַ�����
figure;plot(x, y);
title('Lognormal Distribution');
xlabel('����(��m)');
ylabel('��������ֲ�')
%% �������������Ϳ����������
Vol = (height * pixelsize / Mag) * (width * pixelsize / Mag) * z * picture_num;%�������������m^3
all_particle_Vol = Vol * cv;%�����������m^3
%% ���ݷֲ��������������
%����������������ֲ�
per_bin_all_particle_Vol = all_particle_Vol * y;%�õ�ÿ�������������
per_bin_sigle_particle_Vol = (pi * x.^3)/6;%�õ�ÿ�������ĵ����������m^3
per_bin_all_particle_num = round(per_bin_all_particle_Vol ./ per_bin_sigle_particle_Vol);%�õ�ÿ�������Ŀ���������   ��������ֵ��˼���������С��������һ�����ȵȼ����ӣ��Դ����ơ�
% per_pic_particle_num = round(per_bin_all_particle_num / picture_num);%ÿ��ͼƬ�Ŀ�������
%% �õ����������������У���󲢽��о��ȷֲ����ң�����ÿ��ͼ���У�ÿ��ͼ��Ŀ���������һ����ȣ�
% D_all_particle = zeros(sum(per_bin_all_particle_num),1);%��������
D_all_particle = [];
[haveparticle,ind_no_zero] = find(per_bin_all_particle_num>0);
for p = haveparticle(1):haveparticle(1) + length(haveparticle)
    DD = ones(per_bin_all_particle_num(p),1) * x(p);
    D_all_particle = [D_all_particle;DD]; %��������������һ������
end
D_all_particle_rand = D_all_particle(randperm(length(D_all_particle)));%����˳��
num_rows = ceil(length(D_all_particle_rand) / picture_num);% ����ÿ�е�������Ŀ��ʹ�� ceil ȷ������������������    �����ÿ��ͼƬ�Ŀ������������ȣ����Ը�����ط�
D_per_pic = NaN(num_rows, picture_num);% ����һ��������󣬳�ʼΪ NaN ������λ,ÿ��ͼƬ�Ŀ�����Ӧ������
D_per_pic(1:length(D_all_particle_rand)) = D_all_particle_rand;% �����ݰ�����������  ��m
%% ����������Ȳ����ݦ�=belta*z���㣬ÿ�������������٦�
depth_particle = randi([-z/2 z/2],length(D_all_particle_rand),1);%ÿ�������ڲ��������������һ����ȣ���ƽ��0���м�
depth_per_pic = NaN(num_rows, picture_num);% ����һ��������󣬳�ʼΪ NaN ������λ,ÿ��ͼƬ�Ŀ�����Ӧ�����
depth_per_pic(1:length(D_all_particle_rand)) = depth_particle;% �����ݰ�����������  ��m
sigma_particle_all = depth_per_pic * belta;%�����٦�  ��m
dimensionless_sigma_particle_all = sigma_particle_all ./ D_per_pic;%�����٦�

[dimensionless_sigma_particle_all_sorted,ind_sorted] = sort(abs(dimensionless_sigma_particle_all), 1,'descend');  % ����ÿ�н��н�������
D_per_pic_sorted = zeros(size(D_per_pic));  % ����һ���µ��������ڴ洢������ B
for col = 1:size(D_per_pic, 2)
    D_per_pic_sorted(:, col) = D_per_pic(ind_sorted(:, col), col);%�õ���������Ӧ������
end

for kk = 1:picture_num  %����ÿ��ͼ��
D2 = D_per_pic_sorted(:,kk);
sigma_D2 = dimensionless_sigma_particle_all_sorted(:,kk);
D2 = D2(~any(isnan(D2), 2), :);
sigma_D2 = sigma_D2(~any(isnan(sigma_D2), 2), :);
J1 = ones(height,width);  %Set the background brightness to 0
for k = 1:length(D2)
radius2 = round(D2(k)/2);
particleX = randi([radius2 height-radius2], 1);
particleY = randi([radius2 height-radius2], 1);
%% ����ͼ��
J2 = zeros(height,width);   %Set the background brightness to 0
for i = 1:height
    for j = 1:height
        distance = sqrt((i - particleX)^2 + (j - particleY)^2);
        % Copy the defocused image of particle 1
        if distance <= radius2
            J2(i, j) = J1(i, j); %����         
        end
    end
end
sigma2 = sigma_D2(k)*radius2*2;%�����ٻ�Ϊ���٣���Ϊ��˹ģ�����������ٵ�
Ij3 = imgaussfilt(J2, sigma2);%���и�˹ģ��
J1 = J1-Ij3;
end
imwrite(J1,strcat('E:\Project file\DFD\70��80��90��35��40��45_Nor_Align\Synthetic image\images\',num2str(kk),'.png'));
% figure;imshow(J1);
end


