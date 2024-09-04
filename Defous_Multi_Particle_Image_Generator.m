
clc; clear; close all
%  This method has been patented and is only used for academic research. If used, please use the paper
%  OLR concept can refer to Xu, R., Z. Huang, W. Gong, W. Zhou, and C. Tropea. 2024. Depth from defocus technique for high number densities and non-spherical particles. Measurement, accepted for publication .
%% 图片参数
height = 1080;%图片分辨率
width = 1080;
Mag = 1.6;%放大倍率
pixelsize = 3.45;%像元大小μm
belta = 0.0118;%根据系统标定
cv = 0.01/100;%体积浓度 0.001%，根据需求设定
% max_sigma = 80;%最大直径对应的有量纲sigma
z = 20000;%测量体深度μm  （估计）开放系统可以大点或者根据光源前端面到镜头前端面的距离；管道按照尺寸，这个参数会决定在什么范围内的的颗粒被成像，后面需要转换为有量纲的σ。
picture_num = 10;%图片数量
%% 颗粒分布
mu = 100 * Mag /  pixelsize;  % 分布的均值μm
sigma = 0.1;  % 分布的标准差、分布参数
x = (5 * Mag /  pixelsize)  : (500 * Mag /  pixelsize) ;%粒径范围μm
x = x';
%给定分布类型和计算分布的概率密度函数
y = DistributionGeneration(x, mu, sigma,'LogNormal');% 'LogNormal'  'Normal'  'RR'  'Equality' 四种分布类型
figure;plot(x, y);
title('Lognormal Distribution');
xlabel('粒径(μm)');
ylabel('体积粒径分布')
%% 计算测量体体积和颗粒的总体积
Vol = (height * pixelsize / Mag) * (width * pixelsize / Mag) * z * picture_num;%测量体总体积μm^3
all_particle_Vol = Vol * cv;%颗粒总体积μm^3
%% 根据分布计算颗粒的数量
%假如粒径按照体积分布
per_bin_all_particle_Vol = all_particle_Vol * y;%得到每档粒径的总体积
per_bin_sigle_particle_Vol = (pi * x.^3)/6;%得到每档粒径的单颗粒体积μm^3
per_bin_all_particle_num = round(per_bin_all_particle_Vol ./ per_bin_sigle_particle_Vol);%得到每档粒径的颗粒总数量   四舍五入值得思考，如果有小数则向下一个粒度等级叠加，以此类推。
% per_pic_particle_num = round(per_bin_all_particle_num / picture_num);%每张图片的颗粒数和
%% 得到测量体内粒度序列，最后并进行均匀分布打乱，放在每张图像中（每张图像的颗粒数量不一定相等）
% D_all_particle = zeros(sum(per_bin_all_particle_num),1);%颗粒总数
D_all_particle = [];
[haveparticle,ind_no_zero] = find(per_bin_all_particle_num>0);
for p = haveparticle(1):haveparticle(1) + length(haveparticle)
    DD = ones(per_bin_all_particle_num(p),1) * x(p);
    D_all_particle = [D_all_particle;DD]; %将所有粒径放在一个序列
end
D_all_particle_rand = D_all_particle(randperm(length(D_all_particle)));%打乱顺序
num_rows = ceil(length(D_all_particle_rand) / picture_num);% 计算每列的数据数目，使用 ceil 确保可以容纳所有数据    如果想每张图片的颗粒数量不均匀，可以改这个地方
D_per_pic = NaN(num_rows, picture_num);% 创建一个结果矩阵，初始为 NaN 来填充空位,每个图片的颗粒对应的粒径
D_per_pic(1:length(D_all_particle_rand)) = D_all_particle_rand;% 将数据按列填充进矩阵  μm
%% 给定颗粒深度并根据σ=belta*z计算，每个颗粒的无量纲σ
depth_particle = randi([-z/2 z/2],length(D_all_particle_rand),1);%每个颗粒在测量体内随机给定一个深度，焦平面0在中间
depth_per_pic = NaN(num_rows, picture_num);% 创建一个结果矩阵，初始为 NaN 来填充空位,每个图片的颗粒对应的深度
depth_per_pic(1:length(D_all_particle_rand)) = depth_particle;% 将数据按列填充进矩阵  μm
sigma_particle_all = depth_per_pic * belta;%有量纲σ  μm
dimensionless_sigma_particle_all = sigma_particle_all ./ D_per_pic;%无量纲σ

[dimensionless_sigma_particle_all_sorted,ind_sorted] = sort(abs(dimensionless_sigma_particle_all), 1,'descend');  % 按照每列进行降序排列
D_per_pic_sorted = zeros(size(D_per_pic));  % 创建一个新的数组用于存储排序后的 B
for col = 1:size(D_per_pic, 2)
    D_per_pic_sorted(:, col) = D_per_pic(ind_sorted(:, col), col);%得到与排序后对应的粒径
end

for kk = 1:picture_num  %生成每张图像
D2 = D_per_pic_sorted(:,kk);
sigma_D2 = dimensionless_sigma_particle_all_sorted(:,kk);
D2 = D2(~any(isnan(D2), 2), :);
sigma_D2 = sigma_D2(~any(isnan(sigma_D2), 2), :);
J1 = ones(height,width);  %Set the background brightness to 0
for k = 1:length(D2)
radius2 = round(D2(k)/2);
particleX = randi([radius2 height-radius2], 1);
particleY = randi([radius2 height-radius2], 1);
%% 生成图像
J2 = zeros(height,width);   %Set the background brightness to 0
for i = 1:height
    for j = 1:height
        distance = sqrt((i - particleX)^2 + (j - particleY)^2);
        % Copy the defocused image of particle 1
        if distance <= radius2
            J2(i, j) = J1(i, j); %裁切         
        end
    end
end
sigma2 = sigma_D2(k)*radius2*2;%无量纲化为量纲，因为高斯模糊σ是有量纲的
Ij3 = imgaussfilt(J2, sigma2);%进行高斯模糊
J1 = J1-Ij3;
end
imwrite(J1,strcat('E:\Project file\DFD\70、80、90、35、40、45_Nor_Align\Synthetic image\images\',num2str(kk),'.png'));
% figure;imshow(J1);
end


