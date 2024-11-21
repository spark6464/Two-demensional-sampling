%计算衍射受限单透镜系统的相⼲传递函数CTF/光学传递函数OTF，比较不同孔径尺寸成像

clc; clear;

% 读取图像并处理
U0 = imread('分辨率板USAF-1951.bmp'); % 读取图像
U0 = double(U0(:, :, 1)); % 提取第一层，转为双精度

% 获取图像尺寸
[c, r] = size(U0);

% 光学参数
lambda = 6328e-10; % 波长
k = 2 * pi / lambda; % 波数
D = 0.01; % 透镜孔径
f = 0.4; % 透镜焦距

Lo = 0.005; % 物面尺寸
do = 1.2; % 物瞳距，物到透镜距离
di = do * f / (do - f); % 透镜到观察屏距离
cutoff_frequency = D / 2 / lambda / di; % 截止频率
Li = Lo * di / do; % 像面尺寸

% 像面频谱网格
kethi = linspace(-1 / 2 / Li, 1 / 2 / Li, r) .* r;
nenta = linspace(-1 / 2 / Li, 1 / 2 / Li, c) .* c;
[kethi, nenta] = meshgrid(kethi, nenta); % 像的二维频谱网格

% 传递函数
H = zeros(c, r);
for n = 1:c
    for m = 1:r
        if kethi(n, m)^2 + nenta(n, m)^2 <= cutoff_frequency^2
            H(n, m) = 1;
        end
    end
end

% 绘图
figure;
tiledlayout(2, 3, 'TileSpacing', 'compact', 'Padding', 'compact');

% 子图 1: 相干传递函数
nexttile;
surfl(H);
shading interp;
colormap(gray);
title('相干传递函数');

% 子图 2: 算法1得到的光学传递函数OTF
nexttile;
h = fftshift(fft2(H)); % 计算相干照明下的脉冲响应
HH = abs(fftshift(fft2(h .* conj(h)))); % 计算CTF的自相关运算
OTF1 = HH ./ max(max(HH)); % 对自相关运算结果归一化
surfl(OTF1);
shading interp;
colormap(gray);
title('算法1得到的光学传递函数OTF');

% 子图 3: 算法2得到的光学传递函数OTF
nexttile;
[phi, rou] = cart2pol(kethi, nenta);
OTF2 = zeros(c, r);
for n = 1:c
    for m = 1:r
        if rou(n, m) <= 2 * cutoff_frequency
            OTF2(n, m) = 2 * (acos(rou(n, m) / 2 / cutoff_frequency) - rou(n, m) / 2 / cutoff_frequency .* sqrt(1 - (rou(n, m) / 2 / cutoff_frequency).^2)) / pi;
        end
    end
end
surfl(OTF2);
shading interp;
colormap(gray);
title('算法2得到的光学传递函数OTF');


% 子图 4: 相干照明下的像的光强分布
nexttile;
Gg = fftshift(fft2(U0)); % 理想像的频谱
Gic = Gg .* H; % 相干照明下的像的频谱
Uic = ifft2(Gic); % 相干照明下的，逆傅立叶变换得到像的光场分布
Iic = Uic .* conj(Uic); % 相干照明下的，像的光强分布
imshow(Iic, []);
title('相干照明下的像的光强分布');

% 子图 5: OTF1得到的非相干成像
nexttile;
Gii1 = Gg .* OTF1; % 非相干照明下像的频谱（OTF1做的）
Iii1 = abs(ifft2(Gii1)); % 非相干照明下的成像（OTF1做的）
imshow(Iii1, []);
title('OTF1得到的非相干成像');
colormap(gray);

% 子图 6: OTF2得到的非相干成像
nexttile;
Gii2 = Gg .* OTF2; % 非相干照明下像的频谱（OTF2做的）
Iii2 = abs(ifft2(Gii2)); % 非相干照明下的成像（OTF2做的）
imshow(Iii2, []);
title('OTF2得到的非相干成像');
colormap(gray);

