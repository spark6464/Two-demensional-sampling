%阿贝观点进行衍射受限成像

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
cutoff_frequency = D / 2 / lambda / do; % 截止频率

% 像面频谱网格
kethi = linspace(-1 / 2 / Lo, 1 / 2 / Lo, r) .* r;
nenta = linspace(-1 / 2 / Lo, 1 / 2 / Lo, c) .* c;
[kethi, nenta] = meshgrid(kethi, nenta); % 像的二维频谱网格

% 传递函数
H = zeros(c, r);
for n = 1:c
    for m = 1:r
        if kethi(n, m).^2 + nenta(n, m).^2 <= cutoff_frequency.^2
            H(n, m) = 1;
        end
    end
end

% 绘图
figure;
tiledlayout(1, 3, 'TileSpacing', 'compact', 'Padding', 'compact');

% 子图 1: 物面光强
nexttile;
imshow(U0, []);
title('物面光强');

% 子图 2: 传递函数
nexttile;
imshow(H, []);
title('传递函数');

% 子图 3: 像的光强分布
Gg = fftshift(fft2(U0)); % 理想像的频谱
Gi = Gg .* H; % 像的频谱
Ui = ifft2(Gi); % 逆傅立叶变换得到像的光场分布
Ii = Ui .* conj(Ui); % 像的光强分布
nexttile;
imshow(Ii, []);
title('像的光强分布');
