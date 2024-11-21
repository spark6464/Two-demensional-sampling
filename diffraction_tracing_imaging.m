%衍射追迹实现衍射受限透镜成像
clc; clear;

% 读取图像并处理
U0 = imread('分辨率板USAF-1951.bmp'); % 读取图像
U0 = double(U0(:, :, 1)); % 提取第一层，转为双精度
[c, r] = size(U0); % 获取物面采样数
lambda = 6328e-10; % 波长
k = 2 * pi / lambda; % 波数
D = 0.06; % 透镜孔径
f = 0.4; % 透镜焦距

% 物光到透镜的衍射传递过程
L0 = 0.005; % 物面尺寸
x0 = linspace(-L0/2, L0/2, r);
y0 = linspace(-L0/2, L0/2, c); % 物面坐标
[x0, y0] = meshgrid(x0, y0);

d1 = 1.2; % 物面到透镜的距离
L = r * lambda * d1 / L0; % 衍射光在透镜前表面的尺寸
p = linspace(-L/2, L/2, r);
q = linspace(-L/2, L/2, c); % 透镜前表面坐标
[p, q] = meshgrid(p, q);
F00 = exp(1i * k * d1) / (1i * lambda * d1) * exp(1i * k / 2 / d1 * (p.^2 + q.^2));
Fpq = exp(1i * k / 2 / d1 * (x0.^2 + y0.^2));

a = U0 .* Fpq;
FUpq = fft2(a); % FFT 变换
Ffpq = fftshift(FUpq); % FFT 移位
Fufpq = F00 .* Ffpq; % 透镜前表面的光场复振幅分布
I = Fufpq .* conj(Fufpq); % 光强分布

% 计算通过透镜后的光场
DD = round(D * r / L); % 孔径对应的采样数
pxy = zeros(c, r); % 生成孔径函数
for n = 1:c
    for m = 1:r
        if (n - c / 2)^2 + (m - r / 2)^2 <= (DD / 2)^2
            pxy(n, m) = 1;
        end
    end
end

Fufpqyp = Fufpq .* pxy .* exp(-1i * k .* (p.^2 + q.^2) / 2 / f); % 透过透镜后的光场

% 计算从透镜到观察面的衍射过程
d2 = d1 * f / (d1 - f) - 0.001; % 由物像公式给出相距 d2
Lyp = r * lambda * d2 / L; % 观察面的尺寸
x = linspace(-Lyp/2, Lyp/2, r);
y = linspace(-Lyp/2, Lyp/2, c); % 观察面坐标
[x, y] = meshgrid(x, y);
F0 = exp(1i * k * d2) / (1i * lambda * d2) * exp(1i * k / 2 / d2 * (x.^2 + y.^2));
F = exp(1i * k / 2 / d2 * (p.^2 + q.^2));

% 计算再现像
re_image = fft2(Fufpqyp .* F);
re_image = re_image .* F0;
if Lyp < 0 % 成虚像时倒像
    re_image = flipud(re_image); % 上下数组颠倒
    re_image = fliplr(re_image); % 左右数组颠倒
end

% 绘制 2x2 子图
figure;
tiledlayout(2, 2, 'TileSpacing', 'compact', 'Padding', 'compact');

% 子图 1: 原始物面光强
nexttile;
imshow(U0, []);
title('物面光强');

% 子图 2: 透镜前的光强分布
nexttile;
imshow(I, []);
colormap(pink);
title('透镜前光强分布');

% 子图 3: 孔径函数
nexttile;
imshow(pxy, []);
title('透镜的孔径函数');

% 子图 4: 观察面上的光强分布
nexttile;
imshow(re_image .* conj(re_image), []);
title('观察面上的光强分布');
