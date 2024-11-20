% %二维抽样实验
clc; clear;

% 构建连续带限函数
fxy = cos(peaks(256) .* 4 + pi) + 1; % 连续带限函数
[rr, cc] = size(fxy); % 连续函数的大小

% 计算带限函数频谱
F = fftshift(fft2(fxy));

% 生成 comb 函数
combxy = zeros(rr, cc); % 初始化 comb 函数
X = 4; Y = 4; % 抽样间隔
for n = 1:Y:rr
    for m = 1:X:cc
        combxy(n, m) = 1;
    end
end

% comb 函数频谱
C = fftshift(fft2(combxy));

% 生成抽样函数
gxy = fxy .* combxy; % 抽样后的函数
Gs = fftshift(fft2(gxy)); % 抽样函数频谱

% 二维矩函数滤波器
By = round(rr / 2 / Y);
Bx = round(cc / 2 / X);
H = zeros(rr, cc);
H(round(rr/2)+1-By:round(rr/2)+1+By-1, round(cc/2)+1-Bx:round(cc/2)+1+Bx-1) = 1;

% 用二维矩函数进行滤波
Gsyp = H .* Gs;
gxyyp = X * Y .* abs(ifft2(Gsyp)); % 逆傅里叶变换还原原函数

% 绘制所有图像为子图
figure;

% 子图 1: 带限函数
subplot(3, 4, 1);
imshow(fxy, []);
title('带限函数 f(x, y)');

% 子图 2: 带限函数频谱（横向带宽）
subplot(3, 4, 2);
plot(abs(F(round(rr/2)+1, :)));
title('f(x, y) 频谱 (x 方向)');

% 子图 3: 带限函数频谱（纵向带宽）
subplot(3, 4, 3);
plot(abs(F(:, round(cc/2)+1)));
title('f(x, y) 频谱 (y 方向)');

% 子图 4: 带限函数频谱 3D
subplot(3, 4, 4);
surfl(abs(F)), shading interp, colormap("gray");
title('f(x, y) 频谱 3D');

% 子图 5: comb 函数
subplot(3, 4, 5);
imshow(combxy, []);
title('comb 函数');

% 子图 6: comb 函数频谱 3D
subplot(3, 4, 6);
surfl(abs(C)), shading interp, colormap(gray);
title('comb 函数频谱 3D');

% 子图 7: 抽样函数
subplot(3, 4, 7);
imshow(gxy, []);
title('抽样函数 g(x, y)');

% 子图 8: 抽样函数频谱 3D
subplot(3, 4, 8);
surfl(abs(Gs)), shading interp, colormap(gray);
title('g(x, y) 频谱 3D');

% 子图 9: 抽样函数频谱切片
subplot(3, 4, 9);
plot(abs(Gs(:, cc/2+1)));
title('g(x, y) 频谱切片');

% 子图 10: 二维矩函数滤波器
subplot(3, 4, 10);
imshow(H, []);
title('二维矩函数滤波器');

% 子图 11: 滤波后频谱 3D
subplot(3, 4, 11);
surfl(abs(Gsyp)), shading interp, colormap(gray);
title('滤波后频谱 3D');

% 子图 12: 还原的原函数
subplot(3, 4, 12);
imshow(gxyyp, []);
title('还原的原函数');

 