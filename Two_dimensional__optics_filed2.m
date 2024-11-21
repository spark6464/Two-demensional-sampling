% %二维光场分析-平面波

clc; clear;

% 参数设置
lambda = 6328e-10; % 波长
k = 2 * pi / lambda; % 波数
alpha = pi / 2.005; % 光与 x 轴夹角
beita = pi / 2.005; % 光与 y 轴夹角
L = 0.004; % 观察面尺寸 (m)
x = linspace(-L/2, L/2, 512);
y = x; % 构建 x, y 坐标
[x, y] = meshgrid(x, y); % 构建二维坐标网格

% 入射平面波
U = exp(1i * k * (x .* cos(alpha) + y .* cos(beita))); % 入射平面波
ph = k .* (x .* cos(alpha) + y .* cos(beita)); % 实际相位
phyp = angle(U); % 包裹相位

% 干涉图
diff = U + 1; % 入射平面光与垂直照射平面光干涉
I = diff .* conj(diff); % 光强

% 频谱
UFuv = fftshift(fft2(U)); % 入射光场频谱
IFuv = fftshift(fft2(I)); % 干涉条纹频谱

% 使用 tiledlayout 管理布局
tiledlayout(4, 2, 'TileSpacing', 'compact', 'Padding', 'compact');

% 子图 1: 实际相位
nexttile;
surfl(ph), shading interp, colormap(gray);
title('实际相位');
xlabel('X 方向');
ylabel('Y 方向');

% 子图 2: 包裹相位
nexttile;
imshow(phyp, []);
title('包裹相位');

% 子图 3: 实际相位与包裹相位剖线
nexttile; 
plot(ph(257, :), '--', 'DisplayName', '实际相位');
hold on;
plot(phyp(257, :), 'r', 'DisplayName', '包裹相位');
hold off;
legend();
title('相位剖线对比');
xlabel('像素位置');
ylabel('相位值');
grid on;

% 子图 4: 干涉光强分布
nexttile;
imshow(I, []);
title('干涉光强分布');

% 子图 5: 入射平面光波频谱
nexttile;
imshow(abs(UFuv), [0, max(max(abs(UFuv))) / 50]);
title('平面光波频谱');

% 子图 6: 干涉光强频谱
nexttile;
imshow(abs(IFuv), [0, max(max(abs(IFuv))) / 50]);
title('干涉光强频谱');
