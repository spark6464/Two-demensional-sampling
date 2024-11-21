% %二维光场分析-球面波

clc; clear;

% 参数设置
lambda = 6328e-10; % 波长
k = 2 * pi / lambda; % 波数
x0 = 0.001; % 光源 x 坐标 (m)
y0 = 0.001; % 光源 y 坐标 (m)
z = 0.05; % 光源 z 坐标 (m)
L = 0.005; % 观察面尺寸 (m)
x = linspace(-L/2, L/2, 512);
y = x; % 构建 x, y 坐标
[x, y] = meshgrid(x, y); % 构建二维坐标网格

% 发散球面波
U1 = exp(1i * k * z) .* exp(1i * k .* ((x - x0).^2 + (y - y0).^2) / 2 / z);
ph1 = k .* ((x - x0).^2 + (y - y0).^2) / 2 / z; % 实际相位
phyp1 = angle(U1); % 包裹相位

% 会聚球面波
U2 = exp(-1i * k * z) .* exp(-1i * k .* ((x - x0).^2 + (y - y0).^2) / 2 / z);
ph2 = -k .* ((x - x0).^2 + (y - y0).^2) / 2 / z; % 实际相位
phyp2 = angle(U2); % 包裹相位

% 干涉图
diff1 = U1 + 1; % 发散球面波与平面光干涉
I1 = diff1 .* conj(diff1); % 光强
diff2 = U2 + 1; % 会聚球面波与平面光干涉
I2 = diff2 .* conj(diff2); % 光强

% 使用 tiledlayout 管理布局
tiledlayout(4, 2, 'TileSpacing', 'compact', 'Padding', 'compact');

% 子图 1: 发散球面波实际相位
nexttile;
surfl(ph1), shading interp, colormap(gray);
title('发散球面波实际相位');

% 子图 2: 会聚球面波实际相位
nexttile;
surfl(ph2), shading interp, colormap(gray);
title('会聚球面波实际相位');

% 子图 3: 发散球面波包裹相位
nexttile;
imshow(phyp1, []);
title('发散球面波包裹相位');

% 子图 4: 会聚球面波包裹相位
nexttile;
imshow(phyp2, []);
title('会聚球面波包裹相位');

% 子图 5: 发散球面波干涉光强
nexttile;
imshow(I1, [0, max(max(I1))]);
title('发散球面波干涉光强');

% 子图 6: 会聚球面波干涉光强
nexttile;
imshow(I2, [0, max(max(I2))]);
title('会聚球面波干涉光强');

% 子图 7: 相位剖线对比，占据最后一整行
nexttile([1 2]); % 占据两列
plot(ph2(257, :), '--', 'DisplayName', '实际相位');
hold on;
plot(phyp2(257, :), 'r', 'DisplayName', '包裹相位');
hold off;
legend();
title('相位剖线对比');
xlabel('像素位置');
ylabel('相位值');
grid on;
