%菲涅耳衍射积分T-FFT算法

clc;clear

r=512;c=r;%衍射面抽样数
a=zeros(r,c);%预设一个衍射面
a(r/2-r/4:r/2+r/4,c/2-c/4:c/2+c/4)=1;%设置衍射孔
lambda=6328e-10;%波长
k=2*pi/lambda;%波数
L0=5e-3;%设置衍射面尺寸/m
d=0.1;%设置观察屏到衍射面的距离/m
x0=linspace(-L0/2,L0/2,c);%衍射面x轴坐标
y0=linspace(-L0/2,L0/2,r);%衍射面y轴坐标
[x0,y0]=meshgrid(x0,y0);%衍射面二维网格

F0=exp(1i*k*d)/(1i*lambda*d);%赋值$exp(jkd)/(j\lambda d)$ 
F1=exp(1i*k/2/d*(x0.^2+y0.^2));%赋值$exp[jk(x_{0}^2+y_{0}^2)/2d]$

fa=fft2(a);%光场$U_{0}(x_{0},y_{0})$傅立叶变换
fF1=fft2(F1);%$exp[jk(x_{0}^2+y_{0}^2)/2d]$傅立叶变换
Fuf=fa.*fF1;%频谱相乘
U=F0.*fftshift(ifft2(Fuf));%逆傅立叶变换得到观察屏光场分布$U(x,y)$
I=U.*conj(U);%计算观察屏光强分布
figure,imshow(I,[0,max(max(I))]),colormap(gray),title('实际计算得到的衍射光强分布')

