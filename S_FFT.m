%菲涅耳衍射积分S-FFT算法

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
L=r*lambda*d/L0;%观察屏的尺寸
x=linspace(-L/2,L/2,c);%观察屏x轴坐标
y=linspace(-L/2,L/2,r);%观察屏y轴坐标
[x,y]=meshgrid(x,y);%观察屏二维网格

%计算衍射积分
F0=exp(1i*k*d)/(1i*lambda*d)*exp(1i*k/2/d*(x.^2+y.^2));%赋值$exp(jkd)/(j\lambda d)exp[jk(x^2+y^2)/2d]$ 
F=exp(1i*k/2/d*(x0.^2+y0.^2));%赋值$exp[jk(x_{0}^2+y_{0}^2)/2d]$ 
a=a.*F;%赋值$U_{0}(x_{0},y_{0})exp[jk(x_{0}^2+y_{0}^2)/2d]$ 
Ff=fftshift(fft2(a));%对$U_{0}(x_{0},y_{0})exp[jk(x_{0}^2+y_{0}^2)/2d]$ 傅立叶变换
Fuf=F0.*Ff;%观察屏上的光场分布U(x,y)
I=Fuf.*conj(Fuf);
figure,imshow(I,[0,max(max(I))]),colormap(gray),title('实际计算得到的衍射光强分布')




