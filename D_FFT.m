%菲涅耳衍射积分D-FFT算法

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


kethi=linspace(-1./2./L0,1./2./L0,c).*c;%频域坐标
nenta=linspace(-1./2./L0,1./2./L0,r).*r;%频域坐标
[kethi,nenta]=meshgrid(kethi,nenta);%频域坐标网格

H=exp(1i*k*d.*(1-lambda.*lambda.*(kethi.*kethi+nenta.*nenta)./2));%传递函数H
fa=fftshift(fft2(a));%衍射面光场的傅立叶变换
Fuf=fa.*H;%光场频谱与传递函数相乘
U=ifft2(Fuf);%逆傅立叶变换得到观察屏光场分布
I=U.*conj(U);%观察屏光强分布
figure,imshow(I,[0,max(max(I))]),colormap(gray),title('实际计算得到的衍射光强分布')












