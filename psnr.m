%%PSNR的Matlab实现 　
function PSNR = psnr(f1, f2)
%计算两幅图像的峰值信噪比 　　
k = 8; 
%k为图像是表示地个像素点所用的二进制位数，即位深。 　
fmax = 2.^k - 1; 
a = fmax.^2;
e = double(f1) - double(f2);
[m, n] = size(e);
b = sum(sum(e.^2)); 
PSNR = 10*log10(m*n*a/b);