%%PSNR��Matlabʵ�� ��
function PSNR = psnr(f1, f2)
%��������ͼ��ķ�ֵ����� ����
k = 8; 
%kΪͼ���Ǳ�ʾ�ظ����ص����õĶ�����λ������λ� ��
fmax = 2.^k - 1; 
a = fmax.^2;
e = double(f1) - double(f2);
[m, n] = size(e);
b = sum(sum(e.^2)); 
PSNR = 10*log10(m*n*a/b);