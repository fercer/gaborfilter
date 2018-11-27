fp = fopen('fouriertransform.bin', 'r');

height = fread(fp, 1, 'int');
width = fread(fp, 1, 'int');

ft_real_bin = fread(fp, height*(width/2+1), 'double');
ft_imag_bin = fread(fp, height*(width/2+1), 'double');

fclose(fp);

ft_real = reshape(ft_real_bin, [width/2+1, height]);
ft_imag = reshape(ft_imag_bin, [width/2+1, height]);
ft_img = complex(ft_real, ft_imag);

imshow(log(real(ft_img).^2 + imag(ft_img).^2), [])

img = ifft2(ft_img);
size(img)
imshow(img, [])

%% Read the recupered image:
fp = fopen('D:/Documentos/test/ift.bin', 'r');

height = fread(fp, 1, 'int');
width = fread(fp, 1, 'int');
ifimg_bin = fread(fp, height*width, 'double');
fclose(fp);

ifimg_linux = reshape(ifimg_bin / 300 / 300, [height, width]);
imshow(ifimg_linux, []);

max(ifimg_linux(:))
min(ifimg_linux(:))

%% Read the recupered image:
fp = fopen('ift.bin', 'r');

height = fread(fp, 1, 'int');
width = fread(fp, 1, 'int');
ifimg_bin = fread(fp, height*width, 'double');
fclose(fp);

ifimg = reshape(ifimg_bin / 300 / 300, [height, width]);
imshow(ifimg, []);

max(ifimg(:))
min(ifimg(:))

dif_rec = (ifimg_linux - ifimg).^2;
sum(dif_rec(:))

%% Read the original image:
img_or = imread('D:/test_data/Angios_134/21_ascii.pgm');
max(img_or(:))
min(img_or(:))
img_or = 1.0 - double(img_or) / 255.0;

dif_img = (img_or - ifimg).^2;
sum(dif_img(:))

ft_img_or = fft2(img_or');
ift_img_or = ifft2(ft_img_or);
max(ift_img_or(:))
min(ift_img_or(:))

%% Verify the results from matlab and c
ift_diff = ft_img_or(1:151, :) - ft_img;
sum(sum(ift_diff.^2))