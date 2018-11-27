%% Read the recupered image:
fp = fopen('gabor_resp.bin', 'r');

height = fread(fp, 1, 'int');
width = fread(fp, 1, 'int');
gab_bin = fread(fp, height*width, 'double');
fclose(fp);

min(gab_bin/300/300)
max(gab_bin/300/300)

gab_resp = reshape(gab_bin, [height, width]) / 300/300;

%% Verify the High-Pass filter
fp = fopen('hp_filter.bin', 'r');

par_T = fread(fp, 1, 'int');
par_L = fread(fp, 1, 'double');
height = fread(fp, 1, 'int');
width = fread(fp, 1, 'int');

hpf_bin = fread(fp, 2*height*(width/2 + 1), 'double');

fclose(fp);

hpf = reshape(hpf_bin, [2, (width/2+1), height]);
hpf_real = squeeze(hpf(1, :, :))';
hpf_imag = squeeze(hpf(2, :, :))';

subplot(1, 2, 1)
imshow(hpf_real, [])
subplot(1, 2, 2)
imshow(hpf_imag, [])

%% Verify the Gabor kernels
fp = fopen('gabor_kernels.bin', 'r');
par_T = fread(fp, 1, 'int');
par_L = fread(fp, 1, 'double');
par_K = fread(fp, 1, 'int');
height = fread(fp, 1, 'int');
width = fread(fp, 1, 'int');
gk_bin = fread(fp, par_K*2*height*(width/2 + 1), 'double');
fclose(fp);

gk = reshape(gk_bin, [2, (width/2+1), height, par_K]);
gk_real = permute(squeeze(gk(1, :, :, :)), [2, 1, 3]);
gk_imag = permute(squeeze(gk(2, :, :, :)), [2, 1, 3]);

%% Read the individual responses to the Gabor kernels:
gab_resp2kernel = cell(par_K, 1);

for k=1:par_K
    fp = fopen(sprintf('resp2kernel_%03i.bin', k-1), 'r');

    height = fread(fp, 1, 'int');
    width = fread(fp, 1, 'int');
    par_K_temp = fread(fp, 1, 'int');
    k_temp = fread(fp, 1, 'int');
    
    gab_resp2kernel{k} = reshape(fread(fp, height*width, 'double'), [height, width]);
    fclose(fp);
end

%% Test the origional code:
addpath('D:\Apps\Detection_tools\Gabor_tools')
img = 1.0 - double(imread('D:\test_data\Angios_134\21_ascii.pgm')) / 255.0;
ft_img = fft2(img);

hp_filter = getGaborHPF(300, 300, 13, 2.65);
gabor_kernels = getGaborKernels(300, 300, 13, 2.65, 45);

[resp, angs, resp2kernels] = gabor(ft_img, 13, 2.65, 45);
imshow(resp{1}, [])
min(resp{1}(:))
max(resp{1}(:))

%% Verify the gabor kernels:
all_gk_diffs = zeros(par_K, 1);
for k = 1:par_K;
    gk_diff = (gabor_kernels{1}{k}(:, 1:151) - squeeze(gk_real(:, :, k))).^2;
    all_gk_diffs(k) = sum(gk_diff(:));
end
sum(all_gk_diffs)

subplot(1, 3, 1)
imshow(squeeze(gk_real(:, :, k)), [])
subplot(1, 3, 2)
imshow(gabor_kernels{1}{k}(:, 1:151), []);
subplot(1, 3, 3)
imshow(gk_diff, [])


%% Verify the High pass filter:
hpf_diff = (hp_filter{1}(:, 1:151) - hpf_real).^2;
imshow(hpf_diff, [])
sum(gk_diff(:))


%% Verify the response to each Gabor kernel
for k = 1:par_K
    diff_resp2kernels = (resp2kernels{1}{k}-gab_resp2kernel{k}').^2;
    sum(diff_resp2kernels(:))
end

subplot(2, 2, 1)
imshow(gab_resp2kernel{k}', [])
subplot(2, 2, 2)
imshow(resp2kernels{1}{k}, [])
subplot(2, 2, 3)
imshow(diff_resp2kernels, []);
subplot(2, 2, 4)
imshow(log(diff_resp2kernels), []);

%% Verify the response to the Gabor filter:
dif_resp = (resp{1} - gab_resp').^2;
sum(dif_resp(:))
subplot(2, 2, 1);
imshow(gab_resp', []);
subplot(2, 2, 2);
imshow(resp{1}, []);
subplot(2, 2, 3);
imshow(dif_resp, []);
subplot(2, 2, 4);
imshow(log(dif_resp), []);
