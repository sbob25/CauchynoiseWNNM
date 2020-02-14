% choose image by number
% All images are in the folder "test_images".
img_num = 3;

% solve WNNM model of the Cauchy noise denoising by ADMM algorithm.

if img_num == 1; clear; imname = "barbara"; img = double(imread('barbara.png')); img = imresize(img, 0.5);
elseif img_num == 2; clear; imname = "boat"; img = double(imread('boat.png')); img = imresize(img, 0.5);
elseif img_num == 3; clear; imname = "cameraman"; img = double(imread('cameraman.tif'));
elseif img_num == 4; clear; img_name = "couple"; img = double(imread('couple.png')); img = imresize(img, 0.5);
elseif img_num == 5; clear; imname = "house"; img = double(rgb2gray(imread('house.tiff'))); img = imresize(img, 0.5);
elseif img_num == 6; clear; imname = "man"; img = double(imread('man.tiff')); img = imresize(img, 0.25);
elseif img_num == 7; clear; imname = "mandrill"; img = double(imread('mandrill.png')); img = imresize(img, 0.5);
elseif img_num == 8; clear; imname = "peppers"; img = double(imread('peppers256.png'));
elseif img_num == 9; clear; imname = "plane"; img = double(imread('jetplane.png')); img = imresize(img, 0.5);
elseif img_num == 10; clear; imname = "synthetic"; img = double(rgb2gray(imread('synthetic1.bmp'))); img = imresize(img, 1.01);
end

% Generate Cauchy noisy image.
rng('default');
eta1 = randn(size(img));
eta2 = randn(size(img));
gamma = 5;
noise = gamma.*eta1./eta2;
Y = img + noise;                                % noisy image

fprintf('noisy image PSNR, SSIM : %.2f, %.4f\n', psnr(uint8(Y), uint8(img)), ssim(uint8(Y), uint8(img)));

maxiter = 100;
lambda = 1;
beta = 2*lambda/gamma^2 + eps;
energy = zeros(1,maxiter+1);
W = zeros(size(img));
E = zeros(size(img));
medY = medfilt2(Y);
X = max(min(medY(:)),min(Y,max(medY(:))));      % set initial image

% sigma = 2.7 for gamma = 5
% sigma = 4.0 for gamma = 10
sigma = 2.7;
Par = Params(sigma);

[NL_mat, refIdx, refNoisyPat] = InitialBlockMatching(medfilt2(X), Y, Par);

for i = 1:maxiter
    if i ~= 1
        X = XSubproblem( Y+W./beta-E, Par, NL_mat );
    end
    energy(i) = ComputeEnergy( refIdx, refNoisyPat, X, Par, lambda, gamma );
    if i ~= 1
        De = abs(energy(i-1)-energy(i))/energy(i);
        p2 = psnr(uint8(X), uint8(img)); s2 = ssim(uint8(X), uint8(img));
        fprintf("iteration : %d, PSNR : %.2f, SSIM : %.4f\n", i, p2, s2);
        if De < Par.threshold
            p1 = psnr(uint8(Y), uint8(img)); s1 = ssim(uint8(Y), uint8(img));
            figure;
            subplot(1,2,1); imshow(uint8(Y)); title({'Noisey image';['PSNR = ',num2str(p1)];['SSIM = ', num2str(s1)]});
            subplot(1,2,2); imshow(uint8(X)); title({'Recovery image';['PSNR = ',num2str(p2)];['SSIM = ', num2str(s2)]});
            break
        end
    end
    for j = 1:5
        numer = lambda*E./(gamma^2+E.^2) + beta*(E-Y-W./beta+X);
        denomi = lambda*(gamma^2 - E.^2)./(gamma^2+E.^2).^2 + beta;
        E = E - numer./(denomi+eps);
    end
    W = W + beta*(Y-X-E);
end
