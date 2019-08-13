% choose image by number
% All images are in the folder "test_images"
for globalit = 10

% solve WNNM model of the Cauchy noise denoising by ADMM algorithm.

if globalit == 1; clear; imname = "barbara"; img = double(imread('barbara.png')); img = imresize(img, 0.5);
elseif globalit == 2; clear; imname = "boat"; img = double(imread('boat.png')); img = imresize(img, 0.5);
elseif globalit == 3; clear; imname = "cameraman"; img = double(imread('cameraman.tif'));
elseif globalit == 4; clear; imname = "house"; img = double(rgb2gray(imread('house.tiff'))); img = imresize(img, 0.5);
elseif globalit == 5; clear; imname = "man"; img = double(imread('man.tiff')); img = imresize(img, 0.25);
elseif globalit == 6; clear; imname = "mandrill"; img = double(imread('mandrill.png')); img = imresize(img, 0.5);
elseif globalit == 7; clear; imname = "peppers"; img = double(imread('peppers256.png'));
elseif globalit == 8; clear; imname = "plane"; img = double(imread('jetplane.png')); img = imresize(img, 0.5);
elseif globalit == 9; clear; imname = "synthetic"; img = double(rgb2gray(imread('synthetic1.bmp'))); img = imresize(img, 1.01);
elseif globalit == 10; clear; imname = "walkbridge"; img = double(imread('walkbridge.png')); img = imresize(img, 0.5);
end

rng('default');
eta1 = randn(size(img));
eta2 = randn(size(img));
gamma = 10;
noise = gamma.*eta1./eta2;
Y = img + noise;

fprintf('noisy image PSNR, SSIM : %.2f, %.4f\n', psnr(uint8(Y), uint8(img)), ssim(uint8(Y), uint8(img)));

maxiter = 100;
lambda = 1;
beta = 2*lambda/gamma^2 + eps;
W = zeros(size(img));
medY = medfilt2(Y);
X = max(min(medY(:)),min(Y,max(medY(:))));

sigma = 4;
Par = ParSet(sigma);

[NL_mat, refIdx, refNoisyPat] = InitialBlockMatching(medfilt2(X), Y, Par);

energy = zeros(1,maxiter+1);

E = zeros(size(img));
for i = 1:maxiter
    if i ~= 1
        X = XSubproblem( Y+W./beta-E, Par, NL_mat );
    end
    energy(i) = ComputeEnergy( refIdx, refNoisyPat, X, Par, lambda, gamma );
    if i ~= 1
        De = abs(energy(i-1)-energy(i))/energy(i);
        p2 = psnr(uint8(X), uint8(img)); s2 = ssim(uint8(X), uint8(img));
        fprintf("i, current PSNR, SSIM : %d, %.2f, %.4f\n", i, p2, s2);
        if De < 5e-4
            p1 = psnr(uint8(Y), uint8(img)); s1 = ssim(uint8(Y), uint8(img));
            fprintf("i, PSNR of recovered image, SSIM : %d, %.2f, %.4f\n", i, psnr(uint8(X), uint8(img)), ssim(uint8(X), uint8(img)));
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
end