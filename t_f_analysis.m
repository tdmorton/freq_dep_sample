%% Wavelet Plotter
clearvars
close all

[a, Fs] = audioread('speaking.wav');
x = a;
rng default

SNR = 40;
y = randn(size(x))*std(x)/db2mag(SNR);

s = x + y;
% a = a(1:Fs);

% a = s;

No = 6;
Nv = 32;

figure()
[wt, f] = cwt(a,Fs);%,'NumOctaves', No, 'VoicesPerOctave',Nv);
cwt(a,Fs);%,'NumOctaves', No, 'VoicesPerOctave',Nv)
title('CWT Spectrogram')

%% Simple Thresholding Technique
cutoff_num = 20;
cutoff_Freq = f(cutoff_num);
d_samp_amount = 4;


thresh = 0.1; % less than 10 percent of the signal's power is above the cutoff_Freq
tot_samples = size(a,1);
flag = zeros([tot_samples,1]);

for i = 1:tot_samples
    if sum(abs(wt(1:cutoff_num,i))) < thresh*sum(abs(wt(:,i)))
        flag(i) = 1;
    end
end

one_new_signal = downsample_flag(a,flag);
% 
% figure()
% plot(new_signal)

one_new_new_signal = upsample_flag(one_new_signal,flag);

% figure()
% plot(new_new_signal)

        

% 
% b = linspace(0, length(a)/Fs, length(a));
% figure()
% plot(b,a)
b = linspace(0, length(a)/Fs, length(a));
figure()
plot(a)

%% Candy Wine test

% [k, F] = audioread('c_wine.mp3');
% figure()
% k = k(1:Fs*10,1);
% % sound(k,F)
% cwt(k,F)


% % [sin, Fsin] = audioread('test_sin.wav');
% Fsin = 44100;
% ff = 1.2E3;
% T = 1;
% t = 0:1/Fsin:T;
% 
% sin_wave = sin(2.*pi.*ff.*t);
% 
% figure()
% plot(t,sin_wave)
% 
% figure()
% cwt(sin_wave,Fsin,'NumOctaves', No, 'VoicesPerOctave',Nv);
% [wtsin, fsin] = cwt(sin_wave,Fsin,'NumOctaves', No, 'VoicesPerOctave',Nv);
% 

%% Multi-Level Thresholding
factors = factor(Fs);
threshold_levels = 7;
% 
for c = 2:threshold_levels+1
    levels(c-1) = Fs/lin_multiply(factors(1:c));
end

V = zeros([length(f), 1]);
V(1:length(levels)) = levels.';
N = f;

A = repmat(N,[1 length(V)]);
[minValue,closestIndex] = min(abs(A-V'));
closestValue = N(closestIndex) ;

cutoff_freqs = closestValue(1:length(levels));
cutoff_nums = closestIndex(1:length(levels));

d_samp_amounts = Fs./levels/4;

% flag = zeros([length(a),length(levels)]);

thresh = 0.02; % less than 10 percent of the signal's power is above the cutoff_Freq
tot_samples = size(a,1);
flag = zeros([tot_samples,1]);
mag_thresh = 0.1;

for qq = 1:length(levels)
for i = 1:tot_samples
    if  sum(abs(wt(:,i))) < mag_thresh || sum(abs(wt(1:cutoff_nums(qq),i))) < thresh*sum(abs(wt(:,i)))
        flag(i) = qq;
    end
end
end

figure()
plot(flag)
title('Flag')

new_signal = downsample_flag(a,flag);

figure()
plot(new_signal)
title('Compressed Signal')

new_new_signal = upsample_flag(new_signal,flag);

figure()
plot(new_new_signal)
title('Upsampled Compressed Signal')


%% Noise Evaluation
a = a.';
% new_new_signal = a;

Haar_99 = load('haar_5_95.76.mat'); % 5th order
db_99 = load('db_1_5_95.76.mat'); % 5th order , 1
sym = load('sym_2_5_96.72.mat'); % 



tests = [one_new_new_signal; new_new_signal; Haar_99.haar_5_95; db_99.db_1_5_95; sym.sym_2_5_96 ];


SNR = zeros([size(tests,1),1]);
PSNR = zeros([size(tests,1),1]);
NRMSE = zeros([size(tests,1),1]);
CR = zeros([size(tests,1),1]);

for i = 1:size(tests,1)

squared_error = abs(a-tests(i,:)).^2;
MSE(i) = sum(squared_error(:))./length(a);

%SNR 10log10(meanorig^2/meanorig-recon^2)
SNR(i) = 10.*log10(rms(a).^2./rms(tests(i,:)).^2);

%PSNR 10log10(NX^2/|x-r|^2)
PSNR(i) = 10*log10(length(a)*max(a)^2/MSE(i));

%NRMSE
NRMSE(i) = sqrt((a-tests(i,:)).^2/(a-rms(a)).^2);

%CR
CR(i) = length(a)/length(new_signal);

end

Compression_Method = {'1-Bit Freq-Dep Downsample'; '3-Bit Freq-Dep Downsample'; 'Haar Wavelet'; 'DB Wavelet'; 'Sym Wavelet'};

T = table(Compression_Method, SNR, PSNR, NRMSE);
T

