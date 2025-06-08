% Code package purpose: Accoplish the proposed method in the paper 
%                       entitled "Identification of Higher-mode Numbers 
%                         in Dispersion Curves for Rayleigh Wave Inversion".
% About inversion codes: The inversion program is hosted on the author's
%                          GitHub repository (https://github.com/X-H-Yang).
%
% paper status: major revision - submitted 
%               (journal: IEEE Transactions on Geoscience and Remote Sensing)
% Xiao-Hui Yang, Peng Han, Jiancang Zhuang, Yuanyuan Zhou, Gexue Bai, 
%                              Ruidong Li, Wuhu Zhang, Xiaofei Chen (2025). 
% Identification of Higher-mode Numbers in Dispersion Curves 
%                                              for Rayleigh Wave Inversion. 
% IEEE Transactions on Geoscience and Remote Sensing, major revision.
%
% software version: MATLAB R2017a
%
% Acknowledgement: The forward modeling program used to generate 
%                  theoretical Rayleigh wave dispersion curves in this 
%                  code package was obtained from the  website 
%                  (https://github.com/eespr/MuLTI) provided by 
%                  Killingbeck et al. (2018); the broad learning network 
%                  codes available on the website 
%                  (https://broadlearning.ai/, Chen and Liu 2017 and Chen 
%                  et al. 2018) were also applied for the accomplishment of
%                  this code package.
%
% Killingbeck et al. (2018): Killingbeck, S. F., Livermore, P. W., 
%                            Booth, A. D., & West, L. J. (2018). Multimodal 
%                            layered transdimensional inversion of seismic 
%                            dispersion curves with depth constraints. 
%                            Geochemistry, Geophysics, Geosystems, 19(12), 
%                            4957-4971.
%
% Chen and Liu 2017: Chen, C. P., & Liu, Z. (2017). Broad learning system: 
%                    An effective and efficient incremental learning system
%                    without the need for deep architecture. IEEE 
%                    transactions on neural networks and learning systems, 
%                    29(1), 10-24.
%
% Chen et al. 2018: Chen, C. P., Liu, Z., & Feng, S. (2018). Universal 
%                   approximation capability of broad learning system and 
%                   its structural variations. IEEE transactions on neural 
%                   networks and learning systems, 30(4), 1191-1204.
%
% Input Data:
%   Known fundamental mode data:
%    - File: 'numerical_fun.xls'
%    - Contains: Frequency (Hz) and Phase Velocity (m/s) columns
%    - Assumption: This data represents the fundamental mode dispersion curve
%
%  Higher-order mode data to be identified:
%    - Case 1: 'numerical_A.xls' (known to be 1st higher mode in validation)
%    - Case 2: 'numerical_B.xls' (known to be 3rd higher mode in validation)
%
% Actual earth model parameters of this numerical example:
%     h_true = [3 2 3 4]; % layer thickness (all layers)
%     Vs_true = [200 130 220 300 450]; % shear wave velocity (all layers)
%     Vp_true = [400 260 440 600 900]; % primary wave velocity (all layers)
%     den_true = [1.75 1.70 1.75 1.80 1.90]; % density (all layers)
%
% Variables about training samples:
%       samples_x: input samples containing training and validation samples
%       samples_y: output samples containing training and validation samples
%       train_x_norm: normalized input training dataset
%       train_y_norm: normalized output training dataset
%       validation_x_norm: normalized input validation dataset
%       validation_y_norm: normalized output validation dataset
%       test_x: test dataset (only one sample, combining fundamental mode 
%                                      and a higher mode  to be identified)
%       test_x_norm: normalized input test dataset
%
% Date: 2025/06/08
%
% Developed by: Xiao-Hui Yang, Currently working at 
%               Chengdu University of Information Technology
%
% Email: yangxh@cuit.edu.cn / xiao-hui.yang@hotmail.com
%
% Note: the specific descriptions of the proposed mode identification method
%    for Rayleigh wave inversion can refer to the paper entitiled
%   "Identification of Higher-mode Numbers in Dispersion Curves for Rayleigh
%    Wave Inversion."; the users can cite this paper for scientific research.
% 
% Xiao-Hui Yang, Peng Han, Jiancang Zhuang, Yuanyuan Zhou, Gexue Bai, 
%                              Ruidong Li, Wuhu Zhang, Xiaofei Chen (2025). 
% Identification of Higher-mode Numbers in Dispersion Curves 
%                                              for Rayleigh Wave Inversion. 
% IEEE Transactions on Geoscience and Remote Sensing, major revision.
%

clear;
clc;
close all;

% Seeds the random number generators using system time
rand('seed',sum(100*clock))
randn('seed',sum(100*clock))

%% Raw data
curve_00 = xlsread('numerical_fun.xls'); % read fundamental mode data 
curve_A = xlsread('numerical_A.xls'); % Case 1: a higher mode to be identified
curve_B = xlsread('numerical_B.xls'); % Case 2:  a higher mode to be identified

% fundamental mode
f_00_original = curve_00(:,1)'; f_00_original = f_00_original(:)';
dispersion_00_original = curve_00(:,2); dispersion_00_original = dispersion_00_original(:)';
% Case 1: a higher mode to be identified
f_A_original = curve_A(:,1)'; f_A_original = f_A_original(:)';
dispersion_A_original = curve_A(:,2); dispersion_A_original = dispersion_A_original(:)';
% Case 2: a higher mode to be identified
f_B_original = curve_B(:,1); f_B_original = f_B_original(:)';
dispersion_B_original = curve_B(:,2); dispersion_B_original = dispersion_B_original(:)';

%% Interpolation process - frequency and dispersion values
% the following frequency range should be modified for new data
f_00_min = 5; f_00_max = 50.5; % fundamental mode (frequency range)
f_A_min = 46.5; f_A_max = 85; % case 1: to be identified (frequency range)
f_B_min = 54; f_B_max = 100; % case 2: to be identified (frequency range)

df = 0.5; % frequency interval

f_00 = f_00_min:df:f_00_max; % fundamental mode (frequency points)
f_A = f_A_min:df:f_A_max; % case 1: to be identified (frequency points)
f_B = f_B_min:df:f_B_max; % case 2: to be identified (frequency points)
% fundamental mode (dispersion values)
dispersion_00 = interp1(f_00_original,dispersion_00_original,f_00);
% case 1: to be identified (dispersion values)
dispersion_A = interp1(f_A_original,dispersion_A_original,f_A);
% case 2: to be identified (dispersion values)
dispersion_B = interp1(f_B_original,dispersion_B_original,f_B);

%% Define a target higher mode to be identified

%%%%%%%%%%%%%%%%%%% a higher mode to be identified  %%%%%%%%%%%%%%%%%%%%%%%
f_modeUnknown = f_A; % a higher mode to be determined
test_x = [dispersion_00 dispersion_A]; % same mode with f_modeUnknown
% comment this section when you need to identify another higher mode
% uncomment this section when you need to identify this higher mode
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%% a higher mode to be identified  %%%%%%%%%%%%%%%%%%%%%%%
% f_modeUnknown = f_B; % a higher mode to be determined
% test_x = [dispersion_00 dispersion_B]; % same mode with f_modeUnknown
% % comment this section when you need to identify another higher mode
% % uncomment this section when you need to identify this higher mode
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Input parameters required before network training
higherModeNum = 3; % the maximum possible mode number of higher modes
                   % (e.g., only 1st, 2nd, and 3rd higher modes considered)
                   % note: current foward modeling package allows for the
                   %       maximum possible mode number being 4

% combine frequencies: fundamental mode and a higher mode to be identified
combined_2_freq = [f_00 f_modeUnknown]; % combine frequencies
unique_2_freq = unique(combined_2_freq); % remove the repeated frequencies
f_all = {unique_2_freq; f_00; f_modeUnknown}; % use for easy sampling

Vp_true = [400 260 440 600 900]; % primary wave velocities, assumed known
den_true = [1.75 1.70 1.75 1.80 1.90]; % densities, assumed known
% search space for random sample generation
Vs_profile_lower = [1 1 1 1 100 100 100 200 300];
Vs_profile_upper = [4 4 5 5 300 300 300 400 550];

% hyperparameter setting for network complexity selection (BL network)
Fea_vec = 4:2:10; % hyperparameter 1
Win_vec = 4:2:20; % hyperparameter 2
Enhan_vec = 4:2:30; % hyperparameter 3

% define the number of samples in training and validation datasets
train_samples_N = 500; % number of training samples
validation_samples_N = train_samples_N*0.2; % number of validation samples
samples_N = train_samples_N + validation_samples_N; % number of total samples

%% Sample generation and normalization
tic % Start timing execution

% generate random samples
disp('-------------------------------------------------------------------')
disp('Waiting about one minute for generating one-hot encoding samples ...')
[samples_x,samples_y] = getModeSamples_oneHot(samples_N,higherModeNum,...
    f_all{1},f_all{2},f_all{3},Vp_true,den_true,Vs_profile_lower,...
    Vs_profile_upper); % sample generation: one-hot encoding
disp('training and validations samples were obtained!')

% Normalize data (mean=0, std=1)
[mean_x,std_x,samples_x_norm] = normalized_fun(samples_x);
samples_y_norm = samples_y; % one-hot encoding, omitting normalization

% training and validation datasets
train_x_norm = samples_x_norm(1:train_samples_N*higherModeNum,:);
validation_x_norm = samples_x_norm(train_samples_N*higherModeNum+1:end,:);
train_y_norm = samples_y_norm(1:train_samples_N*higherModeNum,:);
validation_y_norm = samples_y_norm(train_samples_N*higherModeNum+1:end,:);

% test dataset - observed higher mode to be identified
test_x_norm = zeros(size(test_x,1),size(test_x,2));
for i = 1:1:size(test_x,2)
    test_x_norm(:,i) = (test_x(:,i)-mean_x(i))/std_x(i);
end

%% Train BLS networks and perform complexity selection for predicion
disp('-------------------------------------------------------------------')
disp('Performing network training for mode number identification ...')
disp('Network complexity selection space:')
disp('             -- hyperparameter 1 -> [4:2:10]');
disp('             -- hyperparameter 2 -> [4:2:20]');
disp('             -- hyperparameter 3 -> [4:2:30]');
[Y_hat_norm,result,all_index,NumFea_hat,NumWin_hat,NumEnhan_hat] = ...
    bls_classification_Y_sub_withTest(train_x_norm,train_y_norm,...
    validation_x_norm,validation_y_norm,Fea_vec,Win_vec,Enhan_vec,test_x_norm);
disp('Network modeling finished!')
disp('-------------------------------------------------------------------')

Y_hat = Y_hat_norm; % no re-scale because of omitting normalization

% softmax output
exp_values = exp(Y_hat - max(Y_hat));
Y_hat_softmax = exp_values / sum(exp_values);

toc % terminate timing measurement

%% Output: identified results
% display the raw output of the network
disp('Network output:');
disp(Y_hat);
% display the softmax-processed output, representing class probabilities
disp('Softmax-processed output (class probabilities):');
disp(Y_hat_softmax);
% output identified mode number 
[maxVal,maxIdx] = max(Y_hat_softmax);
if maxIdx == 1
    disp(['The identified result of the input data: ',...
        '1st higher mode']);
end
if maxIdx == 2
    disp(['The identified result of the input data: ',...
        '2nd higher mode']);
end
if maxIdx == 3
    disp(['The identified result of the input data: ',...
        '3rd higher mode']);
end
disp('-------------------------------------------------------------------')

% plot loss function curve for BL Network complexity selection
myFontSize = 20;
figure()
plot(all_index(:,1),all_index(:,2),'g','Linewidth',1.5)
xlabel('Time [s]','FontSize',myFontSize);
ylabel('Cross entropy','FontSize',myFontSize);
set(gca,'FontName','Times New Roman','FontSize',myFontSize);
xlim_min = 0;
xlim_max = ceil(max(all_index(:,1))/10)*10 + 10;
ylim_min = 0.8;
ylim_max = 1;
axis([xlim_min xlim_max ylim_min ylim_max]);
disp('Work done!');