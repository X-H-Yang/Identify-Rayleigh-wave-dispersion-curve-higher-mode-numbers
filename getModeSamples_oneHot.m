function [samples_x,samples_y] = getModeSamples_oneHot(samples_N,higherModeNum,...
    unique_2_freq,f_00,f_0i,Vp_true,den_true,Vs_profile_lower,Vs_profile_upper)

%% generate samples
% lower and upper boundary of samples
x_low = Vs_profile_lower(:)';
x_up = Vs_profile_upper(:)';

X_raw = zeros(higherModeNum,length(f_00)+length(f_0i));
Y_raw = zeros(higherModeNum,higherModeNum); % one-hot encoding

[~, match_indices_fun] = ismember(f_00, unique_2_freq);
[~, match_indices_high] = ismember(f_0i, unique_2_freq);

% targetHigerModes = 1:1:higherModeNum;

num_unknown = length(Vs_profile_lower);
N_true = (num_unknown - 1)/2;

i = 1;
while i <= samples_N
    for j = 1:num_unknown
        Vs_profile_raw(i,j) = x_low(j) + (x_up(j)-x_low(j))*rand(1,1);
    end

    h = [Vs_profile_raw(i,1:N_true) 0];
    Vp = Vp_true;
    Vs = Vs_profile_raw(i,N_true+1:2*N_true+1);
    den = den_true;
    
    model_dispersion_R = [];
    
    try
        out = gpdc(h,Vp,Vs,den,'fV',unique_2_freq);
        out2 = rdivide(1, out(:, 2:end));
        
        for kk = 1:1:higherModeNum
            temp = [out2(match_indices_fun,1); out2(match_indices_high,kk+1)];
            temp(isnan(temp)) = 0;
            X_raw(kk,:) = temp';
            Y_raw(kk,kk) = 1; % one-hot encoding
        end
        
        if i == 1
            samples_x = X_raw;
            samples_y = Y_raw;
        else
            samples_x = [samples_x; X_raw];
            samples_y = [samples_y; Y_raw];
        end
        
        i = i+1;

    catch
            %
            %
    end
    
end


end

