clc; clear all; close all;

addpath(genpath('..\shared'))

global  ratio file_path_rgb256;
global count time curr_d gt256;
global sensor sate num mu tag ;
global L Qblocks_size flag_cut_bounds dim_cut thvalues;

curr_d = pwd;
sate = 'ik';  % geo,ik,pl, qb, wv2£¬wv3  pl-OK

method = 'PCRF';
%% start initialization
initialize_sate_RS_MTF();

%% start process
% for  num= 1 : length(file_path_rgb256)
for  num= 15:15
    count = 1 + count
    %% read images and preprocess
    [rgb256, gt256, mul64, pan256 ] = read_image_sate_RS(num);
    
    %% tic start
    tic;
    ORM = gt256;
    gt_dou = im2double(gt256);
    img_mul = im2double(mul64);
    P = im2double(pan256);
    %% get I component
    [m,n] = size(P);
    M =imresize(img_mul,size(P)); 
    
    I = get_I(sate,M, P);
    P=(P-mean(P(:)))*std2(I)/std(P(:)) + mean2(I);   % histogram matching

    %% culculate frequency doman laplacian
    PQ = size(P); % size(f)=[256 256]
    switch sate
        case 'ik'

          u_crf = 3;  
            u_flt = 0.1; 
            gama_lagr = 0.00005;           
            alpha = 0.95;
            sigma = 15; 
            H = lpfilter('gaussian',PQ(1),PQ(2),sigma);  

        case 'pl'
            u_crf = 5;
            u_flt = 0.01;
            gama_lagr = 0.01;
            alpha = 0.9;  
        case 'wv2'
            %            u = get_u(P,H,I);
            u_crf = 6;
            u_flt = 0.1;
            gama_lagr = 0.003;
            alpha = 1.4;      
            sigma = 20;  % wv2  26, ik 17
            H = lpfilter('gaussian',PQ(1),PQ(2),sigma);      % in WV is OK,not in iK
    % *************original      
%             u_crf = 5;
%             u_flt = 0.1;
%             gama_lagr = 0.01;
%             alpha = 1.5;      
        case 'wv3'
            u_crf = 6;
            u_flt = 0.1;
            gama_lagr = 0.003;
            alpha = 1.4;    % original 1.3
            sigma = 20;  % wv2  26, ik 17
            H = lpfilter('gaussian',PQ(1),PQ(2),sigma);      % in WV is OK,not in iK
            % *******original
%              u_crf = 5;
%             u_flt = 0.1;
%             gama_lagr = 0.001;
%             alpha = 1.3;
        case 'qb'
            u_crf = 3;
            u_flt = 0.01;
            gama_lagr = 0.01;
            alpha = 1;
        case 'geo'
            alpha = 1.2;
            u_crf = 50;
            u_flt = 0.1;
            gama_lagr = 0.05;
    end
    
    
    %% start
    lap = fspecial('laplacian',0);
%     lap = fspecial('prewitt');
    lap_f = freqz2(lap,size(P)); % laplase matrix
    lap_f = fftshift(lap_f);

    %% initial
    u2_flt = 0;         
    sig = 0;        
%     not in WV2 WV3
    %% CRFs model
    eps = 1*10^-3;
    iter = 0;
    lagr_m = ones(m);                % lagrange factor
    u_aug = 1;             % augment lagrange parameter   
    rho = 1.01;
    X = 0;
    IH3_old = I;
    IH3_new = P;
    N_I = norm(IH3_new-IH3_old)/norm(IH3_old);
    while (N_I>eps)
        IH3_old = IH3_new;
        iter = iter + 1
        A = H.*fft2(I) + u_crf*lap_f.*lap_f.* fft2(P)+lap_f.*fft2(lagr_m)+u_aug*lap_f.*fft2(X);
        B = H.*H + u_crf*lap_f.*lap_f + u_aug*lap_f.*lap_f;        
        C = A./B;
        IH3_new = real(ifft2(C));
%         H = filter_est(IH3_new,I,u2_flt,u_flt);
        X = shrink(real(ifft2(lap_f.*C))-1/u_aug*lagr_m, gama_lagr/u_aug);
        lagr_m = lagr_m +u_aug*(X - real(ifft2(lap_f.*C)));
        u_aug = rho*u_aug;
        N_I = norm(IH3_new-IH3_old)/norm(IH3_old);
        N_I_all(iter) = N_I;
    end
    De = IH3_new - I;
    
    M_sum = M(:,:,1)+M(:,:,2)+M(:,:,3)+M(:,:,4);
    M_sum(find(M_sum==0))=eps;
    M_rate = zeros(size(M));
    %     F1 = zeros(size(M));
    for i=1:4
        M_rate(:,:,i) = 4*M(:,:,i)./M_sum;
    end       

    
    %% final image
    F_final = zeros(size(M));
    for i = 1:4

        F_final(:,:,i) = M(:,:,i)+ alpha*M_rate(:,:,i).*De;
    end
    
    toc
    time(num) = toc;    

    F_rgb = zeros(PQ(1),PQ(2),3);
    for i = 1:3
        F_rgb(:,:,i) =  F_final(:,:,i);
    end

    
    %% save images
%     save_sate_RS(F_final, M(:,:,1:3), F_rgb, num);
end
%% show image
Eval = Evaluation4(ORM,pan256,uint8(F_final*255))
%     saveData_RS(sig_all,co_PLI, co_PI);
[Q_avg_Segm, SAM_Segm, ERGAS_Segm, SCC_GT_Segm, Q_Segm] = indexes_evaluation(...
F_final,im2double(ORM),ratio,L,Qblocks_size,flag_cut_bounds,dim_cut,thvalues);
Eval2 = [Q_avg_Segm, SAM_Segm, ERGAS_Segm, SCC_GT_Segm, Q_Segm]

figure, imshow(F_rgb)
figure, imshow(rgb256)
% cor
% figure, imshow(I)
% P_L = ifft2(fft2(P).*H);
% figure, imshow(P_L)
% corr2(P,I)
% corr2(P_L,I)
% figure;

%% write result
T = quality_eval_sate_RS();
cd(curr_d);
writetable(T,strcat('final_result/quality_',method,'_4C_',sate,'_RS_sig.csv'),'WriteRowNames',true);

