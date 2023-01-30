clc; clear all; close all;
%% only for plaidas
addpath(genpath('E:\remote sense image fusion\shared'))

global   thvalues  ratio L;
global   im_tag sensor;
global   file_path_rgb_noR;
global  count time num ;
global  sate curr_d alpha;

curr_d = pwd;
sate = 'wv2';   % geo,ik,pl, qb, wv2，wv3

method = 'PCRF';
S_block = 128;
%% start initialization
initialize_sate_FS();

%% start process
% for  num= 1 : length(file_path_rgb_noR)t
for  num= 5:5
    count = 1 + count
    %% read images and preprocess
    [mul_noR, pan_noR ] = read_image_sate_FS(num);
    
    %% tic start
    tic()
    %     ORM = gt256;
    img_mul = im2double(mul_noR);
    P = im2double(pan_noR);
    %% 加权平均求I分量
    [m,n] = size(P);
    M =imresize(img_mul,size(P),'bicubic');   % 双三次插值算法
%     alpha =impGradDes(M,P);   % gradient desend  method to calculate coefficient alpha with min||P-sigmoid alpha*M ||
        alpha(1:4) = 1/4;
    I=alpha(1)*M(:,:,1)+alpha(2)*M(:,:,2)+alpha(3)*M(:,:,3)+alpha(4)*M(:,:,4);
    %% histogram matching
    P=(P-mean(P(:)))*std2(I)/std(P(:)) + mean2(I);   % histogram matching
    %     I_mul = 1/4*(img_mul(:,:,1)+img_mul(:,:,2)+img_mul(:,:,3)+img_mul(:,:,4));
    
    M_sum = M(:,:,1)+M(:,:,2)+M(:,:,3)+M(:,:,4);
    M_sum(find(M_sum==0))=eps;
    M_rate = zeros(size(M));
    %     F1 = zeros(size(M));
    for i=1:4
        M_rate(:,:,i) = 4*M(:,:,i)./M_sum;
    end
    %     g = max(0.86, corr2(P,I));
    
    %% block processing
    in3 = cat(3,P,I);
    fun_crf = @(bs)D_CRFs(bs.data);
    De = blockproc(in3,[S_block, S_block],fun_crf,'BorderSize',[2,2]);
    %%  Computing mean of std. devs.
    %%% Histogram matching of each MS band to Pan
    %     msexp_hm = zeros(N,M,L);
    %     for k=1:4
    %         b = M(:,:,k);
    %         b = (b - mean2(b) + mean2(P)/std2(P)*std2(b)) * std2(P)/std2(b);
    %         b(b<0) = 0;
    % %         msexp_hm(:,:,k) = b;
    %     end
    %     PL = ifft2(fft2(P).*H);
    %     HL = ifft2(fft2(IH3_new).*H);
    %     aux3 = zeros(1,4);
    %     for k=1:4
    %         aux3(k) = std2(M(:,:,k));
    %     end
    %     aux3 = mean(aux3);
    %
    %     beta = 1; % for 11-bit data
    %     % beta = 1.95; % for 8-bit data
    %     %%% Computing weights
    %     w = zeros(1,4);
    %     for k=1:4
    %         aux1 = PL;
    %         b = M(:,:,k);
    %         w(k) = beta.* corr2(aux1(:),b(:))*std(b(:))/aux3;%std(aux2(:));
    %     end
    %% get edge weight
%     lamda = 10^-9; eps = 10^-10; bata = 0.25;
%     Ep = expEdge(P,lamda,eps);
%     for i=1:4
%         Em(:,:,i) = expEdge(M(:,:,i),lamda,eps);
%     end
%     
%     for i=1:4
% %         W(:,:,i) = 3*M(:,:,i).*(bata*Ep+(1-bata)*Em(:,:,i))./(M(:,:,1)+M(:,:,2)+M(:,:,3));
%         W(:,:,i) = 4*M(:,:,i).*(bata*Ep+(1-bata)*Em(:,:,i))./(M(:,:,1)+M(:,:,2)+M(:,:,3)+M(:,:,4));
%     end
%% get final image
    F_final = zeros(size(M));
    for i = 1:4
        %     F_final(:,:,i) = mul_re(:,:,i)+M_rate(:,:,i).*De;
        %                 F_final(:,:,i) = M_rate(:,:,i).*IH3;
        F_final(:,:,i) = M(:,:,i)+M_rate(:,:,i).*De;

%         F_final(:,:,i) = M(:,:,i)+W(:,:,i).*De;
        %          F_final(:,:,i) = M(:,:,i)+w(i)*M_rate(:,:,i).*De;
        %          F_final(:,:,i) = mul_re(:,:,i)+alpha(i)/mean(alpha)*M_rate(:,:,i).*I_H;
    end
    
    toc
    time(num) = toc;
    %     IH2_L = ifft2(fft2(IH2).*H);
    F_rgb = zeros(m,n,3);
    for i = 1:3
        F_rgb(:,:,i) =  F_final(:,:,i);
    end
    
    %     sig_all(num) = sig;
    %     co_PLI(num) = cor;
    %     co_PI(num) = corr2(P,I);
    %      HL = ifft2(fft2(IH3_new).*H);
    %       PL = ifft2(fft2(P).*H);
    %       cor = corr2(PL,I);
    %       cor2 = corr2(HL,I);
    
    save_sate_FS(F_final, M,F_rgb, num);
end
%% show image
figure, imshow(F_rgb);
[Dl,Ds,QNR_index,SAM_index,sCC] = indexes_evaluation_FS(F_final,img_mul,P,...
    L,thvalues,M,sensor,im_tag,ratio);
% [D_lambda,D_S,QNR_index,SAM_index,sCC] = indexes_evaluation_FS(Fused,I_MS_LR,I_PAN);
Eval = [Dl,Ds,QNR_index,SAM_index,sCC]
% corr2(P,I)

% saveData_FS(sig_all,co_PLI, co_PI);
% figure, imshow(I)
% P_L = ifft2(fft2(P).*H);
% % corr2(I,P_L)
% figure, imshow(P_L)

%% write result
T = quality_eval_sate_FS();
cd(curr_d);
writetable(T,strcat('final_result/quality_',method,'_4C_',sate,'_FS_block_5.csv'),'WriteRowNames',true);

%%  function of getting details
function De = D_CRFs(in3)

global sate alpha;
%% culculate frequency doman laplacian
P = in3(:,:,1);
I = in3(:,:,2);
PQ = size(P); % size(f)=[256 256]
%  [H, sig, cor] = get_gau_H_RS(P,I);
P=(P-mean(P(:)))*std2(I)/std(P(:)) + mean2(I);   % histogram matching
%  u = (1-cor)*1000;
[m,n] = size(P);
switch sate
    case 'geo'
        %         u_crf = 5000;
        %         u_flt = 1;
        %         gama_lagr =0.8;
%                 alpha = 1;
        
        u_crf = 5;
        u_flt = 0.1;
        gama_lagr = 0.1;
        alpha = 1;  
    case 'ik'
%         u_crf = 10000;  %% 越小越模糊
%         u_flt = 1;  %% 越小越模糊
%         gama_lagr = 1;
%         alpha = 0.9;
        
        u_crf = 5;
        u_flt = 0.1;
        gama_lagr = 0.08;
        alpha = 1;  
    case 'pl'
        u_crf = 5;
        u_flt = 0.1;
        gama_lagr = 0.05;
        alpha = 1;
    case 'qb'
        u_crf = 5;
        u_flt = 0.1;
        gama_lagr = 0.08;
        alpha = 1;
    case 'wv2'
%         u_crf = 5000;
%         u_flt = 1;
%         gama_lagr = 1;
%         alpha = 1;
        
        u_crf = 6;
        u_flt = 0.1;
        gama_lagr = 0.07;
        alpha = 1.2;
    case 'wv3'
%         u_crf = 5000;
%         u_flt = 1;
%         gama_lagr = 0.1;
%         alpha = 1;
        
        u_crf = 5;
        u_flt = 0.1;
        gama_lagr = 0.1;
        alpha = 1.2;
        
end

%% start
lap = fspecial('laplacian',0);
lap_f = freqz2(lap,size(P)); % laplase matrix
lap_f = ifftshift(lap_f);

%% initial
u2_flt = 0;
sigma =  80;   % changed
H = lpfilter('gaussian',PQ(1),PQ(2),sigma);

%% CRFs model
eps = 2*10^-3;
iter = 0;
lagr_m = ones(m);                % lagrange factor
u_aug = 1;             % augment lagrange parameter
rho = 1.01;
X = 0;
IH3_old = I;
IH3_new = P;
while (norm(IH3_new-IH3_old)/norm(IH3_old)>eps)
    IH3_old = IH3_new;
    iter = iter + 1
    A = H.*fft2(I) + u_crf*lap_f.*lap_f.* fft2(P)+lap_f.*fft2(lagr_m)+u_aug*lap_f.*fft2(X);
    B = H.*H + u_crf*lap_f.*lap_f + u_aug*lap_f.*lap_f;
    C = A./B;
    IH3_new = real(ifft2(C));
    H = filter_est(IH3_new,I,u2_flt,u_flt);
    X = shrink(real(ifft2(lap_f.*C))-1/u_aug*lagr_m, gama_lagr/u_aug);
    lagr_m = lagr_m +u_aug*(X - real(ifft2(lap_f.*C)));
    u_aug = rho*u_aug;
end
De = alpha*(IH3_new - I);
end
