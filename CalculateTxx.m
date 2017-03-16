function Txx = CalculateTxx(img_rs,BlockSize)
% Wrapper for the almost_cyclo_test function which is provided by 
% D. Vazquez-Padin
% 

N=BlockSize;
P=21;
K=9;
Pf=1e-5;

B=CentralCrop(double(img_rs),N+sqrt(K)) ;
% Apply the Laplacian operator
X=del2(B);
% X=B;

% Run the statistical test
[Txx,Tgamma]=almost_cyclo_test(X,N,P,K,Pf);
Txx = abs(Txx);
end

function [Txx,Tgamma]=almost_cyclo_test(X,N,P,K,Pf)

%--------------------------------------------------------------------------
% Core function of the method proposed in:
%
% D. Vazquez-Padin, C. Mosquera and F. Perez-Gonzalez.
% "Two-dimensional statistical test for the presence of almost
% cyclostationarity on images",
% 2010 17th IEEE International Conference on Image Processing (ICIP),
% pp.1745-1748, 26-29 Sept. 2010
%
%--------------------------------------------------------------------------
%
% Inputs:
% -------
% X:  image block
% N:  specifies the N-by-N size of the 2D-FFT grid
% P:  specifies the P-by-P size of the spectral window (P must be odd)
% K:  number of elements of the set of lags (square root of K must be an
%     integer)
% Pf: probability of false alarm
%
% Outputs:
% --------
% Txx:    test statistic
% Tgamma: test statistic after thresholding (Tgamma(Txx<Gamma)=NaN)
%
% Example of use:
% ---------------
% [Txx,Tgamma]=almost_cyclo_test(imresize(randn(128),1.5),128,11,4,1e-4);
%
%--------------------------------------------------------------------------
% This code is provided only for research purposes.
%--------------------------------------------------------------------------


% Check square root of K
if (mod(sqrt(K),1)==0)
    max_ind_tau=sqrt(K);
else
    disp('[error]: The square root of K must be an integer value.');
    Txx=0; Tgamma=0;
    return;
end

% Check size of spectral window
if (mod(P,2)==1)
    beta=1;
    W = sqrt(beta)*kaiser(P,beta)*kaiser(P,beta)';
else
    disp('[error]: The spectral window size P must be odd.');
    Txx=0; Tgamma=0;
    return;
end

% Check size of the image block according to the set of lags
if ((size(X,1)>=N+max_ind_tau)&&(size(X,2)>=N+max_ind_tau))
    x=X(1:N+max_ind_tau-1,1:N+max_ind_tau-1);
else
    disp('[error]: The size of block X is too small for the specified set of lags.');
    Txx=0; Tgamma=0;
    return;
end

% Obtain a zero mean image block
x=x-mean(mean(x));
% Rearrange image blocks according to the set of lags
xrow=im2col(x,[N N],'sliding');

% Preallocate memory for variables
F_tau=zeros(N,N,K); F_tau_s=zeros(N+P-1,N+P-1,K);

% Compute cyclic correlations
for tauk=1:K
    % Cyclic correlation (as vectors)
    cxx_tauk=xrow(:,1).*xrow(:,tauk);
    % 2D Fourier Transform of cyclic correlation
    F_tau(:,:,tauk)=fft2(reshape(cxx_tauk,[N N]),N,N);
    % Shifted version for windowing
    temp=repmat(F_tau(:,:,tauk),[3 3]);
    index=N-(P-1)/2+1;
    F_tau_s(:,:,tauk)=temp(index:index+(N-1)+(P-1),index:index+(N-1)+(P-1));
end

% Clear variables
clear xrow; clear x; clear temp;
% Preallocate memory for variables
S_kl=zeros(N,N,K*K); S_conj_kl=zeros(N,N,K*K);

% Compute the elements of the matrix Sigmaxx
for tauk=1:K
    for taul=1:K
        % Cyclic cross-spectrum estimation
        WFkFl=filter2(W,F_tau_s(:,:,tauk).*F_tau_s(:,:,taul),'full')/(P^2);
        S_kl(:,:,tauk+K*(taul-1))=WFkFl(P:end-(P-1),P:end-(P-1));
        % Cyclic cross-spectrum estimation
        WFkFl=filter2(W,F_tau_s(:,:,tauk).*conj(F_tau_s(:,:,taul)),'full')/(P^2);
        S_conj_kl(:,:,tauk+K*(taul-1))=WFkFl(P:end-(P-1),P:end-(P-1));
    end
end

% Clear variable
clear F_tau_s
% Preallocate memory for variables
Txx=zeros(N); Tgamma=NaN*ones(N);

% For each frequency pair (alpha1,alpha2) apply the algorithm
for alpha1=1:N
    for alpha2=1:N
        % Build up the vector of cyclic correlations
        cxx(:,1)=(1/sqrt(2))*[(reshape(F_tau(alpha1,alpha2,:),1,K)) conj(reshape(F_tau(alpha1,alpha2,:),1,K))];
        % Estimators of covariance matrix
        S_conj=reshape(S_conj_kl(alpha1,alpha2,:),K,K);
        S=reshape(S_kl(alpha1,alpha2,:),K,K);
        % Build up the estimated covariance matrix
        Sigmaxx=(1/2)*[S_conj, S; conj(S), conj(S_conj)];
        % Calculate the test statistic
        Txx(alpha1,alpha2)=cxx'*pinv(Sigmaxx)*cxx;
    end
end

% Set the threshold according to Pf
Gamma = chi2inv(1-Pf,2*K);
% Test statistic Txx after thresholding
Tgamma(abs(Txx)>=Gamma)=abs(Txx(abs(Txx)>=Gamma));

end