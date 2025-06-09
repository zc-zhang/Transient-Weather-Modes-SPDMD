function [ KModes,KEv,KAmps,Norms,XCdmd ] = CompanionMatrix_DMD( Data )

% References:
% Dynamic Mode Decomposition as presented by 
% "Spectral analysis of nonlinear flows" by Rowley et al., Journal of FLuid
% Mechanics, 2009
% "Ergodic theory, dynamic mode decomposition, and computation of spectral properties of the Koopman operator"
%  by Hassan Arbabi and Ignor Mezic, SIADS'17

% inputs : 
% Data - the data matrix: each column is a a set of measurements done at
% each instant - the sampling frequency is assumed to be constant

% outputs: 
% 1 - KModes- Koopman or Dynamic Modes: structures that scale exponentially
% with time - the exponent is given by Koopman or Dynamic Eigenvalues and
% could be imaginary as well

% 2- KEv - Koopman or Dynamic Eigenvalues: the exponents and frequencies
% that make up the time-varaiation of data in time

% 3- Norms - Euclidean (vector) norm of each mode, used to sort the data

% 4- XCdmd - Reconstruction or superposition of the Companion DMD or KMD.




X=Data(:,1:end-1);

c = pinv(X)*Data(:,end);

m = size(Data,2)-1;
C = spdiags(ones(m,1),-1,m,m);
C(:,end)=c;                  %% companion matrix

[P,D]=eig(full(C));                           %% Ritz evalue and evectors

KEunsrtd = diag(D);                     %% Koopman Evalues
KMunsrtd = X*P;                         %% Koopman Modes

 Cnorm = sqrt(sum(abs(KMunsrtd).^2,1));  %% Euclidean norm of Koopman Modes
 [Norms, ind]=sort(Cnorm,'descend');      %% sorting based on the norm KM
 %[Norms, ind]=sort(real(KEunsrtd),'descend');          %% sorting based on real part of Koopman eigenvalues
 KEv= KEunsrtd(ind);           %% Koopman Eigenvalues

KModes = KMunsrtd(:,ind);                %% Koopman Modes

% Compute modal amplitudes (initial projection onto Koopman modes)
 b = KModes\X(:,1);  % Least-squares solution 
% b = pinv(KModes) * X(:,1);   % or consider (pseudo-inverse of KModes)
 KAmps=b; % rename amplitudes

% Time evolution (you can modify the time steps `t` as needed)
t = 0:(size(Data,2)-2);  %% Time steps corresponding to the snapshots
Vandermonde = KEv(:).^t;  %% Vandermonde matrix for time evolution

% Reconstruct the data using the Koopman modes and time evolution
XCdmd = KModes * (Vandermonde .* b);  % Element-wise multiply by amplitudes
%XCdmd = KModes * diag(b)*Vandermonde;  % Element-wise multiply by amplitudes 
end



% UC Santa Barbara
% arbabiha@gmail.com
%=========================================================================%
