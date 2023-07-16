            %% Section 1
    %% Part 1 - Alef)
M = 2;
N = 2;
fs = 20;
WL = 1;
c = [0.2 0.4 0.6 -0.1 -0.3];
d = [0.1 0.3 -0.2 0.5 -0.3];
k = [1 2 3 4 5];
K = length(k);
T = WL * fs * K;
S = zeros(N,T);
for i = 1:length(k)
    t_temp = k(i)-1:WL/fs:k(i)-WL/fs;
    s1_temp = c(i) * sin(2*pi*t_temp);
    s2_temp = d(i) * sin(4*pi*t_temp);
    s_temp = [s1_temp;s2_temp];
    S(:,T/K*(i-1)+1:T/K*i) = s_temp;
end
% % ***
% S(1,:) = S(1,:)/norm(S(1,:));
% S(2,:) = S(2,:)/norm(S(2,:));
% % ***
A = [0.8 -0.6;0.6 0.8];
X = A*S;
figure(1);
t = 0:WL/fs:5-WL/fs;
subplot(4,1,1);
plot(t,S(1,:),'b');
xlabel('t');
title("s_1(t)");
grid on;
subplot(4,1,2);
plot(t,S(2,:),'b');
xlabel('t');
title("s_2(t)");
grid on;
subplot(4,1,3);
plot(t,X(1,:),'r');
xlabel('t');
title("x_1(t)");
grid on;
subplot(4,1,4);
plot(t,X(2,:),'r');
xlabel('t');
title("x_2(t)");
grid on;

    %% Part 2 - be)
X1 = X(:,1:20);
X2 = X(:,21:40);
Rx1 = X1 * transpose(X1);
Rx2 = X2 * transpose(X2);
[Q,D] = eig(inv(Rx2)*Rx1);
B = transpose(Q);
S_pred = B * X;
E = norm(S_pred-S,"fro")^2/norm(S,"fro")^2;
% Rx = X * transpose(X);
% [U,D] = eig(Rx);
% Z = D^(-1/2) * transpose(U) * X;
% Z1 = Z(:,1:20);
% Z2 = Z(:,21:40);
% Rz1 = Z1 * transpose(Z1);
% Rz2 = Z2 * transpose(Z2);
% [Q,D] = eig(inv(Rz2)*Rz1);
% B = transpose(Q);
% S_pred = B * Z;
% E = norm(S_pred-S,"fro")^2/norm(S,"fro")^2;

    %% Part 3 - jim)
B_init = [0.7 0.8;-0.5 0.9];
B_JD = JD(B_init,X,T,K,K,M,N);
S_pred_JD = B_JD * X;
E_JD = norm(S_pred_JD-S,"fro")^2/norm(S,"fro")^2;

    %% Part 4 - dal)
W = randn(M,T);
W = W/norm(W,"fro");
SNR = 20;
sigma = sqrt(norm(X,"fro")/10^(SNR/10));
Y = X + sigma*W;
B_y = JD(B_init,Y,T,K,K,M,N);
S_pred_y = B_y * Y;
E_y = norm(S_pred_y-S,"fro")^2/norm(S,"fro")^2;

    %% Part 5 - he)
SNR = 20;
sigma = sqrt(norm(X,"fro")/10^(SNR/10));
E_avg = zeros(4,1);
n_trials = 100;
K_values = [2,3,4,5];
for k = 2:K
    E_k = zeros(n_trials,1);
    for n=1:n_trials
        W = randn(M,T);
        W = W/norm(W,"fro");
        Y = X + sigma*W;
        B_y = JD(B_init,Y,T,K,k,M,N);
        S_pred_n = B_y * Y;
        E_n = norm(S_pred_n-S,"fro")^2/norm(S,"fro")^2;
        E_k(n) = E_n;
    end
    E_avg(k-1) = mean(E_k);
end
figure(2);
plot(K_values,E_avg,'o-');
xlabel('K');
ylabel("Error");
title("Source Signal Estimation Error");
grid on;

    %% Part 6 - ye)
SNR_Vals = [5,10,15,20];
E_avgs = zeros(length(SNR_Vals),1);
for i = 1:length(SNR_Vals)
   snr = SNR_Vals(i);
   sigma = sqrt(norm(X,"fro")/10^(snr/10));
   E_k = zeros(n_trials,1);
   for n=1:n_trials
        W = randn(M,T);
        W = W/norm(W,"fro");
        Y = X + sigma*W;
        B_y = JD(B_init,Y,T,K,k,M,N);
        S_pred_n = B_y * Y;
        E_n = norm(S_pred_n-S,"fro")^2/norm(S,"fro")^2;
        E_k(n) = E_n;
   end
   E_avgs(i) = mean(E_k);
end
figure(3);
plot(SNR_Vals,E_avgs,'o-');
xlabel('SNR [dB]');
ylabel("Error");
title("Source Signal Estimation Error K=5");
grid on;

    %% Local Necessary Functions
function B = JD(B_init,X,T,K_r,K,M,N)
    trials = 20;
    B = B_init;
    for n = 1:trials
        for i = 1:N
            Ri = Ri_Mat(X,K_r,K,T,B,M,N,i);
            [U,D] = eig(Ri);
            [~,idx] = min(min(D));
            b = U(:,idx);
            B(i,:) = b;
        end
    end
    B = -B;
end

function R = Ri_Mat(X,K_r,K,T,B,M,N,i)
    R = zeros(M,M);
    for k = 1 : K
        idx1 = floor(T/K_r*(k-1)+1);
        idx2 = floor(T/K_r*k);
        X_w = X(idx1:idx2);
        Rxk = X_w*transpose(X_w);
        for j = 1 : N
           if j == i
              continue; 
           end
           R = R + (Rxk * transpose(B(j,:))) * transpose((Rxk * transpose(B(j,:))));
        end
    end
end