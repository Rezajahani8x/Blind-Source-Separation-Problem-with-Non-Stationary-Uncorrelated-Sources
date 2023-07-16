                        %% Section 2 - X1 Analysis
X1 = load("hw4-X1.mat").X1;

            %% Part 1 - alef)
    %% Plotting the observations
[~,T] = size(X1);
fs = 100;
figure(1);
subplot(2,1,1);
plot(1:T,X1(1,:));
xlabel('t');
title("x_{11}(t)");
grid on;
subplot(2,1,2);
plot(1:T,X1(2,:));
xlabel('t');
title("x_{12}(t)");
grid on;

    %% JD with 2 Windows
Rx = X1*transpose(X1);
[U,D] = eig(Rx);
Z1 = D^(-1/2) * transpose(U) * X1;
Zw1 = Z1(:,1:50);
Zw2 = Z1(:,25:74);
Zw3 = Z1(:,35:84);
Rz1 = Zw1 * transpose(Zw2);
Rz2 = Zw1 * transpose(Zw3);
[Q,D1] = eig(inv(Rz2)*Rz1);
B = transpose(Q);
S = B*Z1;
s1 = S(1,:);
s2 = S(2,:);
figure(2);
subplot(2,1,1);
plot(1:T,s1);
xlabel('t');
title('s_1(t)');
grid on;
subplot(2,1,2);
plot(1:T,s2);
xlabel('t');
title('s_2(t)');
grid on;

            %% Part 2 - be)
    %% Fourier Transform Analysis of Sources
S1 = fftshift(fft(s1));
S2 = fftshift(fft(s2));
X_1 = fftshift(fft(X1(1,:)));
X_2 = fftshift(fft(X1(2,:)));
f = -fs/2:fs/T:fs/2 - fs/T;
figure(3);
subplot(2,2,1);
plot(f,abs(S1));
xlabel("f(Hz)");
title("S_1(f)");
grid on;
subplot(2,2,2);
plot(f,abs(S2));
xlabel("f(Hz)");
title("S_2(f)");
grid on;
subplot(2,2,3);
plot(f,abs(X_1));
xlabel("f(Hz)");
title("X_1(f)");
grid on;
subplot(2,2,4);
plot(f,abs(X_2));
xlabel("f(Hz)");
title("X_2(f)");
grid on;