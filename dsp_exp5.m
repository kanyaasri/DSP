
clc;
clear;
close all;

%% IIR FILTER COEFFICIENTS
b = [0.0048 0.0193 0.0289 0.0193 0.0048];
a = [1 -2.3695 2.3140 -1.0547 0.1874];

N = 50;
x = [1 zeros(1,N-1)];
n = 0:N-1;

%% 1. DIRECT FORM-I
y1 = zeros(1,N);

for i = 1:N
    for j = 1:length(b)
        if i-j+1 > 0
            y1(i) = y1(i) + b(j)*x(i-j+1);
        end
    end

    for j = 2:length(a)
        if i-j+1 > 0
            y1(i) = y1(i) - a(j)*y1(i-j+1);
        end
    end

    y1(i) = y1(i)/a(1);
end

%% 2. DIRECT FORM-II (TRANSPOSED)
y2 = filter(b,a,x);

%% 3. CASCADE REALIZATION
[sos,g] = tf2sos(b,a);

yc = g*x;

for i = 1:size(sos,1)
    yc = filter(sos(i,1:3),sos(i,4:6),yc);
end

%% 4. PARALLEL REALIZATION
[r,p,k] = residuez(b,a);

yp = filter(k,1,x);

for i = 1:length(r)
    yp = yp + filter(r(i),[1 -p(i)],x);
end

%% 5. DISPLAY COEFFICIENTS
disp('Numerator coefficients:');
disp(b);

disp('Denominator coefficients:');
disp(a);

disp('Second Order Sections:');
disp(sos);

disp('Cascade Gain:');
disp(g);

disp('Parallel Residues:');
disp(r);

disp('Parallel Poles:');
disp(p);

disp('Parallel Direct Term:');
disp(k);

%% 6. ERROR AND DEVIATION
E = [max(abs(y1-y1));
     max(abs(y1-y2));
     max(abs(y1-yc));
     max(abs(y1-yp))];

reference = max(abs(y1));

if reference > 0
    D = E/reference*100;
else
    D = zeros(size(E));
end

disp(' ');
disp('COMPARISON TABLE');
disp('Structure          Maximum Error       Deviation (%)');

fprintf('Direct Form-I      %12.6g       %12.6g\n',E(1),D(1));
fprintf('Direct Form-II     %12.6g       %12.6g\n',E(2),D(2));
fprintf('Cascade            %12.6g       %12.6g\n',E(3),D(3));
fprintf('Parallel           %12.6g       %12.6g\n',E(4),D(4));

%% 7. IMPULSE RESPONSE COMPARISON
figure;

subplot(2,1,1);
stem(n,y1,'filled');
hold on;
stem(n,y2,'--');
stem(n,yc,':','LineWidth',1.2);
stem(n,yp,'-.');

grid on;
xlabel('Sample Index (n)');
ylabel('Amplitude');
title('Impulse Response Comparison');

legend('Direct Form-I','Direct Form-II', ...
       'Cascade','Parallel','Location','best');

%% 8. ERROR PLOT
subplot(2,1,2);

plot(n,abs(y1-y2),'--','LineWidth',1.2);
hold on;
plot(n,abs(y1-yc),':','LineWidth',1.2);
plot(n,abs(y1-yp),'-.','LineWidth',1.2);

grid on;
xlabel('Sample Index (n)');
ylabel('Absolute Error');
title('Numerical Error Comparison');

legend('DF-I vs DF-II','DF-I vs Cascade', ...
       'DF-I vs Parallel','Location','best');

%% 9. FREQUENCY RESPONSE
figure;
freqz(b,a,1024);
title('Frequency Response of IIR Filter');

%% 10. POLE-ZERO PLOT
figure;
zplane(b,a);
grid on;
title('Pole-Zero Plot of IIR Filter');
