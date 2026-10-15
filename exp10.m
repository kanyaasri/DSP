
clc;
clear;
close all;

%% INPUTS
b = [0.0048 0.0193 0.0289 0.0193 0.0048];
a = [1 -2.3695 2.3140 -1.0547 0.1874];

Bc = 8;                 % Coefficient word length
Bp = 16;                % Product word length
N = 200;

x = [1 zeros(1,N-1)];   % Unit impulse

%% 1. COEFFICIENT QUANTIZATION
S = 2^(Bc-1)-1;

bq = round(b*S)/S;
aq = round(a*S)/S;

CoeffError = [b-bq a-aq];

%% 2. IDEAL AND COEFFICIENT-QUANTIZED RESPONSE
y_ideal = filter(b,a,x);
y_coeff = filter(bq,aq,x);

e_coeff = y_ideal-y_coeff;

MaxCoeffError = max(abs(e_coeff));
SQNR_coeff = 10*log10(sum(y_ideal.^2)/sum(e_coeff.^2));

%% 3. PRODUCT QUANTIZATION
Sp = 2^(Bp-1)-1;

y_prod = zeros(1,N);
v = zeros(1,4);

for n = 1:N

    p = bq*x(n);
    fb = -aq(2:5).*v;

    p = round([p fb]*Sp)/Sp;

    y_prod(n) = sum(p);

    v = [y_prod(n) v(1:3)];

end

e_prod = y_coeff-y_prod;

MaxProdError = max(abs(e_prod));
SQNR_prod = 10*log10(sum(y_coeff.^2)/sum(e_prod.^2));

%% 4. ZERO-INPUT LIMIT-CYCLE TEST
x0 = zeros(1,N);

state = [0.01 0.01 0.01 0.01];
y_zero = zeros(1,N);

for n = 1:N

    acc = -aq(2:5)*state.';

    y_zero(n) = round(acc*Sp)/Sp;

    state = [y_zero(n) state(1:3)];

end

LimitAmplitude = max(abs(y_zero));

LimitSamples = find(abs(y_zero)>1e-12);

%% 5. OVERFLOW / SATURATION TEST
xo = 1.2*ones(1,N);

y_over = filter(bq,aq,xo);

% Normalized saturation
y_sat = min(max(y_over,-1),1);

OverflowCount = sum(abs(y_over)>1);

MaxOverflow = max(abs(y_over));

%% 6. OUTPUT RESULTS
fprintf('\n========== QUANTIZATION RESULTS ==========\n');

fprintf('Maximum coefficient/output error = %.6e\n',MaxCoeffError);
fprintf('Coefficient SQNR = %.2f dB\n',SQNR_coeff);

fprintf('Maximum product quantization error = %.6e\n',MaxProdError);
fprintf('Product SQNR = %.2f dB\n',SQNR_prod);

fprintf('Zero-input maximum amplitude = %.6e\n',LimitAmplitude);
fprintf('Overflow samples = %d\n',OverflowCount);
fprintf('Maximum output magnitude = %.6f\n',MaxOverflow);

%% OUTPUT TABLE
Test = [
    "Coefficient Quantization"
    "Product Quantization"
    "Zero-Input Test"
    "Overflow Test"
];

Maximum_Error_or_Amplitude = [
    MaxCoeffError
    MaxProdError
    LimitAmplitude
    MaxOverflow
];

SQNR_dB = [
    SQNR_coeff
    SQNR_prod
    NaN
    NaN
];

T = table(Test,Maximum_Error_or_Amplitude,SQNR_dB);

disp(T);

%% 7. QUANTIZED COEFFICIENTS
disp('Original b =');
disp(b);

disp('Quantized b =');
disp(bq);

disp('Original a =');
disp(a);

disp('Quantized a =');
disp(aq);

%% 8. IDEAL VS COEFFICIENT-QUANTIZED RESPONSE
figure;

stem(0:N-1,y_ideal,'filled');
hold on;
stem(0:N-1,y_coeff);

grid on;
xlabel('n');
ylabel('Amplitude');
title('Ideal vs Coefficient-Quantized Response');
legend('Ideal','Quantized');

%% 9. COEFFICIENT QUANTIZATION NOISE
figure;

plot(0:N-1,e_coeff,'LineWidth',1);

grid on;
xlabel('n');
ylabel('Error');
title('Coefficient Quantization Noise');

%% 10. PRODUCT QUANTIZATION NOISE
figure;

stem(0:N-1,y_coeff-y_prod,'filled');

grid on;
xlabel('n');
ylabel('Error');
title('Product Quantization Noise');

%% 11. ZERO-INPUT LIMIT-CYCLE TEST
figure;

plot(0:N-1,y_zero,'LineWidth',1);

grid on;
xlabel('n');
ylabel('Amplitude');
title('Zero-Input Limit-Cycle Test');

%% 12. OVERFLOW / SATURATION TEST
figure;

plot(0:N-1,y_over);
hold on;
plot(0:N-1,y_sat,'--');

grid on;
xlabel('n');
ylabel('Amplitude');
title('Overflow / Saturation Test');
legend('Without Saturation','With Saturation');

%% 13. IDEAL FILTER POLE-ZERO PLOT
figure;

zplane(b,a);
grid on;
title('Ideal Filter Pole-Zero Plot');

%% 14. QUANTIZED FILTER POLE-ZERO PLOT
figure;

zplane(bq,aq);
grid on;
title('Quantized Filter Pole-Zero Plot');

%% 15. FINAL RESULT
fprintf('\n========== FINAL RESULT ==========\n');

fprintf('Coefficient quantization effects were analyzed.\n');
fprintf('Product quantization effects were analyzed.\n');
fprintf('Zero-input limit-cycle behavior was tested.\n');
fprintf('Overflow and saturation effects were analyzed.\n');
fprintf('Pole-zero plots were obtained for both filters.\n');
