format long;
D = 1;
L = pi;

errores = zeros(6,1);
deltas = zeros(6,1);
cont_err = 0;
for M = [15,30]
 for N=[5,10,15]    
    cont_err = cont_err+1;
     T = 1; % tambien en 1
     %T=1;
    %discretizacion espacial
    delta_y = L/(M + 1);
    y = linspace(0, L, M+1);
    x = y;
    deltas(cont_err) = delta_y;

    K_N = zeros(M+1,M+1);
    for i = 1:M+1
        for j = 1:M+1
        a_n = exp((-((1:N)*pi).^2*D*T)./L^2).*sin((1:N).*pi.*x(i)/L).*sin((1:N).*pi.*y(j)/L);
        K_N(i, j) = (2/L)*sum(a_n);
        end
    end
    K_N(:,1) = (1/2)*K_N(:,1);
    K_N(:,end) = (1/2)*K_N(:,end);
    K_N = (L/M)*K_N;

    f0 =2/pi*x.*(x<=pi/2)+2/pi*(pi-x).*(x>pi/2);

    u_spectral = K_N*f0';
    %disp(gf);
    %disp(delta_y);
    if M==30 && N==10
        A_Kn = K_N;
    end

    % Comparar con la solución del proyecto 2
    M_p2 = 40;  
    N_p2 = 1000*T+1;  

    %discretizacion espacial
    delta_x_p2 = L/(M_p2 + 1);
    x_p2 = linspace(0, L, M_p2+2);
    %discretizacion del tiempo
    delta_t_p2 = T/(N_p2+1);
    t_p2 = linspace(0,T,N_p2+1);
    d_p2 = D*delta_t_p2/(delta_x_p2)^2;
    Aprin = sparse(1:M_p2, 1:M_p2, (1-2*d_p2), M_p2, M_p2);
    Asup = sparse(1:M_p2-1, 2:M_p2, d_p2, M_p2, M_p2); 
    Ainf = sparse(2:M_p2, 1:M_p2-1, d_p2,M_p2, M_p2);
    A = Aprin + Asup + Ainf;

    f0_p2 = 2/pi*x_p2.*(x_p2<=pi/2)+2/pi*(pi-x_p2).*(x_p2>pi/2);
    u_p2 = zeros(M_p2,N_p2+2);
    u_p2(:,1) = f0_p2(2:M_p2+1);
    for j = 2:N_p2+2 
        u_p2(:,j) = A * u_p2(:,j-1);
    end

    u_diff = spline(x_p2(2:end-1),u_p2(:,end),x);

    errores(cont_err) = norm(u_diff' - u_spectral)/norm(u_spectral);
  end
end
M = [15;15;15;30;30;30];
N = [5;10;15;5;10;15];
Tab = table(M,N,deltas,errores)
%Tab2 =table(M,N,deltas,errores)
%writetable(Tab,'tabla_proy3_T0.1.csv')
%writetable(Tab2,'tabla_proy3_T1.csv')

figure;
plot(x,u_spectral);
hold on
plot(x, f0, 'k--', 'LineWidth', 1);
plot(x,u_diff);
grid on
xlabel('x');
ylabel('Temp');
legend('Solución espectral','Original Exacta', 'Solución Dif. Finitas','Location','best');
title(['Comparacion de soluciones para T=', num2str(T)]);
hold off

% EJERCICIO 2.1
% usar A_Kn para el ejercicio 2.2
%G = A_Kn*f0';
A_inv = A_Kn(2:end-1, 2:end-1);
G_inv = u_diff(2:end-1)';
condi = cond(A_inv);
fprintf('el condicionamiento de la matriz del sistema sin condiciones de contorno es: %e\n', condi);

% Matriz de Tikhonov
phi0 = 0;
epsilon = 0.01;
alphas = logspace(-10, -1, 200);
n = length(G_inv);
lim_inf = epsilon * norm(G_inv);
lim_sup = 1.2*lim_inf;
soluc_found = false;
normas = zeros(size(alphas));

for k= 1:length(alphas)
    alpha = alphas(k);
    A_tk = A_inv'*A_inv + alpha*eye(size(A_inv));
    b_tk = A_inv'*G_inv + alpha*phi0;
    d_alpha = A_tk\b_tk;
    
    residual = norm(A_inv*d_alpha - G_inv);
    normas(k) = residual;
    if residual <= lim_sup && residual >= lim_inf 
        soluc_found = true;
        soluc_final = d_alpha;
        alpha_optk=alpha;
        break;
    end
end

if soluc_found
    disp(['Alpha óptimo encontrado: ', num2str(alpha)]);
    figure;
    plot(x(2:end-1),soluc_final);
    hold on
    plot(x, f0, 'k--', 'LineWidth', 1);
    grid on
    xlabel('x');
    ylabel('Temp');
    legend('Reconstrucción Tikhonov', 'Original Exacta','Location','best');
    title(['Solucion para alpha=',num2str(alpha)]);
    hold off;
    Condi_tk = cond(A_tk);
    fprintf('el condicionamiento de la matriz para el método de Tikhonov es: %.10f\n', Condi_tk);
else
    disp('No se encontró un alpha que cumpla el criterio en este rango.');
end

%Tikhonov 1:
phi0 = soluc_final; 
b_tk1 = A_inv'*G_inv + alpha*phi0;
d_tk1 = A_tk\b_tk1;
%Tikhonov 2:
phi0 = d_tk1;
b_tk2 = A_inv'*G_inv + alpha*phi0;
f_tikhonov2 = A_tk\b_tk2;
%Tikhonov 3:
phi0 = f_tikhonov2;
b_tk2 = A_inv'*G_inv + alpha*phi0;
f_tikhonov3 = A_tk\b_tk2;

figure;
plot(x(2:end-1),soluc_final);
hold on
plot(x(2:end-1), d_tk1);
plot(x(2:end-1),f_tikhonov2);
plot(x(2:end-1),f_tikhonov3);
plot(x, f0, 'k--', 'LineWidth', 1);
grid on
xlabel('x');
ylabel('Temp');
legend('Reconstrucción Tikhonov', 'Reconstrucción Tikhonov 1', 'Reconstrucción Tikhonov 2', 'Reconstrucción Tikhonov 3', 'Original Exacta','Location','best');
title(['Solucion para alpha=',num2str(alpha)]);
hold off

miu_v = 1/norm(A_inv'*A_inv);
num_iter = [10,500,10000];
coef_miu = 0.95;
figure;
plot(x, f0, 'k--', 'LineWidth', 1,'DisplayName','Original Exacta');
hold on
colores = {'r','b','g'};
for i=1:length(num_iter)
    K_max = num_iter(i);
    f_k = zeros(length(G_inv),1);
    Miu = miu_v*coef_miu;
    for k=1:K_max
        f_k = f_k - Miu*A_inv'*(A_inv*f_k - G_inv);
    end
    plot(x(2:end-1), f_k,colores{i}, 'LineWidth',1.2,'DisplayName',['Landweber:',num2str(K_max),'iters']);
end
grid on
xlabel('x');
ylabel('Temp');
legend('Location','best');
title(['Solucion para coeficiente de miu=',num2str(coef_miu)]);
hold off

ruido = 0.01*max(G_inv)*(2*rand(size(G_inv)) - ones(size(G_inv)));
G_noisy = G_inv + ruido;

phi0 = 0;
alphas = logspace(-10, -1, 200);
n = length(G_noisy);
lim_inf = epsilon * norm(G_noisy);
lim_sup = 1.5*lim_inf;
soluc_found = false;
normas = zeros(size(alphas));

for k= 1:length(alphas)
    alpha = alphas(k);
    

    A_tk = A_inv'*A_inv + alpha*eye(size(A_inv));
    b_tk = A_inv'*G_noisy + alpha*phi0;
    d_alpha = A_tk\b_tk;
    
    residual = norm(A_inv*d_alpha - G_noisy);
    normas(k) = residual;
    if residual <= lim_sup && residual >= lim_inf
        soluc_found = true;
         soluc_final_noisy = d_alpha;
        break;
    end
end

if soluc_found
    disp(['Alpha óptimo encontrado: ', num2str(alpha)]);
    figure;
    plot(x(2:end-1),soluc_final_noisy);
    hold on
    plot(x, f0, 'k--', 'LineWidth', 1);
    grid on
    xlabel('x');
    ylabel('Temp');
    legend('Recon Tikhonov Datos Ruidosos', 'Original Exacta','Location','best');
    title(['Solucion para alpha=',num2str(alpha)]);
    hold off;
    Condi_tk = cond(A_tk);
    fprintf('el condicionamiento de la matriz para el método de Tikhonov es: %.10f\n', Condi_tk);
else
    disp('No se encontró un alpha que cumpla el criterio en este rango.');
end

%Tikhonov 1:
phi0 = soluc_final_noisy; 
b_tk1 = A_inv'*G_noisy + alpha*phi0;
d_tk1 = A_tk\b_tk1;
%Tikhonov 2:
phi0 = d_tk1;
b_tk2 = A_inv'*G_noisy + alpha*phi0;
f_tikhonov2 = A_tk\b_tk2;

figure;
plot(x(2:end-1),soluc_final_noisy);
hold on
plot(x(2:end-1), d_tk1);
plot(x(2:end-1),f_tikhonov2);
plot(x, f0, 'k--', 'LineWidth', 1);
grid on
xlabel('x');
ylabel('Temp');
legend('Recons Tikhonov Ruidosos', 'Recons Tikhonov 1 Ruidosos', 'Recons Tikhonov 2 Ruidosos', 'Original Exacta','Location','best');
title(['Solucion para alpha=',num2str(alpha)]);
hold off

%Miu = 0.95/norm(A_inv'*A_inv);
num_iter = [10,500,10000];
coef_miu = 0.95;
figure;
plot(x, f0, 'k--', 'LineWidth', 1,'DisplayName','Original Exacta');
hold on
colores = {'r','b','g'};
for i=1:length(num_iter)
    K_max = num_iter(i);
    f_k = zeros(length(G_noisy),1);
    Miu = miu_v*coef_miu;
    for k=1:K_max
        f_k = f_k - Miu*A_inv'*(A_inv*f_k - G_noisy);
    end
    plot(x(2:end-1), f_k,colores{i}, 'LineWidth',1.2,'DisplayName',['Landweber:',num2str(K_max),'iters']);
end
grid on
xlabel('x');
ylabel('Temp');
legend('Location','best');
title(['Solucion para coeficiente de miu=',num2str(coef_miu)]);
hold off

%A_inv*F=u_spectral
r = rank(A_inv);
%G = u_spectral(2:end-1,:);
[U,S,V] = svd(A_inv);
sigma1 = S(1,1);
sing_val = diag(S);
log_sing_val = (1/sigma1)*sing_val;

figure;
semilogy(1:length(log_sing_val), log_sing_val, 'x', 'LineWidth', 1);
grid on;
xlabel('Index');
ylabel('Singular Values');
title('Valores singulares en escala logaritmica');


m_op = sum(log_sing_val > 1e-16);
fprintf('La cantidad de modos optimos es: %.10f\n', m_op);
fprintf('El rango de A es: %.10f\n', r);

for k=1:r
    figure;
    plot(x, f0, 'k--', 'LineWidth', 1,'DisplayName','Original Exacta');
    hold on
    F_tsvd = zeros(size(V, 1), 1);
    F_wsvd = zeros(size(V, 1), 1);
    for j=1:k
        c(j) = (1/sing_val(j))*U(:,j)'*G_inv;
        F_tsvd = F_tsvd + c(j)*V(:,j);
    end
    for j=1:r
        if j <= k
            c_wsvd = (1/sing_val(j))*U(:,j)'*G_inv;
        else
            c_wsvd = (1/sing_val(k))*U(:,j)'*G_inv;
        end
        F_wsvd = F_wsvd + c_wsvd*V(:,j);
    end
    plot(x(2:end-1), F_tsvd,'b', 'LineWidth',1.2,'DisplayName',['TSVD para k=',num2str(k)]);
    plot(x(2:end-1), F_wsvd,'g', 'LineWidth',1.2,'DisplayName',['WSVD para k=',num2str(k)]);
    grid on
    xlabel('x');
    ylabel('Temp');
    legend('Location','best');
    hold off;
end
ruido2 = 0.01*max(G_inv)*(2*rand(size(G_inv)) - ones(size(G_inv)));
G_ruido = G_inv + ruido2;

for k=1:r
    figure;
    plot(x, f0, 'k--', 'LineWidth', 1,'DisplayName','Original Exacta');
    hold on
    F_tsvd = zeros(size(V, 1), 1);
    F_wsvd = zeros(size(V, 1), 1);
    for j=1:k
        cj = (1/sing_val(j))*U(:,j)'*G_ruido;
        F_tsvd = F_tsvd + cj*V(:,j);
    end
    for j=1:r
        if j <= k
            c_wsvd = (1/sing_val(j))*U(:,j)'*G_ruido;
        else
            c_wsvd = (1/sing_val(k))*U(:,j)'*G_ruido;
        end
        F_wsvd = F_wsvd + c_wsvd * V(:,j);
    end
    plot(x(2:end-1), F_tsvd,'b', 'LineWidth',1.2,'DisplayName',['TSVD con ruido para k=',num2str(k)]);
    plot(x(2:end-1), F_wsvd,'g', 'LineWidth',1.2,'DisplayName',['WSVD con ruido para k=',num2str(k)]);
    grid on
    xlabel('x');
    ylabel('Temp');
    legend('Location','best');
    hold off;
end
