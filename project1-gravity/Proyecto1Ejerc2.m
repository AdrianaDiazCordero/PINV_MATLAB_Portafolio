
% EJERCICIO 2, SEGUNDO SET DE DATOS
disp('----------------------')
fprintf('Ejercicio 2 \n');
disp('----------------------')
F_z2 = importdata("CampoGravitacional2.mat");
rng(42);
t2 = F_z2(1,:);
gravF2 = F_z2(2,:)';
n2 = numel(t2);
ks = 1:5;
res_sol = zeros(size(ks));
sol= zeros(size(ks));
res_sol_ruido = zeros(size(ks));
sol_ruido = zeros(size(ks));
eps = 0.05* range(gravF2);

%perturbando los datos de gravF
ruido2 = (2*rand(size(gravF2)) - 1) * eps;
gravF_ruido2 = gravF2 + ruido2;
x_plot = linspace (0,1,100);

fcomp = figure('Name', 'Comparacion de datos sin y con ruido 2do set de datos');

%for k=[2,4,5,7,10]
for m=1:length(ks)
    k=m;
    A2 = zeros(n2, k+1);
    

    for j=1:(k+1)
        for i=1:n2
            f2 = @(x) x.^(j-1)./((x-t2(i)).^2+1).^(3/2);
            A2(i,j) = integral(f2,0,1);
        end
    end
%calculos de coeficientes sin ruido
coefs2 = (A2.'*A2)\(A2.'*gravF2);
p_inv2 = flip(coefs2);
Lamda_pred2 = polyval(p_inv2,x_plot) ;
res_sol(m) = norm(A2*coefs2-gravF2);
sol(m) = norm(coefs2);

%calculo de coeficientes con ruido
coefs_ruido2 = (A2.'*A2)\(A2.'*gravF_ruido2);
p_invruido2 = flip(coefs_ruido2);
Lamda_predruido2 = polyval(p_invruido2,x_plot) ;
res_sol_ruido(m) = norm(A2*coefs_ruido2-gravF_ruido2);
sol_ruido(m)=norm(coefs_ruido2);

% Gráfica 1: Datos sin perturbar 
subplot(2, 1, 1); % 2 filas, 1 columna, posición 1
plot(x_plot, Lamda_pred2, 'DisplayName', sprintf('masa para k=%d',k))
hold on

% Gráfica 2: Datos perturbados
subplot(2, 1, 2); % 2 filas, 1 columna, posición 2
plot(x_plot, Lamda_predruido2, 'DisplayName', sprintf('masa para k=%d',k))
hold on

% Calculo de otros datos
residuo2 = norm(A2*coefs2 - gravF2);
fprintf('Para k = %d, el residuo es: %e\n', k, residuo2);
residuo_ruido2 = norm(A2*coefs_ruido2 - gravF_ruido2);
fprintf('Para k = %d, el residuo con ruido es: %e\n', k, residuo_ruido2);

%calculo del condicionamiento de A
condici2 = cond(A2);
fprintf('Para k = %d, el número de condición es: %e\n', k, condici2);
fprintf('Los coeficientes encontrados son: \n');
disp(coefs2)
fprintf('Los coeficientes encontrados con ruido son: \n');
disp(coefs_ruido2)
end

subplot(2, 1, 1);
title('Reconstrucción: Datos Originales 2do set de datos');
ylabel('\lambda(x)');
legend('Location', 'best');
grid on;

subplot(2, 1, 2);
title(['Reconstrucción: Datos con Ruido (\epsilon = ', num2str(eps), ')']);
ylabel('\lambda(x)');
legend('Location', 'best');
grid on;
hold off;

% Graficar la Curva L
figure('Name', 'Curva L (datos con ruido) 2do set de datos');
loglog(res_sol_ruido, sol_ruido, '-s', 'LineWidth', 2);
%loglog(norma_res_sinr, norma_sol_sinr, '-o', 'LineWidth', 2);
text(res_sol_ruido, sol_ruido, arrayfun(@(x) num2str(x), ks, 'UniformOutput', false), 'VerticalAlignment','bottom');
xlabel('Norma del Residuo');
ylabel('Norma de la Solución');
title('Método de la Curva L (Escala Logarítmica)');
grid on;

% Graficar la Curva L datos sin ruido
figure('Name', 'Curva L (datos sin ruido) 2do set de datos');
loglog(res_sol, sol, '-o', 'LineWidth', 2);
text(res_sol, sol, arrayfun(@(x) num2str(x), ks, 'UniformOutput', false), 'VerticalAlignment','bottom');
xlabel('Norma del Residuo');
ylabel('Norma de la Solución');
title('Método de la Curva L (datos sin ruido)');
grid on;

disp('------*******************----------')
k_codo = calc_codoCurvaL(res_sol,sol);
fprintf('El k óptimo según la curva L (automático) para datos sin ruido es es: %d\n', k_codo);
disp('------*******************----------')
disp('------*******************----------')
k_codo_ruido = calc_codoCurvaL(res_sol_ruido,sol_ruido);
fprintf('El k óptimo según la curva L (automático) para datos sin ruido es es: %d\n', k_codo_ruido);
disp('------*******************----------')