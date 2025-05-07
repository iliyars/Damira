clear 
clc

%% Параметры маховиков
m_rw = 0.13;       % Масса маховика [кг]
r_rw = 42/1000;    % Радиус маховика [м]
h_rw = 19/1000;    % Высота маховика [м]

%% Максимальные характеристики
rpm = 8000;                   % Максимальная скорость [об/мин]
maxSpeed = rpm*2*pi/60;       % Переводим в [рад/с]
maxTorque = 0.004;            % Максимальный крутящий момент [Н·м]

%% Ориентация маховиков в системе координат спутника
n1 = [1; 0; 0];  % Ось X спутника
n2 = [0; 1; 0];  % Ось Y спутника
n3 = [0; 0; 1];  % Ось Z спутника

%% Смещение маховиков от центра масс [м]
r1 = [4; 0; 0]/1000;
r2 = [0; 4; 0]/1000;
r3 = [0; 0; 4]/1000;

%% Расчет моментов инерции маховиков
% Момент инерции относительно поперечной оси (цилиндр)
ICylyndr = (1/12)*m_rw*(3*r_rw^2 + h_rw^2);
% Момент инерции относительно оси диска
IDisk = (1/2)*m_rw*r_rw^2;

% Тензор инерции в локальной системе координат маховика
I_rw_local = diag([IDisk, ICylyndr, ICylyndr]);

%% Преобразование в систему координат спутника
R1 = Rscrew(n1);
R2 = Rscrew(n2);
R3 = Rscrew(n3);

Ir1B = R1'*I_rw_local*R1;
Ir2B = R2'*I_rw_local*R2;
Ir3B = R3'*I_rw_local*R3;

%% Матрица управления моментами
J = [Ir1B*n1, Ir2B*n2, Ir3B*n3];
Jinv = inv(J'*J)*J';  % Псевдообратная матрица

%% Учет смещения маховиков (теорема Гюйгенса-Штейнера)
sr1 = skew(r1);
Ir1Bcg = Ir1B + m_rw*(sr1')*sr1;
sr2 = skew(r2);
Ir2Bcg = Ir2B + m_rw*(sr2')*sr2;
sr3 = skew(r3);
Ir3Bcg = Ir3B + m_rw*(sr3')*sr3;

%% Параметры для динамического моделирования
% Момент инерции спутника (включая маховики)
I_satellite = diag([0.1, 0.1, 0.1]) + Ir1Bcg + Ir2Bcg + Ir3Bcg;

% Инициализация состояний
omega_wheel = [0; 0; 0];   % Начальная угловая скорость маховиков

% Время моделирования
tspan = [0 60]; % 60 секунд

% Начальные условия: угловая скорость спутника и маховиков
initial_state = [0.1; 0.2; -0.15; omega_wheel];

% Решение ОДУ
[t, state] = ode45(@(t, state) satelliteDynamics(t, state, I_satellite, ...
                    m_rw, maxTorque, maxSpeed, Jinv), tspan, initial_state);

% Извлечение результатов
omega_sat = state(:,1:3);       % Угловая скорость спутника
omega_wheel = state(:,4:6);     % Угловая скорость маховиков

% Визуализация
figure;
subplot(2,1,1);
plot(t, omega_sat);
title('Угловая скорость спутника');
legend('X','Y','Z');
xlabel('Время (с)');
ylabel('рад/с');

subplot(2,1,2);
plot(t, omega_wheel);
title('Угловая скорость маховиков');
legend('X','Y','Z');
xlabel('Время (с)');
ylabel('рад/с');

%% Функция динамики системы
function dstate = satelliteDynamics(t, state, I_sat, m_rw, maxTorque, maxSpeed, Jinv)
    omega_sat = state(1:3);     % Угловая скорость спутника
    omega_wheel = state(4:6);   % Угловая скорость маховиков
    
    % Расчет управляющего момента (PD-регулятор)
    torque_desired = -0.5 * omega_sat - 0.1 * (omega_sat - [0;0;0]);
    
    % Ограничение момента
    torque_applied = sign(torque_desired) .* min(abs(torque_desired), maxTorque);
    
    % Расчет требуемых ускорений маховиков
    alpha_wheel = Jinv * torque_applied;
    
    % Динамика спутника (уравнение Эйлера)
    domega_sat = I_sat \ (-cross(omega_sat, I_sat*omega_sat) - torque_applied);
    
    % Динамика маховиков
    domega_wheel = alpha_wheel;
    
    % Ограничение скорости маховиков
    for i = 1:3
        if (omega_wheel(i) > maxSpeed && domega_wheel(i) > 0) || ...
           (omega_wheel(i) < -maxSpeed && domega_wheel(i) < 0)
            domega_wheel(i) = 0;
        end
    end
    
    dstate = [domega_sat; domega_wheel];
end

%% Вспомогательные функции
function R = Rscrew(nhat)
    % Матрица поворота для выравнивания оси z с вектором nhat
    x = nhat(1);
    y = nhat(2);
    z = nhat(3);
    
    % Углы Эйлера
    psi = atan2(y,x);
    theta = atan2(z,sqrt(x^2+y^2));
    phi = 0;
    
    R = TIB(phi,theta,psi);
end

function out = TIB(phi,theta,psi)
    % Матрица поворота по углам Эйлера (Z-Y-X)
    ct = cos(theta);
    st = sin(theta);
    sp = sin(phi);
    cp = cos(phi);
    ss = sin(psi);
    cs = cos(psi);
    
    out = [ct*cs, sp*st*cs-cp*ss, cp*st*cs+sp*ss;
           ct*ss, sp*st*ss+cp*cs, cp*st*ss-sp*cs;
           -st,    sp*ct,          cp*ct];
end

function out = skew(vec)
    % Кососимметрическая матрица для векторного произведения
    if length(vec) ~= 3
        error('Только 3x1 векторы поддерживаются');
    end
    
    out = [0,      -vec(3),  vec(2);
           vec(3),  0,      -vec(1);
          -vec(2),  vec(1),  0];
end