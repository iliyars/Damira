% Параметры маховиков (примерные значения)
mw.mass = 0.5;     % Масса маховика, кг
mw.radius = 0.05;   % Радиус маховика, м
mw.I = 0.5 * mw.mass * mw.radius^2; % Момент инерции (диск)

% Ограничения маховиков
mw.max_torque = 0.1;   % Максимальный момент, Н·м
mw.max_speed = 6000;   % Максимальная угловая скорость, об/мин -> рад/с
mw.max_speed_rad = mw.max_speed * 2*pi/60;

% Инициализация состояний
omega_wheel = [0; 0; 0];   % Угловая скорость маховиков (x,y,z), рад/с
I_satellite = diag([0.1, 0.1, 0.1]); % Момент инерции спутника, кг·м²

% Время моделирования
tspan = [0 60]; % 60 секунд
dt = 0.1;       % Шаг интегрирования

% Уравнения движения спутника с маховиками
odefun = @(t, state) satelliteDynamics(t, state, I_satellite, mw);

% Начальные условия: угловая скорость спутника и маховиков
initial_state = [0.1; 0.2; -0.15; omega_wheel];

% Решение ОДУ
[t, state] = ode45(odefun, tspan, initial_state);

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

% Функция динамики системы
function dstate = satelliteDynamics(t, state, I_sat, mw)
    omega_sat = state(1:3);     % Угловая скорость спутника
    omega_wheel = state(4:6);   % Угловая скорость маховиков
    
    % Расчет управляющего момента (пример: стабилизация)
    torque_desired = -0.5 * omega_sat; % PD-регулятор
    
    % Ограничение момента и скорости маховиков
    torque_applied = sign(torque_desired) .* ...
        min(abs(torque_desired), mw.max_torque);
    
    % Динамика спутника (уравнение Эйлера)
    domega_sat = I_sat \ ( -cross(omega_sat, I_sat*omega_sat) - torque_applied );
    
    % Динамика маховиков (уравнение вращения)
    domega_wheel = torque_applied / mw.I;
    
    % Ограничение скорости маховиков
    domega_wheel = sign(domega_wheel) .* ...
        min(abs(domega_wheel), mw.max_speed_rad - abs(omega_wheel));
    
    dstate = [domega_sat; domega_wheel];
end