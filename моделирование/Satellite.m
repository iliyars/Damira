function dstatedt = Satellite(t,state)

global BB invI I m mu nextMagUpdate lastMagUpdate lastSensorUpdate maxSpeed
global nextSensorUpdate BfieldMeasured pqrMeasured BfieldNav pqrNav
global BfieldNavPrev pqrNavPrev current Is Ir1Bcg Ir2Bcg Ir3Bcg n1 n2 n3
global maxAlpha Ir1B Ir2B Ir3B ptpMeasured ptpNavPrev ptpNav rwalphas
global fsensor MagFieldBias AngFieldBias EulerBias
global MagFieldNoise AngFieldNoise EulerNoise R Amax lmax CD

x = state(1)
y = state(2)
z = state(3)
%xdot = state(4);
%ydot = state(5);
%zdot = state(6);
q0123 = state(7:10);
ptp = Quaternions2EulerAngles(q0123')';
p = state(11);
q = state(12);
r = state(13);
pqr = state(11:13);
w123 = state(14:16);

%%%Translational Kinematics
vel = state(4:6);
%%%Rotational Kinematics

PQRMAT = [  0   -p    -q   -r;  %уравнение для производной кватернионов:
            p    0     r   -q;
            q   -r     0    p;
            r    q    -p    0];

q0123dot = 0.5*PQRMAT*q0123;


%%%Gravity Model
%planet
r = state(1:3); %% r = [x;y;z]
rho = norm(r);
rhat = r/rho;
Fgrav = -(mu*m/rho^2)*rhat;


%%%Call the magnetic field model
if t >= lastMagUpdate       %% Частота обновлений измерений магниитного поля Земли
    lastMagUpdate = lastMagUpdate + nextMagUpdate;
    %%%Преобразование из декартовых координат в географические. (Необходимо
    %%%для вызова модели магнитного поля IGRF)
    phiE = 0;       % Угол поворота Земли (восточная долгота Гринвича) часто принимается 0 при использовании ECI (инерциальной системы) где Земля не вращается
    thetaE = acos(z/rho);   % Полярный угол (колатитуда) угол от оси Z(Оси вращения Земли) до положения спутника
    psiE = atan2(y,x);      % Азимут (Географическая долгота(
    latitude = 90-thetaE*180/pi;
    longitude = psiE*180/pi;
    rhokm = (rho)/1000;
    [BN,BE,BD] = igrf('01-Jan-2020',latitude,longitude,rhokm,'geocentric');
    %%%Convert NED (North East Down to X,Y,Z in ECI frame)
    %%%First we need to create a rotation matrix from the NED frame to the 
    %%%inertial frame
    BNED = [BN;BE;-BD]; %%PCI has Down as Up
    BI = TIB(phiE,thetaE+pi,psiE)*BNED;
    %BI = eye(3)*BNED;    
    BB = TIBquat(q0123)'*BI;
    %%%Convert to Tesla
    BB = BB*1e-9;
end

if t >= lastSensorUpdate
    %%%%SENSOR BLOCK
    lastSensorUpdate = lastSensorUpdate + nextSensorUpdate;
    [BfieldMeasured,pqrMeasured,ptpMeasured] = Sensor(BB,pqr,ptp); 
    
    %%%NAVIGATION BLOCK
    [BfieldNav,pqrNav,ptpNav] = Navigation(BfieldMeasured,pqrMeasured,ptpMeasured);   
end

end

