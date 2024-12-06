clear
clc

% Set needed parameters
thickness = 1/39.3701; % inches to meters
ro = 140; % density
k = 0.048;
r=5;
ipts=101;           %:number of points in x direction %do again from 201
tsteps=1400;
c = 628;              %::wave speed (m/s)
dx = thickness/(ipts-1);          %::spacing (m)
diffusivity = k/(c*ro);     %::thermal divisivity
ao = 2000;          %:leftside Dirichlet boundary setting
t = 0; %:start time
initial_temp = 300;

%this is my stuff
counter = 0;
x_LIFE = [];
u_LIFE = [];

% Set grid locations
for i=1:ipts
    x(i) = single(i-1)*dx;
end

% Calc time step
dt = r*dx.^2/diffusivity;  %::time step	(s)

% Initialize BCs
u = initial_temp*ones(ipts,1);
u(1) = ao	;	%:apply Dirichlet bc on left side

% Plot solution and pause
plot(x,u);
xlim([0,ipts])
ylim([0,1.2])
drawnow limitrate%pause(.1)

%Find the A matrix variables
a=-r;b=-r;
d=1+2*r;
%Calculate the A matrix
A = a*diag(ones((ipts-1),1),1) + d*diag(ones((ipts),1),0) + b*diag(ones((ipts-1),1),-1);

%set up the nueman boundary conditions:
A(1,2) = 0; % b0 = 0;
A(1,1) = 1; %d0 = +1;

%set up your nuemann boundary conditions:
A(ipts,ipts-1) = A(ipts-1, ipts-1) + 4*b;
A(ipts, ipts) = A(ipts-1, ipts)-3*b;

% START MAIN LOOP IN TIME "n"
for n=1:tsteps
    t = t + dt;  %::increment time

    unew=A\u; %solve for the new unkown
    unew(1) = ao;	%:apply Dirichlet bc on left side

    if unew(ipts) >= 600
        %linearly interpolate to get exactly 600
        fprintf("Last time is %.5f seconds and the earlier time step is %.5f seconds \n",t,(t-dt))
        fprintf("Last temp is %.5f K and the earlier temp step is %.5f K\n",unew(ipts),u(ipts))
        % Calculate the slope of temperature vs. time
        slope = (unew(ipts) - u(ipts)) / dt;

        % Interpolate to find the exact time when temperature hits 600 K
        time_at_600 = (t - dt) + (600 - u(ipts)) / slope;
        % Display the interpolated time
        fprintf("Interpolated time at 600 K is %.5f seconds\n", time_at_600)
        plotTemperatureProfiles(x_100, u_100, x_200, u_200, x_300, u_300, x_500, u_500, x_700, u_700)
        error('DONE')
    end
    u=unew;
    %Plot solution and pause
    plot(x,u);
    xlim([0,thickness])
    ylim([0,2000])
    %pause(0.5)
    drawnow limitrate
    counter = counter +1;
    if counter == 50
        x_100 = x;
        u_100 = u;
    end
    if counter == 100
        x_200 = x;
        u_200 = u;
    end
    if counter == 150
        x_300 = x;
        u_300 = u;
    end
    if counter == 200
        x_500 = x;
        u_500 = u;
    end
    if counter == 250
        x_700 = x;
        u_700 = u;
    end
end


function plotTemperatureProfiles(x_100, u_100, x_200, u_200, x_300, u_300, x_500, u_500, x_700, u_700)
% Function to plot temperature profiles at specific time steps

figure;
hold on;

% Plot each profile with a label
plot(x_100, u_100, 'DisplayName', 'Time Step 50', 'LineWidth', 1.5);
plot(x_200, u_200, 'DisplayName', 'Time Step 100', 'LineWidth', 1.5);
plot(x_300, u_300, 'DisplayName', 'Time Step 150', 'LineWidth', 1.5);
plot(x_500, u_500, 'DisplayName', 'Time Step 200', 'LineWidth', 1.5);
plot(x_700, u_700, 'DisplayName', 'Time Step 250', 'LineWidth', 1.5);

% Formatting the plot
xlabel('Position (x)');
ylabel('Temperature (u)');
title('Temperature Distribution at Different Time Steps');
legend('show');
xlim([0, max([x_100, x_200, x_300, x_500, x_700])]);
ylim([0, 2000]);
grid on;
hold off;
end
