%Modified from code provided in class
 clear
  clc
  

thickness = 1/39.3701; % inches to meters
ro = 140; % density
cp = 628;
k = 0.048;
r=0.5;

% Set needed parameters
  ipts=101;           %:number of points in x direction %do again from 201
  tsteps = 15000;     %::number of time steps %PREDICTED TIME is ball park between ranges, second is with fraction
  scheme = 'DNWIND';  %:type of scheme
  c = r;              %::wave speed (m/s)
  dx = thickness/(ipts-1);          %::spacing (m)
  diffusivity = k/(cp*ro);     %::thermal divisivity
  sigma2 = 0.2;	     %::second order dissipation (when activated)
  sigma4 = 0.04 ;    %::fourth order dissipation (when activated)
  ao = 2000;          %:leftside Dirichlet boundary setting
 
% Initialize velocity arrays with zeros
  u(1:ipts) = 300.;   %:current velocity
  unew(1:ipts) = 300.;%:new velocity
  t = 0.; %:start time

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
  u(1) = ao	;	%:apply Dirichlet bc on left side
  u(ipts) = 4/3*u(ipts-1)-1/3*u(ipts-2);  %:apply Neumann bc on right side

% Plot solution and pause
  plot(x,u);
  xlim([0,ipts])
  ylim([0,1.2])
  drawnow limitrate%pause(.1)
  

  
% START MAIN LOOP IN TIME "n"
  for n=1:tsteps
    t = t + dt;  %::increment time
    % START LOOP IN SPACE "i"
    for i=2:ipts-1
        unew(i) = u(i) + r*(u(i+1) - 2*u(i) + u(i-1));  %::calc new u value at each internal point
    end
    
    % Set Boundary Conditions
    unew(1) = ao;	%:apply Dirichlet bc on left side
    unew(ipts) = 4/3*unew(ipts-1)-1/3*unew(ipts-2);  %:apply Neumann bc on right side
    if unew(ipts) >= 600
        disp(t)
        disp(t-dt) %linearly interpolate to get exactoly 600
        fprintf("Last time is %.5f seconds and the earlier time step is %.5f seconds \n",t,(t-dt))
        fprintf("Last temp is %.5f K and the earlier temp step is %.5f K\n",unew(ipts),u(ipts))
    % Calculate the slope
    slope = (t - (t - dt)) / (unew(ipts) - unew(ipts - 1));
    
    % Interpolated time when temperature hits 600
    time_at_600 = (t - dt) + slope * (600 - unew(ipts - 1));
    
    % Display the interpolated time
    fprintf("Interpolated time at 600 K is %.5f seconds\n", time_at_600)

        error('DONE')
    end
  
    %Update solution to current time level
    u(1:ipts) = unew(1:ipts);  %::update solution
    
    %Plot solution and pause
    plot(x,u);
    xlim([0,thickness])
    ylim([0,2000])
    %pause(0.5)  
    drawnow limitrate
    counter = counter +1;
    if counter >=2000
      counter = 0;
      %x_LIFE = [x_LIFE,x];
      %u_LIFE = [u_LIFE,u];
    end 




  end
  %plot(x_LIFE,u_LIFE)
