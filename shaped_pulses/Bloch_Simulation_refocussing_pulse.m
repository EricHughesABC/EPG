%*********************************************
%************* Bloch simulation **************
%*********************************************

% ############ DEFINE RELEVANT PARAMETERS #################
nmag=512;			% Number of elements of magnetization
hnmag=nmag/2;

%b1 is 2.0544 uT 
%pulse duration is 2.054 ms
%The accompanying slice selection gradient is 5.7088ms long as its plateau 
%at a strength of 0.1.1422 mT/m. The whole gradient including trapezoidal ends is 7.8088ms.
%The pulse duration is 5.7088 ms, the maximum b1 is 2.9667 uT and the 
%effective time (which is how long the pulse would be if it was a %
%block pulse with the same amplitude and same pulse angle is 1.9792 ms.


ppoints=51;		% Number of points in the pulse MUST match size of RF and gradient file

g1 = 1.1422 % mT/m
b1 = 11.414*1 % uT
p1 = 2.0544 % ms
gamma = 2.67522209e8 %rad/s/T
grad=gamma*g1*1e-3/2/pi/100.		% Gradient scaling in Hz/m  405.3 Hz/cm
swid=5;				% Width of sample in cm
b1scale=(b1*1e-6*gamma/2/pi)		% B1 scaling strength (Hz) 126 Hz

pw=p1*1e-3;			% Pulse width (s)

rephase=0;			% FLAG: =1 then does a rephase of slice select following the RF pulse
crush=0;			% FLAG: =1 then crush Mxy prior to pulse to simulate 180 degree pulses and rephase afterwards

tcrush=0.01;        % Crusher gradient duration (seconds)

gcrush=(0.25*nmag-grad*pw)/(2*tcrush);

xxx_axis_mm = linspace(-swid*10/2,+swid*10/2, nmag );

% ############ DEFINE STARTING MAGNETIZATION #################
Mo = 1.0;
Mag = [0,  0, Mo];   % Starting magnetization - change HERE to start from transverse or longitudinal
Ms0=[];

% ############ GET RF AND GRADIENT PULSE SHAPES #################
% Get RF pulse profile to simulate
% fp=fopen('excitation_pulse90.txt', 'r');
% pshape=fscanf(fp, '%f', ppoints);
% fclose(fp);

% ############ read in Philips pulse ############

pshape = [0.
1214.
2536.
3961.
5477.
7074.
8740.
10461.
12223.
14010.
15805.
17592.
19353.
21071.
22728.
24308.
25794.
27170.
28422.
29536.
30500.
31304.
31938.
32397.
32674.
32767.
32674.
32397.
31938.
31304.
30500.
29536.
28422.
27170.
25794.
24308.
22728.
21071.
19353.
17592.
15805.
14010.
12223.
10461.
8740.
7074.
5477.
3961.
2536.
1214.
0.];

% ############ read in Philips pulse ############
%pshape1_eh = load('-mat', 'W:\EPGcode\Excitation pulse 300118\excitation pulse.mat');

%pshape = pshape1_eh.excitation_pulse;
figure(11); plot(pshape);
 
maxpshape=max(pshape);
pshape1=pshape/maxpshape;   % Normalise the RF pulse shape - subsequently scaled by parameter "b1scale"

%pshape1 = cat(1, pshape1, flipud( pshape1(1:end-1,:)))

dt=pw/ppoints;              % Time resolution of simulation

% Get accompanying NORMALISED gradient profile to simulate
% fp=fopen('grad_124.txt', 'r');
% gshape=fscanf(fp, '%f', ppoints);
% fclose(fp);

gshape = ones(ppoints,1);

%gshape = cat( 1, flipud( gshape), gshape(2:end))

figure(1);
plot(pshape1,'.-');hold on
figure(1); plot(gshape, '.-');
legend({'gradient', 'pulse'});
xlabel('points');
ylabel('magnitude');
title('Overlay of Gradient and Pulse');
ylim([-0.2, 1.1]);


%ppoints = 123+124;

% ############ START OF SIMULATION #################
% Set up spatial magnetization
for j=1:nmag
  Ms0(j,:)=Mag;
end


figure(3); plot(Ms0(:,3));
%exit();

% Crush Mxy if pulse is 180 refocussing 
% Starting magnetization should have been [0, 1, 0] or [1, 0, 0]
if(crush==1)
  for j=1:nmag
    temp1(1)=Ms0(j,1);
    temp1(2)=Ms0(j,2);
    alpha=2.0*pi*tcrush*(gcrush*swid*(j-hnmag)/nmag);
    Ms0(j,1) = (temp1(1)*cos(alpha) + temp1(2)*sin(alpha));
    Ms0(j,2) = (-1*temp1(1)*sin(alpha) + temp1(2)*cos(alpha));
  end
end

% Simulate pulse via Bloch rotations
for j=1:nmag % num isochromats
  temp(1,:)=Ms0(j,:);
  
  for k=1:ppoints % points in shaped pulse
    relb1=b1scale*pshape1(k,:);
    domega=gshape(k,:)*grad*swid*(j-hnmag)/nmag;
 
    theta=atan2(relb1,domega);
    beff=sqrt((relb1^2)+(domega^2));
    phi=2*pi*beff*dt;

    ct=cos(theta);
    st=sin(theta);
    cp=cos(phi);
    sp=sin(phi);

    smag1=temp(k,1);
    smag2=temp(k,2);
    smag3=temp(k,3);
    
    xmag = (smag1*((ct*ct*cp)+(st*st))) + (smag2*(ct*sp)) + (smag3*((st*ct)*(1-cp)));
    ymag = (-1*smag1*(sp*ct)) + (smag2*cp) + (smag3*(st*sp));
    zmag = (smag1*((ct*st)*(1-cp))) - (smag2*(st*sp)) + (smag3*((ct*ct)+(st*st*cp)));

    temp(k+1,1) = xmag;
    temp(k+1,2) = ymag;
    temp(k+1,3) = zmag; 
   end

  Ms0(j,:)=temp(k,:);
end

% Rephase Mxy if pulse is selective EXCITATION
if(rephase==1)
  for j=1:nmag
    temp1(1)=Ms0(j,1);
    temp1(2)=Ms0(j,2);
    alpha=-2.0*pi*(grad*swid*(j-hnmag)/nmag)*(pw/2)
    Ms0(j,1) = (temp1(1)*cos(alpha) + temp1(2)*sin(alpha));
    Ms0(j,2) = (-1*temp1(1)*sin(alpha) + temp1(2)*cos(alpha));
  end
end

% Crush Mxy if pulse was a 180 refocussing
if(crush==1)

% Crusher gradient - rephase
  for j=1:nmag
    temp1(1)=Ms0(j,1);
    temp1(2)=Ms0(j,2);
    alpha=2.0*pi*tcrush*(gcrush*swid*(j-hnmag)/nmag);
    Ms0(j,1) = (temp1(1)*cos(alpha) + temp1(2)*sin(alpha));
    Ms0(j,2) = (-1*temp1(1)*sin(alpha) + temp1(2)*cos(alpha));
  end

% Moving average to produce simulate the effect of the crusher
  width=16;
  while(j<=nmag-width)
    mavx=0;
    mavy=0;
    k=1;
    while(k<=width)
      mavx=mavx+Ms0(j+k,1);
      mavy=mavy+Ms0(j+k,2);
      k=k+1;
    end
    Ms0(j,1) = mavx/width;
    Ms0(j,2) = mavy/width;
    j=j+1;
  end

end
  
mx=Ms0(:,1);
my=Ms0(:,2);
mz=Ms0(:,3);
figure(4);
plot(xxx_axis_mm,mx, '.-');
hold on
plot(xxx_axis_mm,my, '.-');
grid();
xlabel('mm');
ylabel('magnitude');
title('XY magnetization');
legend({'x','y'});
ylim([-1.1 1.1]);

figure(5);
plot(xxx_axis_mm,abs(complex(mx,my)));
xlabel('mm');
figure(6);
plot(xxx_axis_mm,90*(1-mz), '.-');
grid();
xlabel('mm');
ylabel('flip angle');
title('Z magnetization');

figure(7);
plot(xxx_axis_mm,(mz), '.-');
grid();
xlabel('mm');
ylabel('magnitude');
ylim([-1.1 1.1]);
title('Z magnetization');

flip_angle180 = 90*(1-mz);
% ############ SAVE RESULT #################
save halfsinc3_junk.txt mx my mz
save flip_angle.txt flip_angle180