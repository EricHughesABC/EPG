%*********************************************
%************* Bloch simulation **************
%*********************************************

% ############ DEFINE RELEVANT PARAMETERS #################
nmag=512;			% Number of elements of magnetization
ppoints=124;		% Number of points in the pulse MUST match size of RF and gradient file
ppoints=101;        % Philips excitation pulse
ppoints=51;        % Philips 180 pulse
hnmag=nmag/2;

grad=500*1;		% Gradient scaling in Hz/cm
swid=5;				% Width of sample in cm
b1scale=430;		% B1 scaling strength (Hz)
pw=2.054e-3;			% Pulse width (s)
rephase=0;			% FLAG: =1 then does a rephase of slice select following the RF pulse
crush=0;			% FLAG: =1 then crush Mxy prior to pulse to simulate 180 degree pulses and rephase afterwards

tcrush=0.01;        % Crusher gradient duration (seconds)

gcrush=(0.25*nmag-grad*pw)/(2*tcrush);

xxx_axis_mm = linspace(-swid*10/2,+swid*10/2, nmag );
% ############ DEFINE STARTING MAGNETIZATION #################
Mo = 1.0;
Mag = [0, 0, Mo];   % Starting magnetization - change HERE to start from transverse or longitudinal
Ms0=[];

% ############ GET RF AND GRADIENT PULSE SHAPES #################
% Get RF pulse profile to simulate
%  fp=fopen('half_sinc3_124.txt', 'r');
%  pshape=fscanf(fp, '%f', ppoints);
%  fclose(fp);
% 
% fp=fopen('excitation_pulse90.txt', 'r');
% pshape=fscanf(fp, '%f', ppoints);
% fclose(fp);

%pshape = ones(124,1);

% pshape = [    0.1991,
%     0.2000,
%     0.1978,
%     0.1925,
%     0.1840,
%     0.1725,
%     0.1580,
%     0.1405,
%     0.1203,
%     0.0976,
%     0.0725,
%     0.0454,
%     0.0167,
%    -0.0134,
%    -0.0445,
%    -0.0761,
%    -0.1078,
%    -0.1392,
%    -0.1698,
%    -0.1992,
%    -0.2269,
%    -0.2525,
%    -0.2754,
%    -0.2954,
%    -0.3120,
%    -0.3248,
%    -0.3336,
%    -0.3380,
%    -0.3378,
%    -0.3328,
%    -0.3229,
%    -0.3079,
%    -0.2879,
%    -0.2628,
%    -0.2328,
%    -0.1980,
%    -0.1586,
%    -0.1148,
%    -0.0670,
%    -0.0156,
%     0.0392,
%     0.0967,
%     0.1566,
%     0.2183,
%     0.2813,
%     0.3449,
%     0.4086,
%     0.4717,
%     0.5338,
%     0.5941,
%     0.6521,
%     0.7073,
%     0.7590,
%     0.8068,
%     0.8501,
%     0.8886,
%     0.9219,
%     0.9496,
%     0.9715,
%     0.9872,
%     0.9967,
%     0.9999,
%     0.9967,
%     0.9872,
%     0.9715,
%     0.9496,
%     0.9219,
%     0.8886,
%     0.8501,
%     0.8068,
%     0.7590,
%     0.7073,
%     0.6521,
%     0.5941,
%     0.5338,
%     0.4717,
%     0.4086,
%     0.3449,
%     0.2813,
%     0.2183,
%     0.1566,
%     0.0967,
%     0.0392,
%    -0.0156,
%    -0.0670,
%    -0.1148,
%    -0.1586,
%    -0.1980,
%    -0.2328,
%    -0.2628,
%    -0.2879,
%    -0.3079,
%    -0.3229,
%    -0.3328,
%    -0.3378,
%    -0.3380,
%    -0.3336,
%    -0.3248,
%    -0.3120,
%    -0.2954,
%    -0.2754,
%    -0.2525,
%    -0.2269,
%    -0.1992,
%    -0.1698,
%    -0.1392,
%    -0.1078,
%    -0.0761,
%    -0.0445,
%    -0.0134,
%     0.0167,
%     0.0454,
%     0.0725,
%     0.0976,
%     0.1203,
%     0.1405,
%     0.1580,
%     0.1725,
%     0.1840,
%     0.1925,
%     0.1978,
%     0.2000,
%     0.1991,
%     0.1953];

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

maxpshape=max(pshape);
pshape1=pshape/maxpshape;   % Normalise the RF pulse shape - subsequently scaled by parameter "b1scale"

dt=pw/ppoints;              % Time resolution of simulation

% Get accompanying NORMALISED gradient profile to simulate
fp=fopen('grad_124.txt', 'r');
gshape=fscanf(fp, '%f', ppoints);
fclose(fp);

% ############ START OF SIMULATION #################
% Set up spatial magnetization
for j=1:nmag
  Ms0(j,:)=Mag;
end

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
for j=1:nmag
  temp(1,:)=Ms0(j,:);
  
  for k=1:ppoints
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

figure(1);
plot(xxx_axis_mm,my);
hold on
plot(xxx_axis_mm,mx);
grid();
xlabel('mm');
legend({'y','x'});
ylim([-1.1 1.1]);

figure(2);
plot(pshape);
title('pulse shape')

figure(7);
plot(xxx_axis_mm,(mz), '.-');
grid();
xlabel('mm');
ylabel('magnitude');
title('Z magnetization');
ylim([-1.1 1.1]);

% ############ SAVE RESULT #################
save halfsinc3_junk.txt mx my mz
