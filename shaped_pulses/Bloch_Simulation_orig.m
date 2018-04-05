%*********************************************
%************* Bloch simulation **************
%*********************************************

% ############ DEFINE RELEVANT PARAMETERS #################
nmag=512;			% Number of elements of magnetization
hnmag=nmag/2;

ppoints=124;		% Number of points in the pulse MUST match size of RF and gradient file
grad=8*4250;		% Gradient scaling in Hz/cm
swid=1;				% Width of sample in cm
b1scale=4*375;		% B1 scaling strength (Hz)
pw=0.0005;			% Pulse width (s)
rephase=0;			% FLAG: =1 then does a rephase of slice select following the RF pulse
crush=0;			% FLAG: =1 then crush Mxy prior to pulse to simulate 180 degree pulses and rephase afterwards

tcrush=0.01;        % Crusher gradient duration (seconds)

gcrush=(0.25*nmag-grad*pw)/(2*tcrush);

% ############ DEFINE STARTING MAGNETIZATION #################
Mo = 1.0;
Mag = [0, 0, Mo];   % Starting magnetization - change HERE to start from transverse or longitudinal
Ms0=[];

% ############ GET RF AND GRADIENT PULSE SHAPES #################
% Get RF pulse profile to simulate
fp=fopen('half_sinc3_124.txt', 'r');
pshape=fscanf(fp, '%f', ppoints);
fclose(fp);
 
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
plot(my);
hold on
plot(mx);

% ############ SAVE RESULT #################
save halfsinc3_junk.txt mx my mz
