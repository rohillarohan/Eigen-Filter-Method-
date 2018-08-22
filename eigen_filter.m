%% Eigen Filter Design
clc
close all
clearvars
%% Specification 
N = 50; % Order of the Filter
fs = 24e3; % Sampling Frequency
L = 500; % Integration Grid Density
fVectPi = [0 1.0472 1.3090 1.8326 1.9635 2.2253 2.3562 pi]; % Normaised Frequency Vector
chord = [0 1 4000 1;5000 0 7000 0;7500 0 8500 1;9000 0 12000 0]; % To Plot the Desired Response
%% First band: 0-1.0472
w1 = 0:1.0472/L:1.0472; 
Q1 = zeros((N/2)+1);
for m = 0:N/2
    for n = 0:N/2
        Q1(m+1,n+1) = trapz(w1,cos(m*w1).*cos(n*w1));
    end
end
p1 = zeros(1,(N/2)+1);
D1(0<=w1 & w1<=1.0472)=1; % Desired Response in Band One
d1 = trapz(w1,D1.^2);
for m = 0:N/2
    p1(m+1) = trapz(w1,D1.*cos(m*w1));
end
%% Second band: 1.3090-1.8326
w2 = 1.3090:(1.8326-1.3090)/L:1.8326;
Q2 = zeros((N/2)+1);
for m = 0:N/2
    for n = 0:N/2
        Q2(m+1,n+1) = trapz(w2,cos(m*w2).*cos(n*w2));
    end
end
p2 = zeros(1,(N/2)+1);
D2(1.3090<=w2 & w2<=1.8326)=0; % Desired Response in Band Two
d2 = trapz(w2,D2.^2);
for m = 0:N/2
    p2(m+1) = trapz(w2,D2.*cos(m*w2));
end
%% Third band: 1.9635-2.2253
w3 = 1.9635:(2.2253-1.9635)/L:2.2253;
Q3 = zeros((N/2)+1);
for m = 0:N/2
    for n = 0:N/2
        Q3(m+1,n+1) = trapz(w3,cos(m*w3).*cos(n*w3));
    end
end
p3 = zeros(1,(N/2)+1);
D3(1.9635<=w3 & w3<=2.2253) = (3.8462.*(w3(1.9635<=w3 & w3<=2.2253)))-7.5577; % Desired Response in Band Three
d3 = trapz(w3,D3.^2);
for m = 0:N/2
    p3(m+1) = trapz(w3,D3.*cos(m*w3));
end
%% Fourth band: 2.3562-pi
w4 = 2.3562:(pi-2.3562)/L:pi;
Q4 = zeros((N/2)+1);
for m = 0:N/2
    for n = 0:N/2
        Q4(m+1,n+1) = trapz(w4,cos(m*w4).*cos(n*w4));
    end
end
p4 = zeros(1,(N/2)+1);
D4(2.3562<=w4 & w4<=pi)=0; % Desired Response in Band Four
d4 = trapz(w4,D4.^2);
for m = 0:N/2
    p4(m+1) = trapz(w4,D4.*cos(m*w4));
end
%% Original Filter Design
d = d1+d2+d3+d4;
Q = Q1+Q2+Q3+Q4;
p = p1+p2+p3+p4;
Qt = [Q p';p d];

[V1, D1] = eig(Qt); % Finding the Eigenvectors and Eigenvalues
[Dis1,Ix1]=sort(diag(D1),'ascend'); % Sorting the Eigenvalues to Find the Eigenvector Corresponding to the Minimum Eigenvalue
aHat1 = V1(:,Ix1(1));
sclFact1 = -1/aHat1(end); % Finding the Scale Factor to Get 'a' from 'aHat'
aTemp1 = sclFact1*aHat1;
a1 = aTemp1(1:end-1);
aFlip1 = flip(a1); % Fliping 'a' to Find 'h' 
hTemp1 = [aFlip1(1:end-1)' a1'];
h1 = 1/2*hTemp1;
h1(26)= 2*h1(26); % Doubling the Middle Sample 
[H1,W1] = freqz(h1); % Finding the Frequency Response
HMag1 = abs(H1);

% Plot for Original Filter
figure;
plot((0:12000/512:12000-(12000/512)),1/max(HMag1)*HMag1,'r')
grid on
hold on
x=[chord(:,1) chord(:,3)]';
y=[chord(:,2) chord(:,4)]';
plot(x,y,'b') % Plotting the Desired Response
legend('Actual Response','Desired Response')
axis([0 12000, 0 1])
title('Arbitrary Frequency Response Using Eigenfilter Method')
xlabel('Frequence in Hz')
ylabel('Magnitude Response')
%% Filter Design with Notches
wNotch1 = (4.5/12)*pi; % Frequency 4500 Hz
wNotch2 = (8/12)*pi; % Frequency 8000 Hz
c1 = cos((0:N/2)*wNotch1);
c2 = cos((0:N/2)*wNotch2);
EHat = [c1 0;c2 0]; % Ehat = [E V], V = 0 in this Case
[Qd,Rd] = qr(EHat'); % QR Decomposition to Find Orthonormal Basis of Null Space of Matrix EHat

B = Qd(:,(3:length(Qd))); % Extracting the Null Spaces of Matrix EHat to Form the Orthonormal Basis
%Null Space of EHat: span{V(L+1),V(L+2),.....V((N?2)+1)} L = # of
%Constraints = 2, Therefore Starting from 3-to-(N/2)+1

W0Mat = B'*Qt*B; % Finding the Matrix B'QtB
[V2, D2] = eig(W0Mat); 
[Dis2,Ix2]=sort(diag(D2),'ascend');
w0 = V2(:,Ix2(1)); % Finding Minimum Eigen Vector w0
aHat2 = B*w0; % Finding aHat
% scalFact2 = -1/aHat2(end); % Finding the Scale Factor to Get 'a' from 'aHat'
% aTemp2 = scalFact2*aHat2;
a2 = aHat2(1:end-1);
aFlip2 = flip(a2);
hTemp2 = [aFlip2(1:end-1)' a2'];
h2 = 1/2*hTemp2;
h2(26)= 2*h2(26); % Doubling the Middle Sample 
[H2,W2] = freqz(h2); % Finding the Frequency Response
HMag2 = abs(H2);

% Plot for the Filter with Notches
figure;
plot((0:12000/512:12000-(12000/512)),1/max(HMag2)*HMag2,'r')
grid on
hold on
plot(x,y,'b') % Plotting the Desired Response
legend('Actual Response','Desired Response')
axis([0 12000, 0 1])
title('Arbitrary Frequency Response With Two Notches Using Eigenfilter Method')
xlabel('Frequence in Hz')
ylabel('Magnitude Response')