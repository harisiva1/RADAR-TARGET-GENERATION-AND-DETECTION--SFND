Radar Target Generation and Detection
Project for Udacity's Sensor Fusion Engineer Nanodegree Program


Project Goal
Configure the FMCW waveform based on the system requirements.
Define the range and velocity of target and simulate its displacement.
For the same simulation loop process the transmit and receive signal to determine the beat signal
Perform Range FFT on the received signal to determine the Range
Towards the end, perform the CFAR processing on the output of 2nd FFT to display the target.
Overview

![OVERVIEW](https://user-images.githubusercontent.com/68550704/121783220-a7354900-cbad-11eb-897d-26a6b5b4e357.png)


 Radar Specifications 

 Frequency of operation = 77GHz
 Max Range = 200m
 Range Resolution = 1 m
 Max Velocity = 100 m/s



rmax = 200;%M
rres = 1;
Maxvel = 100;
c = 3e8;

Operating carrier frequency of Radar 
fc= 77e9;             %carrier freq
Target specifications

User Defined Range and Velocity of target
*%TODO* :
define the target's initial position and velocity. Note : Velocity
remains contant
IR = 110; %initial range
IV = -20; %intial velocity
FMCW Waveform Generation

In this project, we will designing a Radar based on the given system requirements (above).

Max Range and Range Resolution will be considered here for waveform design.

% *%TODO* :
Design the FMCW waveform by giving the specs of each of its parameters.
Calculate the Bandwidth (B), Chirp Time (Tchirp) and Slope (slope) of the FMCW
chirp using the requirements above.

B = c/(2*rres);
Tchirp = (5.5*2*rmax)/c
slope = B/Tchirp;

Then, we need simulate the signal propagation and move target scenario.


![EQUATION](https://user-images.githubusercontent.com/68550704/121783229-b3b9a180-cbad-11eb-9024-fb31c3c56c69.png)

%% Signal generation and Moving Target simulation
% Running the radar scenario over the time. 

for i=1:length(t)         
    
    
    % *%TODO* :
    %For each time stamp update the Range of the Target for constant velocity. 
    at = t(i); %actual time
    
    r_t(i) = IR + IV *at;
    td(i) = (2* r_t(i))/c;
    % *%TODO* :
    %For each time sample we need update the transmitted and
    %received signal.
   
    Tx(i) = cos(2*pi*(fc*at + slope*((at^2)/2)));
    Rx (i) = cos(2*pi*(fc*(at-td(i)) + slope*(((at-td(i))^2)/2)));
    
    % *%TODO* :
    %Now by mixing the Transmit and Receive generate the beat signal
    %This is done by element wise matrix multiplication of Transmit and
    %Receiver Signal
    Mix(i) = Tx(i)* Rx(i);
end
    

Range measurement

The 1st FFT output for the target located at 110 meters
![Range from first FFT](https://user-images.githubusercontent.com/68550704/121782935-625ce280-cbac-11eb-8330-a91dfe6263eb.jpg)
Range and Doppler measurement
2st FFT will generate a Range Doppler Map as seen in the image below and it will be given by variable ‘RDM’.

![Doppler FFT map](https://user-images.githubusercontent.com/68550704/121782941-67219680-cbac-11eb-832c-c5bfab24cd87.jpg)


%%CFAR implementation

The 2D CFAR is similar to 1D CFAR, but is implemented in both dimensions of the range doppler block. The 2D CA-CFAR implementation involves the training cells occupying the cells surrounding the cell under test with a guard grid in between to prevent the impact of a target signal on the noise estimate.

![cfar](https://user-images.githubusercontent.com/68550704/121783244-cfbd4300-cbad-11eb-8255-4f505ff00969.png)

Select the number of Training Cells and Guard Cells in both the dimensions and set offset of threshold

%% CFAR implementation

%Slide Window through the complete Range Doppler Map

% *%TODO* :
%Select the number of Training Cells in both the dimensions.
Tr = 16;
Td = 16;
% *%TODO* :
%Select the number of Guard Cells in both dimensions around the Cell under 
%test (CUT) for accurate estimation
Gr = 4;
Gd = 4;
% *%TODO* :
% offset the threshold by SNR value in dB
off = 3.7;
% *%TODO* :
%Create a vector to store noise_level for each iteration on training cells
noise_level = zeros(1,1);


% *%TODO* :
%design a loop such that it slides the CUT across range doppler map by
%giving margins at the edges for Training and Guard Cells.
%For every iteration sum the signal level within all the training
%cells. To sum convert the value from logarithmic to linear using db2pow
%function. Average the summed values for all of the training
%cells used. After averaging convert it back to logarithimic using pow2db.
%Further add the offset to it to determine the threshold. Next, compare the
%signal under CUT with this threshold. If the CUT level > threshold assign
%it a value of 1, else equate it to 0.


   % Use RDM[x,y] as the matrix from the output of 2D FFT for implementing
   % CFAR
RDM = RDM/max(max(RDM)); 
RDM = db2pow(RDM);
Tsize = (2*Tr+2*Gr+1)*(2*Td+2*Gd+1) - (2*Gr+1)*(2*Gd+1);
for i=(Tr+Gr+1):(Nr/2 -(Tr+Gr))
    for j = (Td+Gd+1):(Nd-(Td+Gd))
        t_r = Tr+Gr+i;
        t_c = Td+Gd+j;
        t_s = sum(sum(RDM((i-(Tr+Gr)):t_r,(j-(Td+Gd)):t_c)));
        s_r = Gr+i;
        s_c = Gd+j;
        subwin = sum(sum(RDM((i-Gr):s_r,(j-Gd):s_c)));
        noise_level = t_s - subwin;
        avg = noise_level/Tsize;
        thresh = pow2db(avg) + off;
        CUT = RDM(i,j);
        if(CUT > thresh)
            RDM(i,j) = 1;
        else
            RDM(i,j) = 0;
        end 
    end
end

RDM(RDM~=0 & RDM~=1) = 0;
The output of the 2D CFAR process，a peak and spread centered at 100m in range direction and -20 m/s in the doppler direction.
![Final](https://user-images.githubusercontent.com/68550704/121782950-7274c200-cbac-11eb-9638-a6972f731c5e.jpg)



