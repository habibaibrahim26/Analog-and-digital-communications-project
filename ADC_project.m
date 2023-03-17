clear;
clc;
close all;
S = input('Enter the amplitudes of the M signals in a matrix with M rows :');
T = input('Enter the duration of Smallest pulse in the signals:');

%Test 1
%{
[1 0 0;1 1 0;0 1 1;1 1 1]
T=1/3;
%}

%Test 2
%{
[2 0 0;-4 -4 0;3 3 3]
T=1;
%}
plot_signal(S,T);
plot_basis(S,T);
constellation(S,T);
%%
%entering S, T for polar NRZ and get there curves
%{
S=[1;-1]
T=1
%}
No_bits=10^5;   %Number of bits stream 
In_stream_bits=randi([0 1],No_bits,1); %random stream generation of zeroes and ones
%Stream of bits encoding into binary Polar NRZ %
for i=1:No_bits
    if In_stream_bits(i)==1
        In_stream_bits(i)=1;
    else
       In_stream_bits(i)=-1;
    end
end
P_NRZ_stream=In_stream_bits;
%% %Generating AWGN with different Levels & calculating noise power spectral density
Eb_No_dB = [-10 -8 -6 -4 -2 0 2 4 6];  %Eb/No
Eb_No = 10.^(Eb_No_dB/10);
SNR = 2*Eb_No;
SNR_dB = pow2db(SNR); %SNR = 2*Eb/No
 %Transmit the encoded stream through AWGN channel
for i=1:length(Eb_No_dB)
  Received_bits(i,:) = awgn(P_NRZ_stream(:) ,SNR_dB(i),'measured');
end

%% %sketching the constellation of the received signal

S_ij=zeros(length(Eb_No_dB),length(P_NRZ_stream));
for i=1:length(Eb_No_dB)
    for j=1:No_bits
    S_ij(i,j)=Received_bits(i,j);
    end
end

y = zeros(1,No_bits);
for i=1:length(Eb_No_dB)
    figure(5)
    subplot(length(Eb_No_dB),1,i);
    scatter(S_ij(i,:),y,'c x','LineWidth',0.5);
    xlim([min(S_ij(1,:))-1   max(S_ij(1,:))+1]);
    xlabel('ϕ_1');
    grid on;
    title(sprintf('Eb/No = %d dB',Eb_No_dB(i)));
 end

%% %Decoding received signal if received bit>0 it will be 1 , if <0 it is -1
P_NRZ_Decoded_Signal=zeros(length(Eb_No_dB),No_bits);
for i=1:length(Eb_No_dB)
   for j=1:No_bits
     if Received_bits(i,j)>0
       P_NRZ_Decoded_Signal(i,j)=1;
          else
       P_NRZ_Decoded_Signal(i,j)=-1;
     end 
  end 
end
%Comparing the Decoded bits by the Original one to get The BER%
for i =1:length(Eb_No_dB) %
 No_error_bit(i) = 0;
  for j=1:No_bits
    if P_NRZ_Decoded_Signal(i,j) ~= P_NRZ_stream(j) % checking for error bits
      No_error_bit(i) = No_error_bit(i)+1;
    end
  end
BER_measured(i) = No_error_bit(i)/No_bits;
end
BER_measured    %Measured bit error rate 
BER_theortical = 1/2.*erfc(sqrt(Eb_No))  %Theortical bit error rate 
%%
%Plotting Measured & Theortical BER vs Eb/No on the same curve
figure
semilogy(Eb_No_dB,BER_measured,'k x -')
grid
hold
     ylabel('Bit error rate');
     xlabel('Eb/No (dB)');
     title('Measured &Theortical BER vs Eb/No')
semilogy(Eb_No_dB,BER_theortical,'g + --')
legend('Measured BER','Theortical BER')
%%
%{
. Comment on the result: 
From the numbers and the curves we can see that measured and theoretical BER are very nearly the same which is because we have big number of symbols
Also, by increasing SNR the BER decreases which is predicted as the signal power increases, the noise effect decreases and also error decreases.

%}
 %%
function plot_signal(S,T)% plot the signals as pulses with pulse duration T
M=size(S,1); %M is the number of signals
for i=1:1:M
    figure 
plot_pulses(S(i,:),T) %plot all signals in case of pulse signals
ylabel(['S ',num2str(i)]);
xlabel('t');
end
end

function plot_basis(S,T) % plot the basis functions as pulses with pulse duration T
phi=basis_func(S,T);
N_basis(phi);
M=size(S,1);
for i=1:1:M
       %if (isnan(phi(i,:))~=true(1,length(phi(i,:)))) 
               figure
               plot_pulses(phi(i,:),T)
               ylabel(['phi ',num2str(i)]);
               xlabel('t');
       %end    
end
end
%%
function plot_pulses(S,T) % it takes a (vector) of pulse amplitudes, turns it into a stream of pulses and plots them
    
    L=length(S); %L counts the number of the signals
    dt=T/1000;  
    fs=T/dt;  
    stream = zeros(1,L*fs); % stream of Zeros along the figure
    
    for i=1:1:L
        if i==1 
            stream(1)=S(i); 
        end
        stream((i-1)*fs+2:(i*fs))=S(i).*ones(1,fs-1);  
        if i==L
            stream(i*fs+1)=0;            
        else
            stream(i*fs+1)=(S(i)+S(i+1))/2;
        end 
    end
   
    t=[0:dt:T*L];
    plot(t,stream) %plot all the signal as a stream of continuous pulses 

end

function E=E_calc(V,T)
    E=sum((V.^2)*T);    %Calculate Energy of the pulse signal
end

function basis_func=basis_func(S,T)% it takes a matrix of the signal (pulse) amplitudes, their pulse duration to generate a matrix of phi amplitudes with the same pulse duration 

L=length(S(1,:));
M= size (S,1);
N=0;
phi=zeros(M,L);   % initialize phi
g=zeros(1,L);   % initialize g
phi(1,:)=S(1,:)/sqrt(E_calc(S(1,:),T));

for i=2:1:M
 
    g=S(i,:);
    for j=1:1:i-1
        
        Sij=sum(S(i,:).*phi(j,:)*T); % Sij=integration(S_i*phi_j)
        g=g-(Sij.*phi(j,:));
  
    end
    
        phi(i,:)=g/sqrt(E_calc(g,T));
end
basis_func=phi;
end

function phies=phies(phi) % it produces a matrix of phies (row1=phi1,row2=phi2) without any NAN (0/0) values by giving it the phi matrix that has NAN values
   M=size(phi,1);
   j=0;
   for i=1:1:M  
       if (isnan(phi(i,:))~=logical(ones(1,length(phi(i,:))))) % don't count all NAN (0/0) rows (unused phies)
              j=j+1;  
              phies(j,:)=phi(i,:);                     
       end   
    end 
end

function N_basis=N_basis(phi) %calculates the number of basis functions N that will be used in the constillation diagram 
    N_basis=size(phies(phi),1);
end

function Sij=Sij(S,phi,T) % generates a matrix of S in terms of phies the 1st row --> S11, the 2nd row --> S12 ,S22  etc

L=length(S(1,:));
N=N_basis(phi);   % the number of the matrix colomns = the number of basis funcs 
M= size (S,1);    % the number of the matrix rows = the number of the signals 
Sij=zeros(M,N);   % initialise the matrix with zeros 

for i=1:1:M
 
    for j=1:1:i
        
        if (isnan(phi(j,:))==logical(ones(1,length(phi(i,:))))) 
            phi(j,:)=zeros(1,L);     %all NAN (0/0) values will be replaced with zeros  
        end  
        
        Sij(i,j)=sum(S(i,:).*phi(j,:)*T);  % Sij=integration(S_i*phi_j)
    end   
end

end

function constellation(S,T)
phi=basis_func(S,T);%Basic functions as a matrix
s = Sij(S,phi,T);    
columnNumbers = find(sum(abs(s)) == 0); %return the coulumn number which has a zero in all rows of it
if columnNumbers >= 1 %if any column in the matrix of phi = 0
L = length(s)-1; %Number of colums will decrease by one
else L = length(s); % the length will be the same if no column in the matrix of phi = 0
end
if L == 1 % if we have only one basic function
     scatter(s(1:end,1),s(1:end,columnNumbers));
     xlabel('phi1(t)');
     title('Constellation Diagram')
else if L == 2  % if we have Two basic functions
     scatter(s(1:end,1),s(1:end,2));
     ylabel('phi2(t)');
     xlabel('phi1(t)');
     title('Constellation Diagram')
else if L == 3   % if we have Three basic functions
     scatter3(s(1:end,1),s(1:end,2),s(1:end,3));
     zlabel('phi3(t)');
     ylabel('phi2(t)');
     xlabel('phi1(t)');
     title('Constellation Diagram')
    else
     scatter(s(1:end,1),zeros(length(s),1));   
     title('Constellation Diagram')
end 

M = length(s); % Number of M signals
E=zeros(1,M);
for i = 1:1:M
  E(i) = sum(s(i,1:end).^2); %Calculating Energy of each symbol from the Matrix of the constillation Diagram
end
Energy=E
end

end
end