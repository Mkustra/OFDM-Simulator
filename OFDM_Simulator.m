clear all; close all; clc

N=128; %liczba podnoœnych
M=100;  %liczba wektorów wejœciowych
l_bit=N*M; %liczba bitów wejœciowych
dl_prefiksu=0.2;
SNR=50;
NA=5; %wspó³czynnik nadpróbkowania


%------ 1) Generowanie sygna³u OFDM---------------

% Generowanie ci¹gu wejœciowego
we=rand(1,l_bit);
for k=1:length(we)
    if we(k)>0.5
        we(k)=1;
    else we(k)=-1;
    end
end

% Podzia³ na podnoœne i ifft
we1=reshape(we,N,[]);
we2=ifft(we1);

%------ 2) Dodawanie prefiksu cyklicznego---------------

a=ceil(N.*dl_prefiksu)-1;
a1=size(we2);
pref1=we2(a1(1)-a:a1(1),1:a1(2));
pref2=[pref1;we2];
we_pref=(pref2(:)).'; %sygna³ po dodaniu prefiksu
d=size(pref1);
pref=d(1); 

%------ 3) Nadpróbkowanie-------------------------------

nad=zeros(NA-1,length(we_pref));
we_nad1=[we_pref;nad];
we_nad=we_nad1(:);

%------ 4) Filtracja dolnopasmowa-----------------------

b=size(we_pref);
widmo1=fft(we_nad);
widmo2=[widmo1(1:ceil(0.5*b(2)));zeros(b(2).*(NA-1),1);widmo1((length(we_nad)+1-floor(0.5*b(2)):length(we_nad)))];
we_filtr=NA.*ifft(widmo2);

%------ 5) Dodawanie szumu-------------------------------
% Aby prze³¹czyæ siê na sk³adowe wielodrgowe nale¿y zakomentowaæ liniê 56
% oraz odkomentowaæ kod znajduj¹cy siê pod nag³ówkiem "Sk³adowe
% wielodrogowe"

% wielodrogowych
 %Jedna sk³adowa
     sygnal = awgn(we_filtr,SNR);
 
 %Skladowe wielodrogowe
%  pr=5; %dl przesuniecia
%  sk1=awgn(we_filtr,SNR);
%  sk2prim=[zeros(pr,1);we_filtr(1:length(we_filtr)-pr)];
%  sk2=awgn(sk2prim,SNR);
%  sygnal = sk1+sk2;
 
 
 %----- 6) Odbiór sygna³u OFDM---------------------------

 %Usuniêcie nadpróbkowania
 a=1;
 wy_nad1=zeros(l_bit,1);
for k=1:length(sygnal)
    if  mod(k,NA) == 1 
        wy_nad1(a)=sygnal(k);
        a = a+1;
    end
end
wy_nad2=(wy_nad1).';

%Zamiana wektora na macierz
 wy_nad=reshape(wy_nad2,N+pref,[]);
 
%Usuniêcie prefiksu cyklicznego
 e=size(wy_nad);
 wy_pref=wy_nad((pref+1):e(1),:);
 
 %FFT
 wy_fft=fft(wy_pref);
 
 %----- 7) Korektor----------------------------------------
 wy1=zeros(N,M);
 kprim=zeros(N,10);
 for c=1:10
    kprim(:,c)=we1(:,c)./wy_fft(:,c);
 end
 k=(sum(kprim.').')./10;
 
 for i=1:M
     wy1(:,i)=wy_fft(:,i).*k;
 end
 
 wy=(wy1(:)).';
 
  %----- 8) Obliczanie metryki EVM----------------------------
  EVM=zeros(N,1);
  for v=1:N
  licznik=(1/M).*(sum((real(wy1(v,:))-real(we1(v,:))).^2+(imag(wy1(v,:))-imag(we1(v,:))).^2));
  EVM(v)=sqrt(licznik)*100;
  end
  
  
  %-----------------------------------------------------------------------
  %-----------------------Badania do wykonania----------------------------
  %-----------------------------------------------------------------------
  
  % 1) Widmo przed i po filtracji
  figure(1)
  plot(1:length(widmo1),abs(widmo1),'b',1:length(widmo2),abs(widmo2),'r') 
  title('Widmo sygna³u przed i po filtracji')
  legend('Przed filtracja','Po filtracji')
  xlabel('Czêstotliwoœæ')
  ylabel('Amplituda')
  
  % 2) Fragment sygna³u przed i po filtracji
    %Czesci rzeczywiste
  figure(2)
  plot(1:100,real(we_nad(1:100)),'b',1:100,real(we_filtr(1:100)),'r') 
  title('Sygna³ przed i po filtracji - czêœæ rzeczywista')
  legend('Przed filtracja','Po filtracji')
  xlabel('Numer próbki')
  ylabel('Wartoœæ sygna³u')
    %Czesci urojone
  figure(3)
  plot(1:100,imag(we_nad(1:100)),'b',1:100,imag(we_filtr(1:100)),'r')
  title('Sygna³ przed i po filtracji - czêœæ urojona')
  legend('Przed filtracja','Po filtracji')
  xlabel('Numer próbki')
  ylabel('Wartoœæ sygna³u')
  
  % 3) Histogram
   %Rzeczywiste
   dh=70;
   figure(4)
   hist(real(we_filtr),dh)
   title('Histogram sygna³u - czêœæ rzeczywista')
   xlabel('Wartoœæ')
   ylabel('Iloœæ próbek')
   %Urojone
   figure(5)
   hist(imag(sygnal),dh)
   title('Histogram sygna³u - czêœæ urojona')
   xlabel('Wartoœæ')
   ylabel('Iloœæ próbek')
   %Modul
   figure(6)
   hist(abs(sygnal),dh)
   title('Histogram sygna³u - modu³')
   xlabel('Wartoœæ')
   ylabel('Iloœæ próbek')
   
   %Obliczenie wsp PAPR
   x_peak=max((abs(sygnal)).^2);
   x_rms=std((abs(sygnal)).^2);
   x_rms1=mean((abs(sygnal)).^2);
   PAPR=x_peak./x_rms;
  
  % 4) Zaleznosc EVM od podnosnej
   figure(7)
   plot(1:N,10*log10(EVM))
   title('Wykres parametru EVM od numeru podnoœnej')
   xlabel('Numer podnoœnej')
   ylabel ('EVM')
   xlim ([1 128]) 
   figure(8)
   wid=fftshift(abs(fft(sygnal)));
   plot(30000:47000,wid(30000:47000))
   xlim ([30000 47000])
   title('Widmo sygna³u')
   xlabel('Czêstotliwoœæ')
   ylabel('Amplituda')
   
  % 5) Diagram konstelacji
   %Przed korekcja
   figure(9)
   plot(real(wy_fft(:)),imag(wy_fft(:)),'.')
   axis([-2, 2, -1, 1])
   title('Diagram konstelacji przed korekcj¹')
   xlabel('Czêœæ rzeczywista')
   ylabel('Czêœæ urojona')
   %Po korekcji
   figure(10)
   plot(real(wy),imag(wy),'.')
   axis([-2, 2, -1, 1])
   title('Diagram konstelacji po korekcji')
   xlabel('Czêœæ rzeczywista')
   ylabel('Czêœæ urojona')