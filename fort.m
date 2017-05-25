clear all;
close all;
%*********INPUT PARAMETERS************* 
input = dlmread('PCBIN.txt');
      MSIZE=100;
      A = zeros(MSIZE,MSIZE);
      CGEN = zeros(MSIZE,MSIZE);
      CAP=zeros(MSIZE,MSIZE);
      IDENT=zeros(MSIZE,MSIZE);
      CAP0=zeros(MSIZE,MSIZE);
      IPIV=zeros(MSIZE,1);
      INDXR=zeros(MSIZE,1);
      INDXC=zeros(MSIZE);
      for i1=1:MSIZE
        IDENT(i1,i1)=1.0;
      end 
      V0=2.997925E8 ;
      V02=V0*V0 ;
      EPS=1.0/(V02*4.E-7*pi) ;
      FAC=1./(2.*pi*EPS) ;
      CMTM=2.54E-5; 
      N = input(1);
      NCDIV = input(2);
      IREF = input(3);
      W = input(4);
      SEP = input(5);
      T = input(6);
      ER = input(7);
%*********ERROR CHECKING AND INITIAL COMPUTATIONS*************             
      if (N==1 || N==0)  
        msg = 'INPUT DATA ERROR';
        error(msg)
      end
      if (IREF > N) 
        msg = 'INPUT DATA ERROR';
        error(msg)
      end
      if (N*NCDIV>MSIZE)
        msg = 'INPUT DATA ERROR';
        error(msg)
      end
      DELTA=W/NCDIV;
%**********WITHOUT DIELECTRIC***********************************       
%**********FILL A1**********************************************      
DIST=0;
      for i2 = 1:NCDIV
        A(1,i2)=PHI(DIST);
        A(i2,1)=A(1,i2); 
        DIST=DIST+DELTA;
      end
      
      for i3=2:NCDIV 
          for J=i3:NCDIV
            A(i3,J)=A(1,J-i3+1);
            A(J,i3)=A(i3,J); 
          end
      end
%******FILL A2--AN************************************************
for ii = 2:N
    DIST=(ii-1)*(W+SEP);
    for jj=1:NCDIV
        XJM1=jj-1;
        A(1,jj+(ii-1)*NCDIV)=PHI(DIST+XJM1*DELTA);
        for jj1=2:NCDIV 
            XJM1=jj-1;
            A(jj1,(ii-1)*NCDIV+1)=PHI(DIST-XJM1*DELTA);
            for jj2=2:NCDIV
                XJM1 = jj2-1;
                A(jj2,(ii-1)*NCDIV+1)=PHI(DIST-XJM1*DELTA);
                for K = 2:NCDIV
                    A(K,M+(I-1)*NCDIV)=A((K-1),(M-1)+(I-1)*NCDIV);
                    for M=K:NCDIV
                        A(M,K+(I-1)*NCDIV)=A((M-1),(K-1)+(I-1)*NCDIV);
                    end
                end
            end
        end
    end
end

%******FILL AT2-ATN***************************************************** 
%matrix transpose
A = A';

%*******FILL A MATRIX**************************************************** 
NM1 = N-1;
for i =1:NM1
    NMI=N-i;
    for j=1:NMI
        for K=1:NCDIV
            for L=1:NCDIV 
                A(K+j*NCDIV,L+j*NCDIV+(i-1)*NCDIV)=A(K,L+(i-1)*NCDIV);
            end
        end
    end
    
end
NM2=N-2;
for i=1:NM2
    NMI=N-(I+1);
    for j=1:NMI
        for K=1:NCDIV
            for L=1:NCDIV
                A(K+(j+I)*NCDIV,L+j*NCDIV)=A(K+i*NCDIV,L);
            end
        end
    end
end
%******DETERMINE GENERALIZED CAPACITANCE MATRIX WITHOUT DIELECTRIC****** 
CAP = A;
CAP = inv(CAP);


for i=1:N
    for j =1:N
        sum = 0;
        for K=1:NCDIV
            for M=1:NCDIV
                SUM=SUM+CAP(K+(i-1)*NCDIV,M+(j-1)*NCDIV);
            end
        SUM=SUM+CAP(K+(i-1)*NCDIV,M+(j-1)*NCDIV);
        end
    end
    CGEN(i,j)=SUM*DELTA;     
end
%******DETERMINE TRANSMISSION-LINE CAPACITANCE AND INDUCTANCE*********** 
%******MATRICES WITHOUT DIELECTRIC************************************** 
sum=0;
for i = 1:N
    for j=1:N
        SUM=SUM+CGEN(i,j);
        SUMCAP=SUM;
    end
end
INDEXR=0;
for i = 1:N
    if i == iref
        break
    else
        INDEXR=INDEXR+1 ;
        INDEXC=0;
        for j =1:N
            if J == IREF 
                break;
            else
                INDEXC=INDEXC+1; 
                SUMC=0; 
                SUMR=0;
                for ii=1:N 
                 SUMC=SUMC+CGEN(ii,j);
                 SUMR=SUMR+CGEN(i,ii);
                end
            end
        end
    end
CAP(INDEXR,INDEXC)=CGEN(i,j)-SUMC*SUMR/SUMCAP; 
CAP0(INDEXR,INDEXC)=CAP(INDEXR,INDEXC);    
end
NM= N-1;
CAP = inv(CAP);
INDUCT=CAP/V02;
%      WRITE(6,21) I,J,INDUCT,I,J 
%21    FORMAT(I3,2X,I3,2X,1PE12.5,8X,'=L(',I3,',',I3,')') 
%20    CONTINUE 

%******************WITH DIELECTRIC************************************** 
%******FILL A1********************************************************** 
dist = 0;
for i = 1:NCDIV
    A(1,i)=PHIB(dist);
    A = A';
    dist=dist+DELTA;
end
for i=2:NCDIV
    for j=i:NCDIV
        A(i,j)=A(1,j-i+1);
        A = A';
    end
end
%******FILL A2--AN******************************************************
for i=2:N
    dist=(I-1)*(W+SEP);
    for j=1:NCDIV
        XJM1=(j-1);
        A(1,j+(i-1)*NCDIV)=PHIB(dist+XJM1*DELTA);
    end
    for j = 2:NCDIV
        XJM1=(j-1);
        A(j,(i-1)*NCDIV+1)=PHIB(DIST-XJM1*DELTA);
    end
    for K=2:NCDIV
        for M=2:NCDIV
            A(K,M+(I-1)*NCDIV)=A((K-1),(M-1)+(i-1)*NCDIV); 
            A(M,K+(I-1)*NCDIV)=A((M-1),(K-1)+(i-1)*NCDIV);
        end
    end
end
%******FILL AT2-ATN***************************************************** 
 A = A';
%*******FILL A MATRIX**************************************************** 
NM1=N-1;
for I=1:NM1
    NMI = N - I;
    for J=1:NMI 
        for K=1:NCDIV
            for L=1:NCDIV
                A(K+J*NCDIV,L+J*NCDIV+(I-1)*NCDIV)=A(K,L+(I-1)*NCDIV);
            end
        end
    end
end
NM2 = N-2;
for I = 1:NM2
    NMI=N-(I+1);
    for J=1:NMI 
        for K=1:NCDIV
            for L=1:NCDIV
                A(K+(J+I)*NCDIV,L+J*NCDIV)=A(K+I*NCDIV,L);
            end
        end
        
    end
    
end
      
%******DETERMINE GENERALIZED CAPACITANCE MATRIX WITH DIELECTRIC********* 
CAP = CAP';
CAP =inv(CAP);
for I=1:N 
    for J=1:N 
        SUM=0;
        for K=1:NCDIV
            for M=1:NCDIV
                SUM=SUM+CAP(K+(I-1)*NCDIV,M+(J-1)*NCDIV);
            end
        end
        CGEN(I,J)=SUM*DELTA;
    end
end
      
      
      
      
      
