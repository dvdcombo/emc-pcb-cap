%C******POTENTIAL WITH THE BOARD ************************************ 
function fib = PHIB(D,FAC,DELTA,ER,T,N)

    DP=D+DELTA/2; 
    DM=D-DELTA/2;
    DM2=DM*DM;
    DP2=DP*DP;
    XK=(ER-1)/(ER+1);
    if D==0
         fib=-DELTA*FAC*(log(DELTA/2)-1);
    else
        fib=FAC*(DM*log(DM2)-DP*log(DP2)+2*DELTA)/(ER+1);
        FIN=0; 
        N=0;
    end
    FIN=0;
    N=N+1; 
    TEMP=FIN;
    XN=(N);
    TNT=2*XN*T;
    TNT2=TNT*TNT; 
    F0=(1+XK)*FAC/(ER+1);
    F1=XK*(2*N-1);
    DMOA=DM2/TNT2; 
    DPOA=DP2/TNT2; 
    F2=D*log((DMOA+1)/(DPOA+1))-(DELTA/2)*(log(DM2+TNT2)/log(DP2+TNT2)); 
    F3=2*DELTA+2*TNT*(atan(DM/TNT)-atan(DP/TNT));
    FIN=F0*F1*(F2+F3); 
    FIN=FIN+TEMP; 
    if (abs(FIN-TEMP)>=1E-5)
        %do nothing
    else
       fib=P0+FIN;
       return;   
    end
      
    P0=DELTA*FAC*(2-log((DELTA/2)*2))/(ER+1);
    FIN=0;   
    N=0;
    
    N=N+1;
    TEMP=FIN;
    XN=N;
    TNT=2*XN*T;
    TNT2=TNT*TNT; 
    F1=XK*(2*N-1);
    F2=2-log((DELTA/2)*2+TNT2)-(4/DELTA)*TNT*atan(DELTA/(2*TNT));
    FIN=FAC*DELTA*(1+XK)*F1*F2/(ER+1);
    FIN=FIN+TEMP;
    if abs(FIN-TEMP)>=1e-5
    %do nothing
    else
    fib=P0+FIN;
    end
end