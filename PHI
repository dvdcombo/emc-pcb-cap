%C******POTENTIAL WITHOUT THE BOARD ************************************ 
function [DP,DM,PHI]=PHI(D,FAC,DELTA,ER,T)

	 DP=D+DELTA/2
	 DM=D-DELTA/2
	 
	 if(D==0) 
   PHI=-DELTA*FAC*(ALOG(DELTA/2)-1) 
   else
	 PHI=FAC*(DM*ALOG(DM)-DP*ALOG(DP)+DELTA);
	 end
end
