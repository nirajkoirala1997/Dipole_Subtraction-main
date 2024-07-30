c~~~~~~~[	NIRAJ FILES WRITE HERE	   ]~~~~~~~~~~~~~~~O

       if (norder.eq.0) then

       write(561, 106) xq,model,norder,
     .  (sub_lo(mq,mf,i,io,is),is=1,3),tot_lo

       elseif (norder.eq.1) then
       write(571,106) xq, model, norder, 
     &  (sub_lo2(mq,mf,i,io,is),is=1,3), tot_lo2

       write(581,108) xq, model, norder, deltas, deltac,
     &  (sub_nlo2(mq,mf,i,io,is),is=1,3), tot_nlo2

       write(591, 108) xq, model, norder, deltas, deltac,
     &  (sub_nlo3(mq,mf,i,io,is),is=1,3), tot_nlo3

       write(801,108) xq, model, norder, deltas, deltac,
     &  tot_lo2, tot_nlo2, tot_nlo3, tot_nlo

       write(851,108)xq, model, norder, deltas, deltac,
     &  tot_qqb, tot_qg, tot_gg, tot_nlo

	endif


