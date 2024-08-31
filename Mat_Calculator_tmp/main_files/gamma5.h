*****************************************************************************
*****************************************************************************
* NO SCHEME
*****************************************************************************

#procedure onlytrace

#do i=1,10;
id,once G(six1?,six1?,?a)=g_(`i',?a) ;
#enddo

#do i=1,10;
tracen `i';
#enddo

#endprocedure

*****************************
*Anticommutation
****************************
#procedure anticommutation

repeat id G(?a,g5,li1?,?c,g5,?d,g5,?e,g5,?b) = -G(?a,li1,g5,?c,g5,?d,g5,?e,g5,?b);
repeat id G(?a,g5,g5,?b)=G(?a,?b);
repeat id G(?a,g5,li1?,?c,g5,?d,g5,?b) = -G(?a,li1,g5,?c,g5,?d,g5,?b);
repeat id G(?a,g5,g5,?b)=G(?a,?b);
repeat id G(?a,g5,li1?,?c,g5,?d) = -G(?a,li1,g5,?c,g5,?d);
repeat id G(?a,g5,g5,?b)=G(?a,?b);
repeat id G(?a,g5,li1?,?c) = -G(?a,li1,g5,?c);

#endprocedure

*****************************************************************************
*LARIN
*****************************************************************************

#procedure larin

#message larin scheme 

#do i=1,6

id,once G(?a,li2?,g5,?b)=(1/6)*eps(li2,gg`i',hh`i',jj`i')*G(?a,gg`i',hh`i',jj`i',?b);
repeat,id eps(?a)= e_(?a);
repeat,contract;

#enddo

***************************************************
*#call gammaprop
*****************************************************

#do i=1,10;
id,once G(six1?,six1?,?a)=g_(`i',?a) ;
#enddo

repeat,id eps(?a)= e_(?a);
repeat,contract;

#do i=1,20;
tracen `i';
#enddo

#endprocedure

******************************************************************************
* T HOOFT
******************************************************************************

#procedure thooft(mat)
#message thooft scheme 


#do i=1,5

id,once G(?a,g5,?b)=(1/24)*eps(gg`i',hh`i',jj`i',kk`i')*G(?a,gg`i',hh`i',jj`i',kk`i',?b);

#enddo

B G;
P +s mat;
 .end

#do i=1,10;
id,once G(six1?,six1?,?a)=g_(`i',?a) ;
#enddo

repeat,id eps(?a)= e_(?a);
repeat,contract;

#do i=1,20;
tracen `i';
#enddo
******************************* GO TO FEYN RULES **********************
*repeat,id eps(?a)= e_(?a);
*repeat,contract;

#endprocedure

******************************************************************************
******************************************************************************
* OWEN
******************************************************************************

#procedure owen

#message owen scheme

repeat id G(?a,g5,li1?,?c,g5,?d,g5,?e,g5,?b) = -G(?a,li1,g5,?c,g5,?d,g5,?e,g5,?b);
repeat id G(?a,g5,g5,?b)=G(?a,?b);
repeat id G(?a,g5,li1?,?c,g5,?d,g5,?b) = -G(?a,li1,g5,?c,g5,?d,g5,?b);
repeat id G(?a,g5,g5,?b)=G(?a,?b);
repeat id G(?a,g5,li1?,?c,g5,?d) = -G(?a,li1,g5,?c,g5,?d);
repeat id G(?a,g5,g5,?b)=G(?a,?b);
repeat id G(?a,g5,li1?,?c) = -G(?a,li1,g5,?c);

*repeat,id G(?a,g5,?b,g5,?c,g5,?d,g5)=G(?a,?b,?c,?d);
*repeat,id G(?a,g5,?b,g5,?c,g5,?d)=G(?a,?b,?c,g5,?d);
*repeat,id G(?a,g5,?b,g5,?c)=G(?a,?b,?c);


#do i=11,20;
id,once G(six1?,six1?,?a,g5,?b)=g_(`i',?a,5_,?b) ;
.sort
#enddo

#do i=1,10;
id,once G(six1?,six1?,?a)=g_(`i',?a) ;
.sort
#enddo

#do i=1,10;
tracen `i';
#enddo

#do i=11,20;
trace4 `i';
#enddo

repeat,contract;

#endprocedure

******************************************************************************
*  MVV 
******************************************************************************

#procedure mvv(mat)

#message mvv scheme 

bracket G;
*print +s mat;
.sort

#do i=1,6
id,once G(?a,li2?,g5,?b,g5,?c)=(1/6)*eps(li2,gg`i',hh`i',jj`i')*G(?a,gg`i',hh`i',jj`i',?b,g5,?c);
*id,once G(?a,g5,?b,g5,?c)=(1/24)*eps(gg`i',hh`i',jj`i',kk`i')*G(?a,gg`i',hh`i',jj`i',kk`i',?b,g5,?c);
#enddo


*Using cyclicity of trace*
id G(six1?,six1?,?a,g5,?b)=G(six1,six1,?b,?a,g5);


repeat,id eps(?a)= e_(?a);
contract;

repeat;
      id,once,G(six1?,six1?,?a,li1?,g5) = distrib_(-2,3,G1,G2,?a)*G3(li1,g5);
      id  G2(?a)*G3(li1?,g5) = eps(?a,li1);
endrepeat;

repeat;
    if ( count(G1,1) );
        id,once,G1(?a) = g_(1,?a);
        Tracen,1;
    endif;
endrepeat;
.sort

*#do i=11,15
*id,once,G1(?a) = g_(`i',?a);
*#enddo
*.sort

*#do i=11,15;
*tracen `i';
*#enddo

id eps(?a) =e_(?a);
repeat, Contract;

#do i=1,10;
id,once G(six1?,six1?,?a)=g_(`i',?a) ;
#enddo

#do i=1,20;
tracen `i';
#enddo

#endprocedure

*****************************************************************************
