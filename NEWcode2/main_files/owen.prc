******************************************************************************

#procedure owen
#message owen scheme

repeat,id G(?a,g5,?b,g5,?c,g5,?d,g5)=G(?a,?b,?c,?d);
repeat,id G(?a,g5,?b,g5,?c,g5,?d)=G(?a,?b,?c,g5,?d);
repeat,id G(?a,g5,?b,g5,?c)=G(?a,?b,?c);

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
