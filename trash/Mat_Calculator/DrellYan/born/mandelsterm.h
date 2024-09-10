#procedure mandelsterm

id p1.p2 = s/2;
id p3.p4 = s/2;
id p1.p3 = -t/2;
id p2.p4 = -t/2;
id p1.p4 = -u/2;
id p2.p3 = -u/2;
id p1.p5 = s15/2;

#endprocedure

#procedure mandelsterm2


#do ii = 1,5
#do jj = 1,5

#if `ii'>`jj'
.sort
L  s`ii'`jj'= s`jj'`ii';
.sort
#endif

#enddo
#enddo


id p1.p1=0;
id p2.p2=0;
id p3.p3=0;
id p4.p4=0;
id p5.p5=0;

#do ii = 1,5
#do jj = 1,5

id p`ii'.p`jj'= s`ii'`jj'/2;

#enddo
#enddo

#do ii = 1,5
#do jj = 1,5

id fprop(p`ii'-p`jj')= -1/s`ii'`jj';
id fprop(p`ii'+p`jj')= 1/s`ii'`jj';

#do kk = 1,5
id fprop(p`ii'+p`jj'+p`kk')= 1/(s`ii'`jj' + s`jj'`kk' + s`kk'`ii');
id fprop(p`ii'+p`jj'-p`kk')= 1/(s`ii'`jj' - s`jj'`kk' - s`kk'`ii');
id fprop(p`ii'-p`jj'+p`kk')= 1/(-s`ii'`jj' - s`jj'`kk' + s`kk'`ii');
id fprop(p`ii'-p`jj'-p`kk')= 1/(-s`ii'`jj' + s`jj'`kk' - s`kk'`ii');

#enddo
#enddo
#enddo

#endprocedure
