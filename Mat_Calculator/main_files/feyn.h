#procedure feynrules

#do i=1,10
id in(upq(-`i',p1?)) =  U(six`i',p1,mu)*df(cifx`i',p1);
id in(UPQ(-`i',p1?)) = VB(six`i',p1,mu)*df(cifx`i',p1);
id in(glu(-`i',p1?)) = epolglu(lix`i',p1,0)*db(cix`i',p1);
#enddo
*id in(glu(-1,p1)) = epolglu(lix1,p1,0)*db(cix1,p1);
*id in(glu(-3,p2)) = epolglu(lix3,p2,0)*db(cix3,p2);

#do i=1,10
id ou(elt(-`i',p1?)) = UB(six`i',p1,me)*1/2;
id ou(ELT(-`i',p1?)) = V(six`i',p1,me)*1/2;
id ou(upq(-`i',p1?)) = UB(six`i',p1,mu)*df(cifx`i',p1);
id ou(UPQ(-`i',p1?)) = V(six`i',p1,mu)*df(cifx`i',p1);
id ou(glu(-`i',p1?)) = epolglu(lix`i',p1,0)*db(cix`i',p1);
id ou(ph(-`i',p1?)) = epolph(lix`i',p1,0);
id ou(Hig(-`i',p1?)) = 1; 
#enddo

**************************
* Electron-Photon Vertex *
**************************
#do i=1,10
#do j=1,10
#do k=1,10

id Vx(x1?,x2?,ELTeltph,`i',`j',`k',p1?,p2?,p3?) = (i_)*qe*G(si`i',si`j',li`k');

id Vx(x1?,x2?,ELTeltph,`i',`j',-`k',p1?,p2?,p3?) = (i_)*qe*G(si`i',si`j',lix`k');

id Vx(x1?,x2?,ELTeltph,`i',-`j',`k',p1?,p2?,p3?) = (i_)*qe*G(si`i',six`j',li`k');

id Vx(x1?,x2?,ELTeltph,`i',-`j',-`k',p1?,p2?,p3?) = (i_)*qe*G(si`i',six`j',lix`k');

id Vx(x1?,x2?,ELTeltph,-`i',`j',`k',p1?,p2?,p3?) = (i_)*qe*G(six`i',si`j',li`k');

id Vx(x1?,x2?,ELTeltph,-`i',`j',-`k',p1?,p2?,p3?) = (i_)*qe*G(six`i',si`j',lix`k');

id Vx(x1?,x2?,ELTeltph,-`i',-`j',`k',p1?,p2?,p3?) = (i_)*qe*G(six`i',six`j',li`k');

id Vx(x1?,x2?,ELTeltph,-`i',-`j',-`k',p1?,p2?,p3?) = (i_)*qe*G(six`i',six`j',lix`k');

#enddo
#enddo
#enddo


**********************************
** Graviton-SM_Particles Vertex  *
**********************************
#do i=1,10
#do j=1,10
#do k=1,10
*
*----------------------------
* 1. Graviton-Gluon Vertex
*----------------------------

id Vx(x1?, x2?, gluglugr, `i', `j', `k', p1?, p2?, p3?) = 
    (i_)*kg*(-1/2)* (Cgr( li`i',li`j',li`k',ji`k')*p1.p2
              +Dgr( li`i',li`j',li`k',ji`k',p1,p2)
              +Egr( li`i',li`j',li`k',ji`k',p1,p2) )*d_(ci`i',ci`j') ;

id Vx(x1?, x2?, gluglugr,-`i',-`j', `k', p1?, p2?, p3?) = 
    (i_)*kg*(-1/2)* (Cgr( lix`i',lix`j',li`k',ji`k')*p1.p2
              +Dgr( lix`i',lix`j',li`k',ji`k',p1,p2)
              +Egr( lix`i',lix`j',li`k',ji`k',p1,p2) )*d_(cix`i',cix`j');
#enddo
#enddo
#enddo
*
**----------------------------
** 2. Graviton-Photon Vertex
**----------------------------
*
#do i=1,10
#do j=1,10
#do k=1,10

id Vx(x1?, x2?, phphgr, `i', `j', `k', p1?, p2?, p3?) = 
    (i_)*kg*(-1/2)* (Cgr( li`i',li`j',li`k',ji`k')*p1.p2
              +Dgr( li`i',li`j',li`k',ji`k',p1,p2)
              +Egr( li`i',li`j',li`k',ji`k',p1,p2) );

id Vx(x1?, x2?, phphgr,-`i',-`j', `k', p1?, p2?, p3?) = 
    (i_)*kg*(-1/2)* (Cgr( lix`i',lix`j',li`k',ji`k')*p1.p2
              +Dgr( lix`i',lix`j',li`k',ji`k',p1,p2)
              +Egr( lix`i',lix`j',li`k',ji`k',p1,p2) );
#enddo
#enddo
#enddo


***************************
* Electron-Z boson Vertex *
***************************
#do i=1,10
#do j=1,10
#do k=1,10

id Vx(x1?,x2?,ELTeltzbos,`i',`j',`k',p1?,p2?,p3?) = (i_)*gew*(cve*G(si`i',si`j',li`k')-cae*G(si`i',si`j',li`k',g5));

id Vx(x1?,x2?,ELTeltzbos,`i',`j',-`k',p1?,p2?,p3?) = (i_)*gew*(cve*G(si`i',si`j',lix`k')-cae*G(si`i',si`j',lix`k',g5));

id Vx(x1?,x2?,ELTeltzbos,`i',-`j',`k',p1?,p2?,p3?) = (i_)*gew*(cve*G(si`i',six`j',li`k')-cae*G(si`i',six`j',li`k',g5));

id Vx(x1?,x2?,ELTeltzbos,`i',-`j',-`k',p1?,p2?,p3?) = (i_)*gew*(cve*G(si`i',six`j',lix`k')-cae*G(si`i',six`j',lix`k',g5));

id Vx(x1?,x2?,ELTeltzbos,-`i',`j',`k',p1?,p2?,p3?) = (i_)*gew*(cve*G(six`i',si`j',li`k')-cae*G(six`i',si`j',li`k',g5));

id Vx(x1?,x2?,ELTeltzbos,-`i',`j',-`k',p1?,p2?,p3?) = (i_)*gew*(cve*G(six`i',si`j',lix`k')-cae*G(six`i',si`j',lix`k',g5));

id Vx(x1?,x2?,ELTeltzbos,-`i',-`j',`k',p1?,p2?,p3?) = (i_)*gew*(cve*G(six`i',six`j',li`k')-cae*G(six`i',six`j',li`k',g5));

id Vx(x1?,x2?,ELTeltzbos,-`i',-`j',-`k',p1?,p2?,p3?) = (i_)*gew*(cve*G(six`i',six`j',lix`k')-cae*G(six`i',six`j',lix`k',g5));

#enddo
#enddo
#enddo

***********************
* Quark-Photon Vertex *
***********************
#do i=1,10
#do j=1,10
#do k=1,10

id Vx(x1?,x2?,UPQupqph,`i',`j',`k',p1?,p2?,p3?) = (i_)*qu*G(si`i',si`j',li`k')*d_(cif`i',cif`j');

id Vx(x1?,x2?,UPQupqph,`i',`j',-`k',p1?,p2?,p3?) = (i_)*qu*G(si`i',si`j',lix`k')*d_(cif`i',cif`j');

id Vx(x1?,x2?,UPQupqph,`i',-`j',`k',p1?,p2?,p3?) = (i_)*qu*G(si`i',six`j',li`k')*d_(cif`i',cifx`j');

id Vx(x1?,x2?,UPQupqph,`i',-`j',-`k',p1?,p2?,p3?) = (i_)*qu*G(si`i',six`j',lix`k')*d_(cif`i',cifx`j');

id Vx(x1?,x2?,UPQupqph,-`i',`j',`k',p1?,p2?,p3?) = (i_)*qu*G(six`i',si`j',li`k')*d_(cifx`i',cif`j');

id Vx(x1?,x2?,UPQupqph,-`i',`j',-`k',p1?,p2?,p3?) = (i_)*qu*G(six`i',si`j',lix`k')*d_(cifx`i',cif`j');

id Vx(x1?,x2?,UPQupqph,-`i',-`j',`k',p1?,p2?,p3?) = (i_)*qu*G(six`i',six`j',li`k')*d_(cifx`i',cifx`j');

id Vx(x1?,x2?,UPQupqph,-`i',-`j',-`k',p1?,p2?,p3?) = (i_)*qu*G(six`i',six`j',lix`k')*d_(cifx`i',cifx`j');

#enddo
#enddo
#enddo


**********************
* Quark-Gluon Vertex *
**********************
#do i=1,10
#do j=1,10
#do k=1,10

id Vx(x1?,x2?,UPQupqglu, `i', `j', `k',p1?,p2?,p3?) = (i_)*gs*G( si`i', si`j', li`k')*T( cif`i', cif`j', ci`k');
id Vx(x1?,x2?,UPQupqglu, `i', `j',-`k',p1?,p2?,p3?) = (i_)*gs*G( si`i', si`j',lix`k')*T( cif`i', cif`j',cix`k');
id Vx(x1?,x2?,UPQupqglu, `i',-`j', `k',p1?,p2?,p3?) = (i_)*gs*G( si`i',six`j', li`k')*T( cif`i',cifx`j', ci`k');
id Vx(x1?,x2?,UPQupqglu, `i',-`j',-`k',p1?,p2?,p3?) = (i_)*gs*G( si`i',six`j',lix`k')*T( cif`i',cifx`j',cix`k');
id Vx(x1?,x2?,UPQupqglu,-`i', `j', `k',p1?,p2?,p3?) = (i_)*gs*G(six`i', si`j', li`k')*T(cifx`i', cif`j', ci`k');
id Vx(x1?,x2?,UPQupqglu,-`i', `j',-`k',p1?,p2?,p3?) = (i_)*gs*G(six`i', si`j',lix`k')*T(cifx`i', cif`j',cix`k');
id Vx(x1?,x2?,UPQupqglu,-`i',-`j', `k',p1?,p2?,p3?) = (i_)*gs*G(six`i',six`j', li`k')*T(cifx`i',cifx`j', ci`k');
id Vx(x1?,x2?,UPQupqglu,-`i',-`j',-`k',p1?,p2?,p3?) = (i_)*gs*G(six`i',six`j',lix`k')*T(cifx`i',cifx`j',cix`k');

#enddo
#enddo
#enddo

************************
* Quark-Z boson Vertex *
************************
#do i=1,10
#do j=1,10
#do k=1,10

id Vx(x1?,x2?,UPQupqzbos,`i',`j',`k',p1?,p2?,p3?) = (i_)*gew*(cvu*G(si`i',si`j',li`k')-cau*G(si`i',si`j',li`k',g5))*d_(cif`i',cif`j');

id Vx(x1?,x2?,UPQupqzbos,`i',`j',-`k',p1?,p2?,p3?) = (i_)*gew*(cvu*G(si`i',si`j',lix`k')-cau*G(si`i',si`j',lix`k',g5))*d_(cif`i',cif`j');

id Vx(x1?,x2?,UPQupqzbos,`i',-`j',`k',p1?,p2?,p3?) = (i_)*gew*(cvu*G(si`i',six`j',li`k')-cau*G(si`i',six`j',li`k',g5))*d_(cif`i',cifx`j');

id Vx(x1?,x2?,UPQupqzbos,`i',-`j',-`k',p1?,p2?,p3?) = (i_)*gew*(cvu*G(si`i',six`j',lix`k')-cau*G(si`i',six`j',lix`k',g5))*d_(cif`i',cifx`j');

id Vx(x1?,x2?,UPQupqzbos,-`i',`j',`k',p1?,p2?,p3?) = (i_)*gew*(cvu*G(six`i',si`j',li`k')-cau*G(six`i',si`j',li`k',g5))*d_(cifx`i',cif`j');

id Vx(x1?,x2?,UPQupqzbos,-`i',`j',-`k',p1?,p2?,p3?) = (i_)*gew*(cvu*G(six`i',si`j',lix`k')-cau*G(six`i',si`j',lix`k',g5))*d_(cifx`i',cif`j');

id Vx(x1?,x2?,UPQupqzbos,-`i',-`j',`k',p1?,p2?,p3?) = (i_)*gew*(cvu*G(six`i',six`j',li`k')-cau*G(six`i',six`j',li`k',g5))*d_(cifx`i',cifx`j');

id Vx(x1?,x2?,UPQupqzbos,-`i',-`j',-`k',p1?,p2?,p3?) = (i_)*gew*(cvu*G(six`i',six`j',lix`k')-cau*G(six`i',six`j',lix`k',g5))*d_(cifx`i',cifx`j');

#enddo
#enddo
#enddo

***********************
* TRIPLE Gluon Vertex *
***********************
#do i=1,10
#do j=1,10
#do k=1,10

id Vx(x1?,x2?,glugluglu,`i',`j',`k',p1?,p2?,p3?) = gs*f(ci`i',ci`j',ci`k')*(d_(li`i',li`j')*(p1(li`k')-p2(li`k'))+
							              d_(li`j',li`k')*(p2(li`i')-p3(li`i'))+
							              d_(li`k',li`i')*(p3(li`j')-p1(li`j'))); 

id Vx(x1?,x2?,glugluglu,`i',`j',-`k',p1?,p2?,p3?) = gs*f(ci`i',ci`j',cix`k')*(d_(li`i',li`j')*(p1(lix`k')-p2(lix`k'))+
									d_(li`j',lix`k')*(p2(li`i')-p3(li`i'))+
									d_(lix`k',li`i')*(p3(li`j')-p1(li`j'))); 

id Vx(x1?,x2?,glugluglu,`i',-`j',`k',p1?,p2?,p3?) = gs*f(ci`i',cix`j',ci`k')*(d_(li`i',lix`j')*(p1(li`k')-p2(li`k'))+
									d_(lix`j',li`k')*(p2(li`i')-p3(li`i'))+
									d_(li`k',li`i')*(p3(lix`j')-p1(lix`j'))); 

id Vx(x1?,x2?,glugluglu,`i',-`j',-`k',p1?,p2?,p3?) = gs*f(ci`i',cix`j',cix`k')*(d_(li`i',lix`j')*(p1(lix`k')-p2(lix`k'))+
									  d_(lix`j',lix`k')*(p2(li`i')-p3(li`i'))+
									  d_(lix`k',li`i')*(p3(lix`j')-p1(lix`j'))); 

id Vx(x1?,x2?,glugluglu,-`i',`j',`k',p1?,p2?,p3?) = gs*f(cix`i',ci`j',ci`k')*(d_(lix`i',li`j')*(p1(li`k')-p2(li`k'))+
									d_(li`j',li`k')*(p2(lix`i')-p3(lix`i'))+
									d_(li`k',lix`i')*(p3(li`j')-p1(li`j'))); 

id Vx(x1?,x2?,glugluglu,-`i',`j',-`k',p1?,p2?,p3?) = gs*f(cix`i',ci`j',cix`k')*(d_(lix`i',li`j')*(p1(lix`k')-p2(lix`k'))+
									  d_(li`j',lix`k')*(p2(lix`i')-p3(lix`i'))+
									  d_(lix`k',lix`i')*(p3(li`j')-p1(li`j'))); 

id Vx(x1?,x2?,glugluglu,-`i',-`j',`k',p1?,p2?,p3?) = gs*f(cix`i',cix`j',ci`k')*(d_(lix`i',lix`j')*( p1(li`k') - p2(li`k') )
									 +d_(lix`j',li`k')*( p2(lix`i') - p3(lix`i') )
 									  +d_(li`k',lix`i')*( p3(lix`j') - p1(lix`j') )
									  ); 

id Vx(x1?,x2?,glugluglu,-`i',-`j',-`k',p1?,p2?,p3?) = gs*f(cix`i',cix`j',cix`k')*(d_(lix`i',lix`j')*(p1(lix`k')-p2(lix`k'))+
									    d_(lix`j',lix`k')*(p2(lix`i')-p3(lix`i'))+
									    d_(lix`k',lix`i')*(p3(lix`j')-p1(lix`j'))); 

#enddo
#enddo
#enddo

***********************
* FOUR Gluon Vertex *
***********************



**************************
* g g H effective vertex *
**************************
#do i = 1,10
#do j = 1,10

id Vx(x1?,x2?,glugluHig,`i',`j',x3?,p1?,p2?,p3?) = -i_*ch*d_(ci`i',ci`j')*(-d_(li`i',li`j')*p1.p2 + p1(li`i')*p2(li`j'));
id Vx(x1?,x2?,glugluHig,`i',-`j',x3?,p1?,p2?,p3?) = -i_*ch*d_(ci`i',cix`j')*(-d_(li`i',lix`j')*p1.p2 + p1(li`i')*p2(lix`j'));
id Vx(x1?,x2?,glugluHig,-`i',`j',x3?,p1?,p2?,p3?) = -i_*ch*d_(cix`i',ci`j')*(-d_(lix`i',li`j')*p1.p2 + p1(lix`i')*p2(li`j'));
id Vx(x1?,x2?,glugluHig,-`i',-`j',x3?,p1?,p2?,p3?) = -i_*ch*d_(cix`i',cix`j')*(-d_(lix`i',lix`j')*p1.p2 + p1(lix`i')*p2(lix`j'));

*id Vx(x1?,x2?,glugluHig,`i',`j',x3?,p1?,p2?,p3?) = -i_*ch*d_(ci`i',ci`j')*(+d_(li`i',li`j')*p1.p2 + p1(li`i')*p2(li`j'));
*id Vx(x1?,x2?,glugluHig,`i',-`j',x3?,p1?,p2?,p3?) = -i_*ch*d_(ci`i',cix`j')*(+d_(li`i',lix`j')*p1.p2 + p1(li`i')*p2(lix`j'));
*id Vx(x1?,x2?,glugluHig,-`i',`j',x3?,p1?,p2?,p3?) = -i_*ch*d_(cix`i',ci`j')*(+d_(lix`i',li`j')*p1.p2 + p1(lix`i')*p2(li`j'));
*id Vx(x1?,x2?,glugluHig,-`i',-`j',x3?,p1?,p2?,p3?) = -i_*ch*d_(cix`i',cix`j')*(+d_(lix`i',lix`j')*p1.p2 + p1(lix`i')*p2(lix`j'));

#enddo
#enddo

*****************************
** g g g H effective vertex *
*****************************
*#do i=1,10
*#do j=1,10
*#do k=1,10
*
*id Vx(x1?,x2?,gluglugluHig,`i',`j',`k',x3?,p1?,p2?,p3?,p4?) = -gs*ch*f(ci`i',ci`j',ci`k')*( (p3(li`i')-p2(li`i'))*d_(li`j',li`k')
*										  + (p1(li`j')-p3(li`j'))*d_(li`i',li`k')
*                                                                              	  + (p2(li`k')-p1(li`k'))*d_(li`i',li`j'));
*
*id Vx(x1?,x2?,gluglugluHig,`i',`j',-`k',x3?,p1?,p2?,p3?,p4?) = -gs*ch*f(ci`i',ci`j',cix`k')*( (p3(li`i')-p2(li`i'))*d_(li`j',lix`k')
*										  + (p1(li`j')-p3(li`j'))*d_(li`i',lix`k')
*                                                                              	  + (p2(lix`k')-p1(lix`k'))*d_(li`i',li`j'));
*
*id Vx(x1?,x2?,gluglugluHig,`i',-`j',`k',x3?,p1?,p2?,p3?,p4?) = -gs*ch*f(ci`i',cix`j',ci`k')*( (p3(li`i')-p2(li`i'))*d_(lix`j',li`k')
*										  + (p1(lix`j')-p3(lix`j'))*d_(li`i',li`k')
*										  + (p2(li`k')-p1(li`k'))*d_(li`i',lix`j'));
*
*id Vx(x1?,x2?,gluglugluHig,`i',-`j',-`k',x3?,p1?,p2?,p3?,p4?) = -gs*ch*f(ci`i',cix`j',cix`k')*( (p3(li`i')-p2(li`i'))*d_(lix`j',lix`k')
*										  + (p1(lix`j')-p3(lix`j'))*d_(li`i',lix`k')
*                                                                              	  + (p2(lix`k')-p1(lix`k'))*d_(li`i',lix`j'));
*
*id Vx(x1?,x2?,gluglugluHig,-`i',`j',`k',x3?,p1?,p2?,p3?,p4?) = -gs*ch*f(cix`i',ci`j',ci`k')*( (p3(lix`i')-p2(lix`i'))*d_(li`j',li`k')
*										  + (p1(li`j')-p3(li`j'))*d_(lix`i',li`k')
*                                                                              	  + (p2(li`k')-p1(li`k'))*d_(lix`i',li`j'));
*
*id Vx(x1?,x2?,gluglugluHig,-`i',`j',-`k',x3?,p1?,p2?,p3?,p4?) = -gs*ch*f(cix`i',ci`j',cix`k')*( (p3(lix`i')-p2(lix`i'))*d_(li`j',lix`k')
*										  + (p1(li`j')-p3(li`j'))*d_(lix`i',lix`k')
*                                                                             	  + (p2(lix`k')-p1(lix`k'))*d_(lix`i',li`j'));
*
*id Vx(x1?,x2?,gluglugluHig,-`i',-`j',`k',x3?,p1?,p2?,p3?,p4?) = -gs*ch*f(cix`i',cix`j',ci`k')*( (p3(lix`i')-p2(lix`i'))*d_(lix`j',li`k')
*										  + (p1(lix`j')-p3(lix`j'))*d_(lix`i',li`k')
*                                                                              	  + (p2(li`k')-p1(li`k'))*d_(lix`i',lix`j'));
*
*id Vx(x1?,x2?,gluglugluHig,-`i',-`j',-`k',x3?,p1?,p2?,p3?,p4?) = -gs*ch*f(cix`i',cix`j',cix`k')*( (p3(lix`i')-p2(lix`i'))*d_(lix`j',lix`k')
*										  + (p1(lix`j')-p3(lix`j'))*d_(lix`i',lix`k')
*                                                                              	  + (p2(lix`k')-p1(lix`k'))*d_(lix`i',lix`j'));
*
*#enddo
*#enddo
*#enddo
*


*********************
* Photon Propagator *
*********************
#do i=1,10
#do j=1,10

id AA(`i',`j',p1?,x1?)=-i_*d_(li`i',li`j')*phprop(p1);

#enddo
#enddo

**********************
* Z-Boson Propagator *
**********************
#do i=1,10
#do j=1,10

*id zz(`i',`j',p1?,x1?)=d_(li`i',li`j')*zprop(p1);

id zz(`i',`j',p1?,x1?)=i_*(d_(li`i',li`j')-(p1(li`i')*p1(li`j'))/mz^2)*zprop(p1);


#enddo
#enddo

********************
* Gluon Propagator *
********************
#do i=1,10
#do j=1,10

id GG(`i',`j',p1?,x1?)=-i_*d_(li`i',li`j')*gprop(p1)*d_(ci`i',ci`j');

#enddo
#enddo


********************
* Quark Propagator *
********************
#do i=1,10
#do j=1,10

id QQ(`i',`j',p1?,x3?)=i_*(G(si`i',si`j',p1)+x3*G(si`i',si`j'))*fprop(p1)*d_(cif`i',cif`j');

#enddo
#enddo
***********************
* Electron Propagator *
***********************
#do i=1,10
#do j=1,10

id EE(`i',`j',p1?,x3?)=i_*(G(si`i',si`j',p1)+x3*G(si`i',si`j'))*fprop(p1);

#enddo
#enddo

***********************
* Graviton Propagator *
***********************
#do i=1,10
#do j=1,10

*id GR(`i',`j',p1?,x1?)= (i_)*(1/2)*Bgr(li`i',ji`i',li`j',ji`j',p1)*grprop(p1);
id GR(`i',`j',p1?,x1?)= Bgr(li`i',ji`i',li`j',ji`j',p1)*grprop(p1);

#enddo
#enddo
.sort

* LINEARIZE 
id G(?a,aa?)=G(?a,aa);
id Bgr(?a,aa?)=Bgr(?a,aa);
id Dgr(?a,aa?)=Dgr(?a,aa);
id Egr(?a,aa?)=Egr(?a,aa);
#endprocedure
