#procedure SOn
*
*	Procedure to compute color traces for the SO(NF) groups
*	We follow the article by Cvitanovic (Phys.Rev.D14(1976)1536
*
*	We use [T(i),T(j)] = i_*f(i,j,k)*T(k)    (f is the C of Cvitanovic)
*
*	We use the indices i in the space of the fundamental representation
*	The indices j are in the space of the adjoint.
*	The dimension should be the dimension of the fundamental representation.
*
*	We need: (C)Function Tr(cyclic). Indicates the traces.
*	CFunction T, Tp, f(antisymmetric);
*	Symbols a,nf,NF,NA;
*	Indices i1=NF,i2=NF,i3=NF,i4=NF;
*	Indices j1=NA,j2=NA,j3=NA;
*	Dimension NF;
*
*	Usually the value of a is taken to be 1/2;
*	nf is the number of flavors in the fundamental representation.
*	NF is the dimension of the fundamental representation.
*	NA is the dimension of the adjoint representation.
*
*	Routine by J.Vermaseren, 7-jan-1997
*
repeat;
	id,once,Tr(?a) = T(?a,i1,i1)*nf;
	sum i1;
	repeat;
		id,once,T(j1?,j2?,?a,i1?,i2?) = T(j1,i1,i3)*T(j2,?a,i3,i2);
		sum i3;
	endrepeat;
endrepeat;
repeat;
	id,once,f(j1?,j2?,j3?) = 2/a/i_*T(j1,i1,i2)*T(j2,i2,i3)*T(j3,i3,i1);
	sum i1,i2,i3;
endrepeat;
id	T(j1?,i1?,i2?)*T(j1?,i3?,i4?) = Tp(i1,i2,i3,i4);
#do i = 1,1
if ( count(Tp,1) ) redefine i "0";
.sort
id,once,Tp(i1?,i2?,i3?,i4?) =
			a/2*(d_(i1,i4)*d_(i2,i3)-d_(i1,i3)*d_(i2,i4));
#enddo
repeat;
	id	T(j1?,?a,i1?,i2?)*T(j2?,?b,i2?,i3?) =  T(j1,?a,j2,?b,i1,i3);
	id	T(j1?,?a,i1?,i2?)*T(j2?,?b,i3?,i2?) = -T(j1,?a,j2,?b,i1,i3);
	id	T(j1?,?a,i2?,i1?)*T(j2?,?b,i2?,i3?) = -T(j1,?a,j2,?b,i1,i3);
	id	T(j1?,?a,i2?,i1?)*T(j2?,?b,i3?,i2?) =  T(j1,?a,j2,?b,i1,i3);
endrepeat;
id	T(?a,i1?,i1?) = Tr(?a)/nf;
id	Tr(j1?,j2?) = a*d_(j1,j2)*nf;
id	Tr(j1?,j2?,j3?) = i_*a/2*f(j1,j2,j3);
.sort
id	cF = a/2*(NF-1);
id	cA = a*(NF-2);
*id	[cF-cA/6] = a*(NF/3-1/6);
id	cF = a*(NF/3-1/6)+cA/6;
.sort
*
#endprocedure
