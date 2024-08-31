#procedure kpositive
*
id Prop(-k1,0) = Prop(k1,0);
id 1/Prop(-k1,0) = 1/Prop(k1,0);
*
id Prop(-k1+p1,0) = Prop(k1 - p1,0);
id 1/Prop(-k1+p1,0) = 1/Prop(k1 - p1,0);
*
id Prop( - k1 + p1 + p2,0) = Prop(  k1 - p1 - p2,0);
id 1/Prop( - k1 + p1 + p2,0) = 1/Prop(  k1 - p1 - p2,0);
*
*
id Prop(k1,0) = propA;
id 1/Prop(k1,0) = 1/propA;
*
id Prop(k1 - p1,0) = propB;
id 1/Prop(k1 - p1,0) = 1/propB;
*
id Prop(  k1 - p1 - p2,0) = propC;
id 1/Prop(  k1 - p1 - p2,0) = 1/propC;
*
#endprocedure
