#sed -i '1,2d' out.m
sed -i 's/p1.n1/(&)/g' out.m
sed -i 's/p2.n1/(&)/g' out.m
sed -i 's/p3.n1/(&)/g' out.m
sed -i 's/p4.n1/(&)/g' out.m
sed -i 's/p5.n1/(&)/g' out.m

sed -i 's/p1.n2/(&)/g' out.m
sed -i 's/p2.n2/(&)/g' out.m
sed -i 's/p3.n2/(&)/g' out.m
sed -i 's/p4.n2/(&)/g' out.m
sed -i 's/p5.n2/(&)/g' out.m
sed -i 's/mat =/mat =(/g' out.m
sed -i 's/;/);/g' out.m
