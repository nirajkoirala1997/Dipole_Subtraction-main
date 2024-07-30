id INT(F1,1,0,0) = 0;

id INT(F1,0,1,0) = 0;

id INT(F1,0,0,1) = 0;

id INT(F1,1,1,0) = 0;

id INT(F1,1,-1,1) = 
  + INT(F1,1,0,1)
    * (-1/2*s);

id INT(F1,0,1,1) = 0;

