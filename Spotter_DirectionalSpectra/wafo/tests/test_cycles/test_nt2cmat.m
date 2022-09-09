function test_suite=test_nt2cmat()
  initTestSuite;
end
function test_nt2cmat_()
    F0 = [0   6   3   8;... 
         0   0   1   6;... 
         0   0   0   3;... 
         0   0   0   0]; 
   NT = cmat2nt(F0); 
   F = nt2cmat(NT); 
    
   assert(NT, [ 0    0    0    0;... 
               17   11    8    0;... 
               24   18   14    0;... 
               27   21   17    0]); 
 
   assert(F, [0   6   3   8;... 
               0   0   1   6;... 
               0   0   0   3;... 
               0   0   0   0]);
end
