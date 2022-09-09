function test_suite=test_hldpi2fft()
  initTestSuite;
end
function test_hldpi2fft_()
   % x = rndnorm(0, 1,5,2); 
  x = [-0.0233845632050972   0.9070186193622006;... 
        0.6529594866766634   1.3689145060433903;... 
        0.4477857310723146  -0.6311953712037597;... 
       -1.9256785038579962   0.5886257667993168;... 
       -0.5290011931824666  -0.3602090880229930]; 
  assert(hldpi2fft(x,'gaus',1), [1.004271397497626, 0.639326941754388], 1e-10); 
  assert(hldpi2fft(x,'epan',1), [2.16092716985540, 1.41860783640250], 1e-10); 
  assert(hldpi2fft(x,'biwe',1), [2.56206246064978, 1.68194557166574], 1e-10); 
  assert(hldpi2fft(x,'triw',1), [2.91217497561157, 1.91178781913982], 1e-10); 
  assert(hldpi2fft(x,'tria',1), [2.37840810167410, 1.56137995683763], 1e-10); 
  assert(hldpi2fft(x,'rect',1), [1.71512903132756, 1.12594978591783], 1e-10); 
  assert(hldpi2fft(x,'lapl',1), [0.749151608134926, 0.491803868625919], 1e-10); 
  assert(hldpi2fft(x,'logi',1), [0.554399884886751, 0.363953043939754], 1e-10);
end
