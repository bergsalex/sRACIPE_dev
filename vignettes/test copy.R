test <- cppFunction() function(){
  i=1
  while(i>0){}
}




cppFunction(plugins=c("cpp11"), '
    int useCpp11() {
            auto x = 10;
            return x;
            }')
