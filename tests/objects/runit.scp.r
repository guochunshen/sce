
#test the scp function
ro <- getOption("RUnit")
ro$silent <- TRUE
options("RUnit"=ro)

test.scp <- function(){
  win=owin(c(0,10),c(0,10))
  sp=sample(c("A","B","C"),50,replace=TRUE)
  x=runif(50,0,50)
  y=runif(50,0,50)
  scpcom=scp(sp,x,y,win)
  checkEquals(50,scpcom$N)
  checkEquals(3,scpcom$S)
  
  #if two points are in the same point, remove one.
  x[1:2]=c(2,2)
  y[1:2]=c(4,4)
  scpcom=scp(sp,x,y,win,forceUnique=TRUE)
  checkEquals(49,scpcom$N)
  checkEquals(3,scpcom$S)
  
}


test.check_ind_mapped_data<-function(){
  
  win=owin(c(0,10),c(0,10))
  
  #check the null value error in species, x, y and win
  checkException(check_ind_mapped_data(c(1:4),c(1:4),c(1:4),NULL))
  checkException(check_ind_mapped_data(c(1:4),c(1:4),NULL,win))
  checkException(check_ind_mapped_data(c(1:4),NULL,c(1:4),win))
  checkException(check_ind_mapped_data(NULL,c(1:4),c(1:4),win))
  
  #check the warning message in species, x, y
  check_ind_mapped_data(c(1:4,NA),c(1:5),c(1:5),c(1:5),win=win)
  checkTrue(any(names(warnings())=="NA values are contained in species"))


  #check the NA and infinite values in x and y
  checkException(check_ind_mapped_data(c(1:5),c(1:4,NA),c(1:5),win=win))
  checkException(check_ind_mapped_data(c(1:5),c(1:5),c(1:4,NA),win=win))
  checkException(check_ind_mapped_data(c(1:5),c(1:4,Inf),c(1:5),win=win))
  checkException(check_ind_mapped_data(c(1:5),c(1:5),c(1:4,Inf),win=win))
  checkException(check_ind_mapped_data(c(1:5),c(1:4,-Inf),c(1:5),win=win))
  checkException(check_ind_mapped_data(c(1:5),c(1:5),c(1:4,-Inf),win=win))
  
  #check length of the species, x and y
  checkException(check_ind_mapped_data(c(1:4),c(1:5),c(1:5),win=win))
  checkException(check_ind_mapped_data(c(1:5),c(1:4),c(1:5),win=win))
  checkException(check_ind_mapped_data(c(1:5),c(1:5),c(1:4),win=win))
  
  #check type of win
  checkException(check_ind_mapped_data(c(1:5),c(1:5),c(1:5),win=c(1:5)))
  
  #check unique locations of individuals, maximum resolution is 6
  checkException(check_ind_mapped_data(c(1:5,1),c(1:5,1.00000001),c(1:5,1.00000002),win=win,FALSE))
  
  #check the correct format and the return value
  checkTrue(check_ind_mapped_data(c(1:5,1),c(1:5,1.00001),c(1:5,1.00002),win=win))
  dui_index=check_ind_mapped_data(c(1:5,1),c(1:5,1.00000001),c(1:5,1.00000002),win=win,TRUE)
  checkEquals(1, attr(dui_index,"del_index"))
  
}