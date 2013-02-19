/* 21.08.03 changes due to dynamic allocation of nearest, points_dat and lm_estimate  */


#include "Klocnon.h"
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <time.h>



double summe_log(int N)
{
  double n;
  int m;
  
  n=0;
  m=N;
  
  while (m>0){
    n=n+log(m);
    m=m-1;
  }
  return n;
}

#ifdef notdef            
int poisson(double mu)
{
  double c, beta, alfa, k, X, U1, U2;
  int N;
  
  if (mu<30){
     printf("error: invalid value of mu (must be greater than 29)");
     exit(1);
  }
    
  c=0.767-3.36/mu; beta=pi/sqrt(3*mu); alfa=beta*mu; k=log(c)-mu-log(beta);

  do{
   do {
      U1=drand48();
      X=(alfa-log((1-U1)/U1))/beta;
   } while (X<=-0.5);
   N=X+0.5;
   U2=drand48();
  } while (alfa-beta*X+log(U2/pow(1+exp(alfa-beta*X),2))>k+N*log(mu)-summe_log(N));
  
  return N;
}
#endif


     
Pointptr palloc(void);  
     
void init(Pointptr *nearest)
{
  int i, k;
  
  for (i=0; i<max_no_point; i++)
     for (k=0; k<max_no_t+max_ext_t; k++)
                nearest[i*(max_no_t+max_ext_t)+k]=NULL;
}

void insert(int i, int j, Pointptr *nearest, int interval_no, Pointptr point)
{ 
  Pointptr temp;
  
  temp=nearest[i*(max_no_t+max_ext_t)+interval_no];
  nearest[i*(max_no_t+max_ext_t)+interval_no]=palloc();
  nearest[i*(max_no_t+max_ext_t)+interval_no]->next=temp;
  nearest[i*(max_no_t+max_ext_t)+interval_no]->x=point->x;
  nearest[i*(max_no_t+max_ext_t)+interval_no]->y=point->y;
  nearest[i*(max_no_t+max_ext_t)+interval_no]->number=j;
}




Pointptr palloc(void)
{
  return (Point *) malloc(sizeof(Point));
}


void empty_nearest(Pointptr *nearest, int ant_pkt)
{
 int i, k;
 Pointptr p1, p2;
 
 for (i=0; i<ant_pkt; i++)
       for (k=0; k<max_no_t+max_ext_t; k++){
             p1=nearest[i*(max_no_t+max_ext_t)+k];
             while (p1 != NULL){
                 p2=p1->next;
                 free(p1);
                 p1=p2; 
             }
       }
}


void delete(Pointptr p)
{
 Pointptr temp;

 temp=p->next;
 free(p);
 p=temp;
}

void nearest_print(Pointptr *nearest, int ant_pkt)
{
  int i, k;
  Pointptr point;
  
  for (i=0; i<ant_pkt; i++){
     printf("pointnr : %d\n",i); 
     for (k=0; k<max_no_t; k++){
           point=nearest[i*(max_no_t+max_ext_t)+k];
           if (point!=NULL)
              printf("k: %d\n",k); 
           while (point!=NULL){
               printf("%f %f \n",point->x,point->y);
               point=point->next;
           }
     }
  }
}
           
double distance(Point x1, Point x2)
{
   return sqrt(pow(x1.x-x2.x,2)+pow(x1.y-x2.y,2));
}

double min(double a1, double a2)
{
   if (a1<a2)
     return a1;
   else
     return a2;
}

double max(double a1, double a2)
{
   if (a1>a2)
     return a1;
   else
     return a2;
}

 


double transinv_edgecorrection(Point p1, Point p2)
{
  /* frames contains lower left and upper right corner points for 5 frames */
  double frames[]={ 60,  0, 90, 20,
                     0, 40, 30, 60,
                    60, 40, 90, 60,
                   120, 40,150, 60,
                    60, 80, 90,100};
  double area, h1, h2,overlap1, overlap2;
  int i,j;
  
  h1=p2.x-p1.x;
  h2=p2.y-p1.y;

  area=0;

  for (i=0; i<5; i++)
    for (j=0; j<5; j++){
      /* add area of overlap between i'th frame and translated j'th frame */
      overlap1=min(frames[i*4+2],frames[j*4+2]+h1)-max(frames[i*4+0],frames[j*4+0]+h1);
      overlap2=min(frames[i*4+3],frames[j*4+3]+h2)-max(frames[i*4+1],frames[j*4+1]+h2);
      area+=max(overlap1,0.0)*max(overlap2,0.0);
    }
  
  if (area==0){
    printf("area of overlap is 0 h1: %lf h2: %lf p1: %lf %lf p2 %lf %lf\n",h1,h2,p1.x,p1.y,p2.x,p2.y);
    exit(1);
  }
  
  return 1/area;
  
}
double transinv_edgecorrection_whale(Point p1, Point p2){
  double frames[]={
    0,-0.0004551905,0.01462204,0.0003115462,
    0,-0.0003407614,0.03942419,0.000987542,
    0,-0.0008386734,0.03374314,0.0003202328,
    0,-0.000324685,0.02651161,0.0006353672,
    0,-0.0004688558,0.01817779,0.0003383692,
    0,-0.0003140789,0.0007019696,0.0003140789,
    0,-0.0003972983,0.005373567,0.0003161298,
    0,-0.0003186418,0.008111414,0.0003226524,
    0,-0.0003834979,0.008183349,0.0003249686,
    0,-0.0003140789,0.005892508,0.0003496286,
    0,-0.0003140789,0.006082297,0.0003140789,
    0,-0.0003140789,0.00614134,0.000344745,
    0,-0.0003640065,0.00598367,0.0003285977,
    0,-0.0003183022,0.00873438,0.0003608674,
    0,-0.000314079,0.002717756,0.0003140789,
    0,-0.0003914846,0.0131783,0.0003212217,
    0,-0.0003686647,0.01406609,0.0004161915,
    0,-0.0004520547,0.01409101,0.0003271478,
    0,-0.0003774098,0.02184279,0.0004844827,
    0,-0.0005014945,0.01856769,0.0003291226,
    0,-0.0003869667,0.00929209,0.0003233476,
    0,-0.0003383306,0.02659544,0.0006318366,
    0,-0.0003229789,0.01524427,0.0004081764,
    0,-0.0007093381,0.02643268,0.0003386341,
    0,-0.0003163379,0.005684322,0.0003356001,
    0,-0.0003209379,0.0266492,0.0005785641,
    0,-0.000314412,0.01751101,0.0004037908,
    0,-0.0003640529,0.008693457,0.0003210513,
    0,-0.0004215428,0.01506726,0.0003171339,
    0,-0.0003213139,0.01487412,0.0003813101,
    0,-0.0003162161,0.005000514,0.0003439008,
  };
  double area, h1, h2,overlap1, overlap2;
  int l;
  
  h1=p2.x-p1.x;
  h2=p2.y-p1.y;

  area=0;

  for (l=0; l<31; l++){
    /* add area of overlap between l'th frame and translated l'th frame in l'th group */
    overlap1=min(frames[l*4+2],frames[l*4+2]+h1)-max(frames[l*4+0],frames[l*4+0]+h1);
    overlap2=min(frames[l*4+3],frames[l*4+3]+h2)-max(frames[l*4+1],frames[l*4+1]+h2);
    /*overlap2=min(4.002/6371.01,4.002/6371.01+h2)-max(0,h2);*/
    area+=max(overlap1,0.0)*max(overlap2,0.0);
  }


  if (area==0){
    printf("area of overlap is 0\n");
    exit(1);
  }
  
  return 1/area;
  
}


double transinv_edgecorrection_rectangle(Point p1, Point p2, double sideEW, double sideNS)
{
  /* frames contains lower left and upper right corner points for 5 frames */
  double frames[4];
  double area, h1, h2,overlap1, overlap2;
  int i,j;
  
  frames[0]=0; frames[1]=0; frames[2]=sideEW; frames[3]=sideNS;
  h1=p2.x-p1.x;
  h2=p2.y-p1.y;

  area=0;

  /* add area of overlap between frame and translated frame */
  i=0; j=0;
  overlap1=min(frames[i*4+2],frames[j*4+2]+h1)-max(frames[i*4+0],frames[j*4+0]+h1);
  overlap2=min(frames[i*4+3],frames[j*4+3]+h2)-max(frames[i*4+1],frames[j*4+1]+h2);
  area+=max(overlap1,0.0)*max(overlap2,0.0);
  
  if (area==0){
    printf("area of overlap is 0 h1: %lf h2: %lf p1: %lf %lf p2 %lf %lf\n",h1,h2,p1.x,p1.y,p2.x,p2.y);
    exit(1);
  }
  
  return 1/area;
  
}

double transinv_edgecorrection_spruce(Point p1, Point p2)
{
  /* frames contains lower left and upper right corner points for 5 frames */
  double frames[]={ 0,  0, 56, 38};
  double area, h1, h2,overlap1, overlap2;
  int i,j;
  
  h1=p2.x-p1.x;
  h2=p2.y-p1.y;

  area=0;

  /* add area of overlap between frame and translated frame */
  i=0; j=0;
  overlap1=min(frames[i*4+2],frames[j*4+2]+h1)-max(frames[i*4+0],frames[j*4+0]+h1);
  overlap2=min(frames[i*4+3],frames[j*4+3]+h2)-max(frames[i*4+1],frames[j*4+1]+h2);
  area+=max(overlap1,0.0)*max(overlap2,0.0);
  
  if (area==0){
    printf("area of overlap is 0 h1: %lf h2: %lf p1: %lf %lf p2 %lf %lf\n",h1,h2,p1.x,p1.y,p2.x,p2.y);
    exit(1);
  }
  
  return 1/area;
  
}

double transinv_edgecorrection_weed(Point p1, Point p2)
{
  /* frames contains lower left and upper right corner points for each 5 frames in the 9 groups of frames in Brix's weed dataset */
  double frames[]={  
    0,39,30,59,
    60,77,90,97,
    119,39,149,59,
    60,0,90,20,
    60,39,90,59,
    300,38.5,330,58.5,
    359,77.5,389,97.5,
    418,38.5,448,58.5,
    359,0,389,20,
    359,38.5,389,58.5,
    597.5,38.5,627.5,58.5,
    657.5,77.5,687.5,97.5,
    716.5,38.5,746.5,58.5,
    657.5,0,687.5,20,
    657.5,38.5,687.5,58.5,
    0,236,30,256,
    59,275,89,295,
    118.5,236,148.5,256,
    59,197,89,217,
    59,236,89,256,
    299,235.5,329,255.5,
    358,274.5,388,294.5,
    416,235.5,446,255.5,
    358,197,388,217,
    358,235.5,388,255.5,
    598,234,628,254,
    658.5,272,688.5,292,
    718,234,748,254,
    658.5,195,688.5,215,
    658.5,234,688.5,254,
    0,433,30,453,
    58.5,472,88.5,492,
    118,433,148,453,
    58.5,394,88.5,414,
    58.5,433,88.5,453,
    298,433,328,453,
    358,473,388,493,
    416,433,446,453,
    358,394,388,414,
    358,433,388,453,
    597,432,627,452,
    656,471,686,491,
    715,432,745,452,
    656,393,686,413,
    656,432,686,452};
  double area, h1, h2,overlap1, overlap2;
  int i,j,l;
  
  h1=p2.x-p1.x;
  h2=p2.y-p1.y;

  area=0;

  for (l=0; l<9; l++)
    for (i=0; i<5; i++)
       for (j=0; j<5; j++){
       /* add area of overlap between i'th frame and translated j'th frame in l'th group */
          overlap1=min(frames[l*20+i*4+2],frames[l*20+j*4+2]+h1)-max(frames[l*20+i*4+0],frames[l*20+j*4+0]+h1);
          overlap2=min(frames[l*20+i*4+3],frames[l*20+j*4+3]+h2)-max(frames[l*20+i*4+1],frames[l*20+j*4+1]+h2);
          area+=max(overlap1,0.0)*max(overlap2,0.0);
       }


  if (area==0){
    printf("area of overlap is 0\n");
    exit(1);
  }
  
  return 1/area;
  
}

double edge(Point x, Point y, int correctype, double sideEW, double sideNS)
{
  switch(correctype){
  case 0: return transinv_edgecorrection_weed(x,y); break;
  case 1: return transinv_edgecorrection_spruce(x,y); break;
  case 2: return 1.0; break; /* ticks  */
  case 3: return transinv_edgecorrection_rectangle(x,y,sideEW,sideNS); break;
  case 4: return transinv_edgecorrection_whale(x,y); break;
  case 5: return ripley_edgecorrection(x,y,sideEW,sideNS); break;
  default: return transinv_edgecorrection_weed(x,y); break;
  }
}

/* Ripleys correction */

double k1(d1,d2,dist)
double d1; double d2; double dist;
{
double arg1,arg2,kk;

if (dist<=d1)
 arg1=0;
else
 arg1=acos(d1/dist);

if (dist<=d2)
 arg2=0;
else
 arg2=acos(d2/dist);

kk=1/(1-((arg1+arg2)/pi));
return(kk);

}

double k2(d1,d2,dist)
double d1; double d2; double dist;
{
double arg1,arg2,kk;
arg1=acos(d1/dist);
arg2=acos(d2/dist);
kk=1/(0.75-(0.5*(arg1+arg2)/pi));
return(kk);
}

double ripley_edgecorrection(Point p1, Point p2, double sideEW, double sideNS){
  double w, h1, h2, d1, d2, d, t1;

  h1=sideEW-p1.x;
  h2=sideNS-p1.y; 
  d1=min(p1.x,h1);
  d2=min(p1.y,h2);  
  d=sqrt(SQR(d1)+SQR(d2));

  t1=distance(p1,p2);
  if (t1<=d)
    w=k1(d1,d2,t1);
  else
    w=k2(d1,d2,t1);

  return w;
}


double neighbours2(int i, int interval_no, Pointptr *nearest, Point *points, double lm_estimate[], int correctype, double sideEW, double sideNS)
{  
  double K_x1;
  Pointptr px2;

  K_x1=0;
  px2=nearest[i*(max_no_t+max_ext_t)+interval_no];

  while (px2!=NULL){
    K_x1+=edge(points[i],*px2,correctype,sideEW,sideNS)/lm_estimate[px2->number]; 
    px2=px2->next;
  }

  return K_x1/lm_estimate[i];
}
 
double lm(Point p, double *par, int no, int parametric){
  if (parametric)
    return exp(par[0]+par[1]*p.y+par[2]*pow(p.y,2)+par[3]*pow(p.y,3));
  else
    return par[no];
}

void Ktrans(int *ant_pkt_dat, double *pointsx, double *pointsy, double *par, double *max_t, double *K_estimatedat, int *correctype, double *sideEW, double *sideNS, int *parametric)
{
  /* sideEW and sideNS only used when correctype is 3 (rectangle) */
  Pointptr *nearest; 
  int i, j, k, l, interval_no;
  Point *points_dat;
  double discretization_t;
  double *lm_estimate;
  double K201[max_no_t];
  /* variable used for check below  
  FILE *udfil;
  */
  
  lm_estimate=(double *) calloc(*ant_pkt_dat,sizeof(double));
  points_dat=(Point *) calloc(*ant_pkt_dat,sizeof(Point));

  for (l=0; l<*ant_pkt_dat; l++){
     points_dat[l].x=pointsx[l];
     points_dat[l].y=pointsy[l];
     points_dat[l].next=NULL;  
  }
  nearest=(Pointptr *) calloc(*ant_pkt_dat*(max_no_t+max_ext_t),sizeof(Pointptr));
  init(nearest);

  discretization_t=*max_t/max_no_t;

  for (i=0; i<*ant_pkt_dat; i++)
     for (j=i+1; j<*ant_pkt_dat; j++){ 
        interval_no= (int) (distance(points_dat[i],points_dat[j])/discretization_t);
        if ((j!=i) && (interval_no<max_no_t)){ 
           insert(i,j,nearest,interval_no,&(points_dat[j]));
           insert(j,i,nearest,interval_no,&(points_dat[i]));
        } 
     }

  /* nearest_print(nearest,*ant_pkt_dat); */
  /* compute value of intensity surface for each point  */
  for (i=0; i<*ant_pkt_dat; i++)
    lm_estimate[i]=lm(points_dat[i],par,i,*parametric);
  for (k=0; k<max_no_t+1; k++)
         K_estimatedat[k]=0;  /* K_estimate[0] er dummy */

  for (k=1; k<max_no_t+1; k++){
    K_estimatedat[k]+=K_estimatedat[k-1];
    for (i=0; i<*ant_pkt_dat; i++) 
       K_estimatedat[k]+=neighbours2(i,k-1,nearest,points_dat,lm_estimate,*correctype,*sideEW,*sideNS);
  }
  empty_nearest(nearest,*ant_pkt_dat);
  free(lm_estimate);
  free(points_dat);
  free(nearest);
  /* printf("K ok\n");*/

  if (CHECK){
    /* checker med hoved under armen version  */
    for (k=0; k<max_no_t; k++){
      K201[k]=0;
      for (i=0; i<*ant_pkt_dat; i++) 
	for (j=0; j<*ant_pkt_dat; j++)
	  if ((i!=j) && (distance(points_dat[i],points_dat[j])<(k+1)*discretization_t))
	    K201[k]+=edge(points_dat[i],points_dat[j],*correctype,*sideEW,*sideNS)/(lm_estimate[i]*lm_estimate[j]);
    }    
    for (k=0; k<max_no_t; k++)
      if (fabs(K201[k]-K_estimatedat[k+1])>0.000001)
	printf("Error in estimation of K: Khua[%d]=%f K[%d]=%f\n",k,K201[k],k+1,K_estimatedat[k+1]);
    printf("check over\n");
  }

}



void dKtrans_dbeta(int *ant_pkt_dat, double *pointsx, double *pointsy, double *par, double *max_t, double *K_estimatedat, int *correctype, double *sideEW, double *sideNS, int *parametric, double *z, int *p)
{
  /* z contains values of covariates at data points  */
  /* sideEW and sideNS only used when correctype is 3 (rectangle) */
  Pointptr *nearest; 
  int i, j, k, l, interval_no, q;
  Point *points_dat;
  double discretization_t;
  double *lm_estimate;
  double K201[max_no_t];
  /* variable used for check below  
  FILE *udfil;
  */
  
  lm_estimate=(double *) calloc(*ant_pkt_dat,sizeof(double));
  points_dat=(Point *) calloc(*ant_pkt_dat,sizeof(Point));

  for (l=0; l<*ant_pkt_dat; l++){
     points_dat[l].x=pointsx[l];
     points_dat[l].y=pointsy[l];
     points_dat[l].next=NULL;  
  }
  nearest=(Pointptr *) calloc(*ant_pkt_dat*(max_no_t+max_ext_t),sizeof(Pointptr));
  init(nearest);

  discretization_t=*max_t/max_no_t;

  for (i=0; i<*ant_pkt_dat; i++)
     for (j=i+1; j<*ant_pkt_dat; j++){ 
        interval_no= (int) (distance(points_dat[i],points_dat[j])/discretization_t);
        if ((j!=i) && (interval_no<max_no_t)){ 
           insert(i,j,nearest,interval_no,&(points_dat[j]));
           insert(j,i,nearest,interval_no,&(points_dat[i]));
        } 
     }

  /* nearest_print(nearest,*ant_pkt_dat); */
  /* compute value of intensity surface for each point  */
  for (i=0; i<*ant_pkt_dat; i++)
    lm_estimate[i]=lm(points_dat[i],par,i,*parametric);
  for (k=0; k<max_no_t+1; k++)
    for (q=0; q<*p; q++)
         K_estimatedat[k**p+q]=0;  /* K_estimate[0] er dummy */

  for (k=1; k<max_no_t+1; k++){
    for (q=0; q<*p; q++)
      K_estimatedat[k**p+q]+=K_estimatedat[(k-1)**p+q];
    for (i=0; i<*ant_pkt_dat; i++) 
      for (q=0; q<*p; q++)
	K_estimatedat[k**p+q]+=z[i**p+q]*neighbours2(i,k-1,nearest,points_dat,lm_estimate,*correctype,*sideEW,*sideNS);
  }
  for (k=1; k<max_no_t+1; k++)
    for (q=0; q<*p; q++)
      K_estimatedat[k**p+q]*=2;
  
  empty_nearest(nearest,*ant_pkt_dat);
  free(lm_estimate);
  free(points_dat);
  free(nearest);
  /* printf("K ok\n"); */

}




double edge_zero_correction(double t, double bw)
{
  double temp;
  
  temp=t/pow(bw,2);
 
  if (temp>1)
    return 1;
  else
    return 0.5+0.75*temp-0.25*pow(temp,3);
}



void gtrans(int *ant_pkt_dat, double *pointsx, double *pointsy, double *par, double *max_t, double *g_estimatedat, double *bandwidth, int *correctype, double *sideEW, double *sideNS, int *parametric, int *adaptive, int *kerneltype)
{
  Pointptr *nearest;
  Pointptr px2; 
  int i, j, k, l, interval_no;
  Point *points_dat;
  double discretization_t, d;
  double *lm_estimate;
  int lower, upper, pairs1;
  /* variable used for check below  
  double K201;
  int pairs2;
  FILE *udfil; */

  lm_estimate=(double *) calloc(*ant_pkt_dat,sizeof(double));
  points_dat=(Point *) calloc(*ant_pkt_dat,sizeof(Point));

  for (l=0; l<*ant_pkt_dat; l++){
     points_dat[l].x=pointsx[l];
     points_dat[l].y=pointsy[l];
     points_dat[l].next=NULL;  
  }

  nearest=(Pointptr *) calloc(*ant_pkt_dat*(max_no_t+max_ext_t),sizeof(Pointptr));

  init(nearest);

  discretization_t=*max_t/max_no_t;
  
  if (discretization_t*(max_no_t+max_ext_t)<*max_t+*bandwidth){
    printf("max_ext_t too small\n");
    exit(1);
  }
  
  for (i=0; i<*ant_pkt_dat; i++)
     for (j=i+1; j<*ant_pkt_dat; j++){ 
        interval_no= (int) (distance(points_dat[i],points_dat[j])/discretization_t);
        if ((j!=i) && (interval_no<max_no_t+max_ext_t)){ 
           insert(i,j,nearest,interval_no,&(points_dat[j]));
           insert(j,i,nearest,interval_no,&(points_dat[i]));
        } 
     }
  /*  nearest_print(nearest,*ant_pkt_dat); */

  /* compute value of intensity surface for each point  */
  for (i=0; i<*ant_pkt_dat; i++)
    lm_estimate[i]=lm(points_dat[i],par,i,*parametric);

  for (k=0; k<max_no_t; k++)
         g_estimatedat[k]=0;  

  for (k=0; k<max_no_t; k++){
    pairs1=0;
    for (i=0; i<*ant_pkt_dat; i++){ 
       lower=k-(int) (*bandwidth/discretization_t);
       if (lower<0) 
		lower=0;
       upper=k+1+(int) (*bandwidth/discretization_t);
       if (upper>max_no_t+max_ext_t)
                upper=max_no_t+max_ext_t;
       for (j=lower; j<upper; j++){
             px2=nearest[i*(max_no_t+max_ext_t)+j];
             while (px2!=NULL){
                 d=distance(points_dat[i],*px2);
                if (fabs((k+1)*discretization_t-d)<*bandwidth){
		  pairs1++;
		  g_estimatedat[k]+=kernel((k+1)*discretization_t,d,*bandwidth,*adaptive,*kerneltype)*edge(points_dat[i],*px2,*correctype,*sideEW,*sideNS)/(2*pi*((k+1)*discretization_t)*lm_estimate[px2->number]*lm_estimate[i]);
		/**edge_zero_correction((k+1)*discretization_t,*bandwidth));*/                }
                px2=px2->next;
             }
        }
    }
  } 
 empty_nearest(nearest,*ant_pkt_dat);
 free(lm_estimate);
 free(points_dat);
 free(nearest);

#ifdef notdef 
 udfil=fopen("g.out","w");
 for (k=0; k<max_no_t; k++){
  K201=0;
  for (i=0; i<*ant_pkt_dat; i++) 
    for (j=0; j<*ant_pkt_dat; j++)
      if (i!=j){
	d=distance(points_dat[i],points_dat[j]);
        if (fabs((k+1)*discretization_t-d)<*bandwidth){
           pairs2++;
           K201+=epker((k+1)*discretization_t-d,*bandwidth)*edge(points_dat[i],points_dat[j],*correctype)/(2*pi*((k+1)*discretization_t)*lm_estimate[j]*lm_estimate[i]);
        }
      }
      fprintf(udfil,"%lf %lf\n",(k+1)*discretization_t,K201); 
 }
 fclose(udfil);
#endif
}


double kernel(double r, double d, double bw, int adaptive, int type)
{
  switch(type){
  case 1: if (adaptive) return uniform(r-d,min(bw,r)); else return uniform(r-d,bw); break;
  case 2: if (adaptive) return epker(r-d,min(bw,r)); else return epker(r-d,bw); break;
  default: if (adaptive) return uniform(r-d,min(bw,r)); else return uniform(r-d,bw); break;
  }
}

double uniform(double t, double h)
{
  if (fabs(t)<h)
    return 1/(2*h);
  else
    return 0;
}

double epker(double t, double h)
{
    double dd;
    
   
   if (fabs(t)<h)
     dd=(0.75/h)*(1-pow(t/h,2));
      else
     dd=0;   

     
   return(dd);

   }


void Ktrans_cross(int *ant_pkt_dati, double *pointsxi, double *pointsyi, int *ant_pkt_datj, double *pointsxj, double *pointsyj, double *pari, double *parj, double *max_t, double *K_estimatedat, int *correctype, double *sideEW, double *sideNS, int *parametric)
{
  /* sideEW and sideNS only used when correctype is 3 (rectangle) */
  Pointptr *nearest;
  int i, j, k, l, interval_no;
  Point *points_dati;
  Point *points_datj;
  double discretization_t;
  double *lm_estimatei;
  double *lm_estimatej;
  double K201[max_no_t], tempK;
  /* variable used for check below  
  FILE *udfil; */

  lm_estimatei=(double *) calloc(*ant_pkt_dati,sizeof(double));
  points_dati=(Point *) calloc(*ant_pkt_dati,sizeof(Point));
  lm_estimatej=(double *) calloc(*ant_pkt_datj,sizeof(double));
  points_datj=(Point *) calloc(*ant_pkt_datj,sizeof(Point));


  for (l=0; l<*ant_pkt_dati; l++){
     points_dati[l].x=pointsxi[l];
     points_dati[l].y=pointsyi[l];
     points_dati[l].next=NULL;  
  }
  for (l=0; l<*ant_pkt_datj; l++){
     points_datj[l].x=pointsxj[l];
     points_datj[l].y=pointsyj[l];
     points_datj[l].next=NULL;  
  }

  nearest=(Pointptr *) calloc(*ant_pkt_dati*(max_no_t+max_ext_t),sizeof(Pointptr));
  init(nearest);

  discretization_t=*max_t/max_no_t;

  for (i=0; i<*ant_pkt_dati; i++)
     for (j=0; j<*ant_pkt_datj; j++){ 
        interval_no= (int) floor(distance(points_dati[i],points_datj[j])/discretization_t);
	if (interval_no<max_no_t+max_ext_t)
	  insert(i,j,nearest,interval_no,&(points_datj[j]));
     }

  /*nearest_print(nearest,*ant_pkt_dati); */
  /* compute value of intensity surface for each point  */
  for (i=0; i<*ant_pkt_dati; i++)
    lm_estimatei[i]=lm(points_dati[i],pari,i,*parametric);
  for (j=0; j<*ant_pkt_datj; j++)
    lm_estimatej[j]=lm(points_datj[j],parj,j,*parametric);

  for (k=0; k<max_no_t+1; k++)
         K_estimatedat[k]=0;  /* K_estimate[0] er dummy */


  for (k=1; k<max_no_t+1; k++){
    K_estimatedat[k]+=K_estimatedat[k-1];
    for (i=0; i<*ant_pkt_dati; i++) 
       K_estimatedat[k]+=neighbours2_cross(i,k-1,nearest,points_dati,lm_estimatei,lm_estimatej,*correctype,*sideEW,*sideNS);
  }
  empty_nearest(nearest,*ant_pkt_dati);  
  free(lm_estimatei);
  free(points_dati);
  free(lm_estimatej);
  free(points_datj);
  free(nearest);
  printf("Kcross ok\n");

  if (CHECK){
    /* checker med hoved under armen version  */
    for (k=0; k<max_no_t; k++){
      K201[k]=0;
      for (i=0; i<*ant_pkt_dati; i++){
	tempK=0;
	for (j=0; j<*ant_pkt_datj; j++)
	  if (distance(points_dati[i],points_datj[j])<(k+1)*discretization_t){
            if (k==82 && (distance(points_dati[i],points_datj[j])>=(k)*discretization_t))
	      printf("k: %d i: %d j: %d\n",k,i,j); 
	    tempK+=edge(points_dati[i],points_datj[j],*correctype,*sideEW,*sideNS)/lm_estimatej[j];
	  }
	K201[k]+=tempK/lm_estimatei[i];
      }
    }
    printf("hua ok\n");
    for (k=0; k<max_no_t; k++)
      if (fabs(K201[k]-K_estimatedat[k+1])>0.000001)
	printf("Error in estimation of K_cross: Khua[%d]=%f Kcross[%d]=%f\n",k,K201[k],k+1,K_estimatedat[k+1]);
    printf("check ok\n");
  }
}
 

  

double neighbours2_cross(int i, int interval_no, Pointptr *nearest, Point *pointsi, double *lm_estimatei, double *lm_estimatej, int correctype, double sideEW, double sideNS)
{  
  double K_x1;
  Pointptr px2;

  K_x1=0;
  px2=nearest[i*(max_no_t+max_ext_t)+interval_no];

  while (px2!=NULL){    
    /* if (interval_no==82)
       printf("k: %d i: %d j: %d\n",interval_no,i,px2->number); */
    K_x1+=edge(pointsi[i],*px2,correctype,sideEW,sideNS)/lm_estimatej[px2->number]; 
    px2=px2->next;

  }

  return K_x1/lm_estimatei[i];
}
