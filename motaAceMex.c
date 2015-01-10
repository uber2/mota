#include "mex.h"

/* 7.8.00 By Henning Voss and Markus Abel */

#include<stdio.h>
#include<math.h>
#include<stdlib.h>
#include<time.h>

void aceMEX(double x[],double Parameter[],double Table_Phi[],double *phi,int *ind,int *iind,double *sum,double *yterm,double* cey,double* yind)
{
/* ----------------------------------------------------- */  
void quicksort(double *vet, int *indvet, int s1, int s2);
int window(double x[], int length, int wl, double out[]);
int cef(double y[], int ind[], int iind[], int length, int wl, 
	 double out[],double* cey,double* yind);
double corr(double *x, double *y, int length);
/* ----------------------------------------------------- */ 

   double mean, variance, stdev, esqo, esqi, oscreal, iscreal;
  /* variables in main program */
  int d, dd, itero, iteri; 
  int i, j;
  /*FILE *fp;*/
 
  
  /* integer parameters */
  int length;
  int dim, wl, oi, ii; 
  /* flags */
 
  int showoloop, showiloop;
  /* double parameters */
  double osc, isc;
    /* data field */
    /*double *x;*/

  
     j=0;
     length     =   (int)Parameter[0];
     dim        =   (int)Parameter[1];
     wl         =   (int)Parameter[2];
     oi         =   (int)Parameter[3];
     ii         =   (int)Parameter[4];
     showoloop  =   (int)Parameter[5];
     showiloop  =   (int)Parameter[6];
     osc        =   (int)Parameter[7];
     isc        =   (int) Parameter[8];




  
  /* calculation of sorted indices ind[] */
  for(d=0;d<=dim;d++) for(i=0;i<length;i++) ind[d*length+i]=i;
  for(d=0;d<=dim;d++) quicksort(&x[d*length],&ind[d*length],0,length-1);

  /* transformation to rank ordered variables */
  for(d=0;d<=dim;d++) 
    for(i=0;i<length;i++) iind[d*length+ind[d*length+i]]=i;

  /* calculation of normalized (mean 0, stdev 1) initial function: */
  /* mean and stdev of field (0,1,...,length-1) by formula */
  mean=(length-1.)/2.; 
  /*  stdev=sqrt((length*length-1.)/12.); */
  stdev=sqrt((double)length)*sqrt((double)length-1.)/sqrt(12.);
  for(d=0;d<=dim;d++) 
    for(i=0;i<length;i++) 
      phi[d*length+i]=((double)iind[d*length+i]-mean)/stdev;

  /* outer loop */
  itero=0;
  oscreal=1; /* to calculate outer stopping crit. */
  while(itero<=oi && osc<oscreal){

    /* inner loop */
    iteri=0; iscreal=1;
    while(iteri<=ii && isc<iscreal){
      
      /* loop over rhs */
      for(d=1;d<=dim;d++){
	for(i=0;i<length;i++) sum[i]=0;
	for(dd=1;dd<=dim;dd++){
	  /* summing up all but the phi(d) on the rhs */
          if(dd != d) for(i=0;i<length;i++) sum[i] += phi[dd*length+i];
	}
	for(i=0;i<length;i++) yterm[i] = phi[i]-sum[i];
	cef(yterm,&ind[d*length],&iind[d*length],length,wl,&phi[d*length],cey,yind);
      } /*end of loop over rhs */
      
    /* calculation of inner iteration criterion */
    iscreal=esqi; /* used as memory of old esqi */
    esqi=0;
    for(i=0;i<length;i++) sum[i]=0;
    for(d=1;d<=dim;d++) for(i=0;i<length;i++) sum[i] += phi[d*length+i]; 
    /* (sum is also rhs for outer loop) */
    for(i=0;i<length;i++) 
      esqi+=(sum[i]-phi[i])*(sum[i]-phi[i])/length;
    iscreal=fabs(iscreal-esqi);
    if (showiloop) 
      printf("iteri=%5d yields e^2: %1.6lf (difference: %e)\n",
	     iteri,esqi,iscreal);

    iteri ++ ;
    }/*end of inner loop*/

    cef(sum,&ind[0],&iind[0],length,wl,&phi[0],cey,yind);
    /*normalization: */
    mean=0; variance=0; 
    for(i=0;i<length;i++) { 
      mean += phi[i]/length; 
      variance += phi[i]*phi[i]/length;   
    }
    stdev=sqrt(variance-mean*mean);
    for(i=0;i<length;i++) phi[i]=(phi[i]-mean)/stdev ;
    
    /* calculation of outer iteration criterion */
    oscreal=esqo;  /* used as memory of old esqo */
    esqo=0;
    for(i=0;i<length;i++) 
      esqo+=(sum[i]-phi[i])*(sum[i]-phi[i])/length;
    oscreal=fabs(oscreal-esqo);
    if (showoloop) 
      printf("itero=%5d yields e^2: %1.6lf (difference: %e)\n",
	     itero,esqo,oscreal);

    itero ++ ;            
  } /* end of outer loop */
  
  /* output */

  /* writing the phis in long form to file */
  for(d=0;d<=dim;d++) 
  {
      for(i=0;i<length;i++)
      {
       Table_Phi[j]=phi[d*length+i];
        j=j+1;
      }
   }
    
  
  

mxFree(phi);
mxFree(ind);
mxFree(iind);
mxFree(sum);
mxFree(yterm);
mxFree(cey);
mxFree(yind);

}

/*---------------------------------------------------------*/

void quicksort(double *vet, int *indvet, int s1, int s2)
{
   int i, j, a, sc, pasind;
   double pas, sep;

   sc = 0;
   i = s1;
   j = s2;
   if (vet[s1] != vet[s2])  sep = ((double)(vet[s1] + vet[s2]))/ 2;
   else 
   {
      sep = vet[s1];
      a = s1 + 1;
      while ((sep == vet[a]) && (a < s2)) a++;
      sep = ((double)(vet[s1] + vet[a])) / 2;
   }
   while (i < j)
   {
      while ((i <= j) && (vet[i] < sep)) i++;
      while ((i <= j) && (vet[j] >= sep)) j--;
      if (i < j)
      {
	 sc = 1;
	 pas = vet[i];
	 vet[i] = vet[j];
	 vet[j] = pas;
	 pasind = indvet[i];
	 indvet[i] = indvet[j];
	 indvet[j] = pasind;
      }
   }
   if ((sc != 0) || (i != s1))
   {
      if (s1 < j) quicksort(vet, indvet, s1, j);
      if (i < s2) quicksort(vet, indvet, i, s2);
   }
}

/*---------------------------------------------------------*/

int window(double x[], int length, int wl, double out[])
{
  int i, norm=2*wl+1;
  double sum=0;
  for(i=0;i<=2*wl;i++) sum+=x[i]; /* mean for the first point */
  for(i=wl;i<length-wl-1;i++){ /* cant count till the end because of sum */
    out[i]=sum/norm;
    sum=sum-x[i-wl]+x[i+wl+1];
  }
  out[length-wl-1]=sum/norm; /* mean for the last point */
  for(i=0;i<wl;i++){ out[i]=out[wl]; out[length-1-i]=out[length-1-wl]; };
/* for(i=0;i<length;i++) printf("x: %lf windowed x: %lf\n",x[i],r[i]); */
  return 0;
}

/*---------------------------------------------------------*/

int cef(double y[], int ind[], int iind[], int length, int wl,double out[],double* cey,double* yind)
{
  int i,n;
  //double *yind, *cey;
  //cey=(double *)calloc(length,sizeof(double)) ;
  //yind=(double *)calloc(length,sizeof(double)) ;

  for(i=0;i<length;i++) yind[i]=y[ind[i]];
  window(yind,length,wl,cey);
  for(i=0;i<length;i++) out[i]=cey[iind[i]];
  
  return 0;
}

/*---------------------------------------------------------*/

double corr(double *x, double *y, int length)
{
  int i;
  double  prod=0, mean1=0, mean2=0, variance1=0, variance2=0;
  for(i=0;i<length;i++) {mean1+=x[i]/length; mean2+=y[i]/length;}
  for(i=0;i<length;i++) {
    prod+=x[i]*y[i]/length;
    variance1 += ((x[i]- mean1) *( x[i]- mean1))/(length-1) ;
    variance2 += ((y[i]- mean2) *( y[i]- mean2))/(length-1) ;
  }
  return ((prod-mean1*mean2)/sqrt(variance1)/sqrt(variance2));
}

/* -------------------------------------------------------- */

void mexFunction( int nlhs, mxArray *plhs[],
                  int nrhs, const mxArray *prhs[] )
{
    
    double *phi, *sum, *yterm, *cey, *yind;
    int *ind, *iind;
    
    int Input_NRow,d, Input_NCol ;
    double *Table_Phi;
    double *x = mxGetPr(prhs[0]);
    double *Para = mxGetPr(prhs[1]); 
    
	if (nrhs != 2) mexErrMsgTxt("False Number of RHS Arguments.\n Usage: Phi = aceMex_fast(X,Parameter)");    
    if (nlhs != 1) mexErrMsgTxt("False Number of LHS Arguments.\n Usage: Phi = aceMex_fast(X,Parameter)");
    
    Input_NRow = mxGetM(prhs[0]); /* Number of Rows */
    Input_NCol = mxGetM(prhs[0]); /* Number of Column */
	  
 
    /* initializing fields */
    phi=(double*)mxCalloc((Input_NCol)*Input_NRow,sizeof(double));
    ind=(int*)mxCalloc((Input_NCol)*Input_NRow,sizeof(int));
    iind=(int*)mxCalloc((Input_NCol)*Input_NRow,sizeof(int)) ;
    sum=(double*)mxCalloc(Input_NRow,sizeof(double));
    yterm=(double*)mxCalloc(Input_NRow,sizeof(double));
    cey=(double *)mxCalloc(Input_NRow,sizeof(double)) ;
    yind=(double *)mxCalloc(Input_NRow,sizeof(double)) ;
    
    
    plhs[0] = mxCreateDoubleMatrix(Input_NRow, 1, mxREAL);
    Table_Phi = mxGetPr(plhs[0]);
    
    aceMEX(x,Para,Table_Phi,phi,ind,iind,sum,yterm,cey,yind);


}
