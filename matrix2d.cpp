/*
File: MATRIX2D.C
Basic 2D-matrix implementation in C++.
Gustavo Rohde, summer 2000.
*/


#include "matrix2d.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

/*
 * Allocate memory for a two-dimensional RxC matrix.
 * Returns 0 on failure.
 */
void matrix_2d::init(int x, int y)
{ 
  ny= y;
  nx= x;
  
  //En_Trg = 0;
  ntot=(JHtype)(ny*nx);

  // m = (JHtype *)malloc(ntot*sizeof(JHtype));
  int i,j;
  m =(JHtype **)malloc(nx*sizeof(JHtype*));
  for(i=0;i<nx;i++){
    
    m[i] = (JHtype *)malloc(ny*sizeof(JHtype));
    
  }
  
  tmp_src_hist = (JHtype *)malloc(nx*sizeof(JHtype));
  tmp_trg_hist = (JHtype *)malloc(ny*sizeof(JHtype));
  
  nl2 = (data_type)log(double(2));
  if ((m==NULL)||(tmp_src_hist==NULL)||(tmp_src_hist==NULL)){
    exit(1);
  }else{
    setzero();
  }
  
    
}


/*
 * Destroy 
 */
void matrix_2d::destroy()
{

  if ((m==NULL)||(tmp_src_hist==NULL)||(tmp_src_hist==NULL)){
	  printf("here is the bug\n");
	  exit(1);
  }else{
	
    //free(m);
    
    int i,j;
 
    for (i=0; i<nx; i++) {
	free(m[i]);
    }

    free(m);
    
    free(tmp_trg_hist);
    free(tmp_src_hist);
  }

}


/*
 * Set all matrix elements to zero.
 */
void matrix_2d::setzero()
{

  JHtype i;
  JHtype j;

  for(j=0;j<nx;j++)
	  tmp_src_hist[j]=0;

   for(j=0;j<ny;j++)
	  tmp_trg_hist[j]=0;

   for(i=0;i<nx;i++){
     for(j=0;j<ny;j++){
	  m[i][j]=0;
     }
   }
  
}


/*
 * computes the mutual information of the array
 */
data_type2 matrix_2d::mi(void){

  int i,j;
  data_type psrc;
  data_type ptrg;
  data_type pjh;
  data_type result = 0;
  data_type Hsrc = 0;
  data_type Htrg = 0;
  data_type JH = 0;
  JHtype total_src = 0;
  JHtype total_trg = 0;
  JHtype total_jh;
  
  //zero the temporary source image histogram
  for(i=0;i<nx;i++)
    tmp_src_hist[i] = 0;
  
  //zero the temporary target image histogram
  for(i=0;i<ny;i++)
    tmp_trg_hist[i] = 0;
  
  //collapse data of src part to calculate the entropy of the instantaneous source image
  for(i=0;i<nx;i++){
    for(j=0;j<ny;j++){
      tmp_src_hist[i] = tmp_src_hist[i]+get(i,j);
    }
  }
  
  //collapse all of trg part to calculate the entropy of the target image
  for(j=0;j<ny;j++){
    for(i=0;i<nx;i++){
      tmp_trg_hist[j] = tmp_trg_hist[j]+get(i,j);
    }
  }
  
  //calculate total of the tmp_src_hist
  for(i=0;i<nx;i++)
    total_src = tmp_src_hist[i] + total_src;
  
  //calculate total of the tmp_trg_hist
  for(i=0;i<ny;i++)
    total_trg = tmp_trg_hist[i] + total_trg;
  
  
  //compute the entropy of the source image
  for(i=0;i<ny;i++){
    psrc = ((data_type)tmp_src_hist[i])/((data_type)total_src);
    if(psrc>0)
      Hsrc = Hsrc - psrc*(((data_type)log(psrc))/nl2);
  }
  
  
  //compute the entropy of the target image
  for(i=0;i<ny;i++){
    ptrg = ((data_type)tmp_trg_hist[i])/((data_type)total_trg);
    if(ptrg>0)
      Htrg = Htrg - ptrg*(((data_type)log(ptrg))/nl2);
  }
  
  //compute the joint entropy
  total_jh = sum();
  
  for(i=0;i<nx;i++){
    for(j=0;j<ny;j++){
      pjh = ((data_type)get(i,j))/total_jh;
      if(pjh>0)
	JH = JH - pjh*(((data_type)log(pjh))/nl2);
    }
  }
  
  //JH = -JH;
  //Htrg = -Htrg;
  //Hsrc = -Hsrc;
  //NOTE return the negative of the mutual information for minimization purposes
  //return -(Hsrc+Htrg-JH);
  return -( (Hsrc+Htrg)/JH );
}

/*
 * computes the entropy of the trg image, needs to be called only once
 */
/*
void matrix_2d::trg_en_cpt(void){

	//ny contains the number of bins for target image histogram
	//the second argument is the one that refers to the trg data

	int i,j;
	JHtype *array;
	JHtype total = 0;
	data_type ptrg;
	data_type Htrg = 0;
	

	//allocate memory and set it to zero
	array = (JHtype *) malloc(ny*sizeof(JHtype));
	for(i=0;i<ny;i++)
		array[i] = 0;

	//collapse all the values into the array
	for(j=0;j<ny;j++){
		for(i=0;i<nx;i++){
			array[j] = array[j]+get(i,j);
		}
	}

	//compute the total of this array
	for(i=0;i<ny;i++)
		total = total + array[i];

	//compute the entropy of the target image
	for(i=0;i<ny;i++){
		ptrg = ((data_type)array[i])/((data_type)total);
		if(ptrg>0)
			Htrg = Htrg + ptrg*(((data_type)log(ptrg))/nl2);
	}
	
	En_Trg = Htrg;
	free(array);
}
*/

/*
 * computes the total of the array
 */
JHtype matrix_2d::sum(void){

  int i,j;
  JHtype result = 0;
  
  for(i=0;i<nx;i++){
    for(j=0;j<ny;j++){
      result = result + m[i][j];
    }
  }
  
  return result;
  
}


/*
 * returns the value stored in indexes i, j
 */
JHtype matrix_2d::get(int i, int j){

  //JHtype kk;

  //kk = j*nx+i;
	
  //return m[kk];

  return m[i][j];

}


/*
 * set the value in indexes i, j to val
 */
void matrix_2d::set(int i, int j, JHtype val){

  //JHtype kk;

  //kk = j*nx+i;
  //m[kk]=val;

  m[i][j]=val;
	
}


/*
 * increment the existing value at i, j by 1
 */
void matrix_2d::increment(int i, int j){

  //JHtype kk;

  //kk = j*nx+i;
  //m[kk]=m[kk]+1;
  m[i][j]+=1;
	
}


/*
 *subtract the existing value at i, j by val
 */
void matrix_2d::subtract(int i, int j, JHtype val){

  JHtype kk;
  
  //kk = j*nx+i;
  //m[kk]=m[kk]-val;
  //if(m[kk]<0)
  //	printf("negative here\n");
  m[i][j]-=val;
  
}


/*
 *add the existing value at i, j by val
 */
void matrix_2d::add(int i, int j, JHtype val){

  //JHtype kk;

  //kk = j*nx+i;
  //m[kk]=m[kk]+val;

  m[i][j]+=val;
	
}
