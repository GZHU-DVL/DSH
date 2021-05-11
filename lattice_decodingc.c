#include "mex.h"
#include <stdio.h>
#include <string.h>
#include <math.h>


void mexFunction (
  int nlhs, mxArray *plhs[],
  int nrhs, const mxArray *prhs[])
  
{

	double *x, *lattice, *dec_point, *aux, *aux_x, *auxh, *w, *auxc, *g, *coset2, *cosets, *result_coseti;
	int rows_x, columns_x, rows_l, columns_l, i, j, aux_index, indexmin; 
	double dist1, dist2;
	mxArray *result_dec, *result_coset2;
	mxArray *arguments[2];


    if ((nrhs!=3)&(nrhs!=2)) {
		mexErrMsgTxt("Sintax: lattice_decodingc(x, lattice) or lattice_decodingc(x, lattice, cosets)\n");
	  }

	/******************************************************************/
	  /* The dimensions of the arguments are obtained */
	  x = mxGetPr(prhs[0]);      
	  rows_x = mxGetM(prhs[0]);
	  columns_x = mxGetN(prhs[0]);
	  lattice = mxGetPr(prhs[1]);
	  rows_l = mxGetM(prhs[1]);
	  columns_l = mxGetN(prhs[1]);
      if (nrhs==3) {
        cosets = mxGetPr(prhs[2]);
      }

	  /* printf("Dimensions of x: %d x %d\n", rows_x, columns_x); */
	  
	  plhs[0] = mxCreateDoubleMatrix(rows_x, columns_x, mxREAL);
	  dec_point = mxGetPr(plhs[0]);   /* x rounded */
	  aux = mxCalloc(1, sizeof(double));
	  aux_x = mxCalloc(1, sizeof(double));
	  auxh = mxCalloc(4, sizeof(double));

	  switch ((int) *lattice) {
		case 0: 
			/* printf("Cubic lattice\n"); */
			for (i=0; i<columns_x; i++) {
				/* printf("%f\n", *(x+i)); */
				rounding(x+i, dec_point+i, aux);
				/* printf("%f, %f\n", *(dec_point+i), *aux); */
			}
			break;
		case 1:
			/* printf("Hexagonal lattice\n");
			 Here, the definition of the hexagonal lattice as the union of two 
     		  rectangular lattices is used */
			*aux_x = (*x)*sqrt(3)/3;
		    rounding(aux_x, auxh, aux);
			*auxh = (*auxh)*3/sqrt(3);
			rounding(x+1, auxh+1, aux);
			*aux_x = (*(x) - 3/(2*sqrt(3)))*sqrt(3)/3;
			rounding(aux_x, auxh+2, aux);
			*(auxh+2) = (*(auxh+2))*3/sqrt(3) + 3/(2*sqrt(3));
			*aux_x = *(x+1) - 0.5;
            rounding(aux_x, auxh+3, aux);
			*(auxh+3) = *(auxh+3) + 0.5;
			/* coset distance computation */
			dist1 = pow(*x - *auxh, 2) + pow(*(x+1) - *(auxh+1), 2);
			dist2 = pow(*x - *(auxh+2), 2) + pow(*(x+1) - *(auxh+3), 2);
	    	if (dist1<dist2) {
				*dec_point = *auxh;
				*(dec_point + 1) = *(auxh+1);
			}
			else {
				*dec_point = *(auxh+2);
				*(dec_point + 1) = *(auxh+3);
			}
			break; 
		case 2:
			/* printf("Spherical Voronoi cell\n"); */
			dist1 = 0;
			for (i=0; i<columns_x; i++) {
				rounding(x+i, dec_point+i, aux);
				dist1 = dist1 + pow(*(x+i) - *(dec_point + i), 2);
			}
			if (dist1 > 0.25) {
				*(dec_point) = 10000;	/* this is for indicating x out of any sphere */
			}
			break;
		case 3: 
			/* printf("Checkerboard lattice Dn (n>=3)\n"); */
			w = mxCalloc(columns_x, sizeof(double));
			auxc = mxCalloc(columns_x, sizeof(double));
			g = mxCalloc(columns_x, sizeof(double));
			for (i=0; i<columns_x; i++) {
				/* printf("%f\n", *(x+i)); */
				rounding(x+i, dec_point+i, w+i);
				*(g+i) = *(dec_point+i);
				*(auxc+i) = fabs(*(dec_point+i) - *(x+i));
				/* printf("%f, %f\n", *(f+i), *(w+i)); */
			}
			/* the components of the decoded point must have an even sum */
			dist1 = 0;
			for (i=0; i<columns_x; i++) {
				dist1 = dist1 + *(dec_point+i);
			}
			if (fmod(dist1,2)!=0) {
				/* in this case, dec_point must be equal to g */
				aux_index = 0;	/* initialization of the maximum */
				for (i=1; i<columns_x; i++)	{
					if (*(auxc+i)>=*(auxc + aux_index)) {
						aux_index = i;
					}
				}
				*(dec_point + aux_index) = *(w + aux_index);
			}
			/* if fmod(dist1,2)==0, dec_point is equal to f */
			break;
			
		case 4: 
			/* printf("Lattice Dn* (dual of Dn)\n"); */
			/* Here, the definition of Dn* as the union of of two cosets of Z^n is used */
			coset2 = mxCalloc(columns_x, sizeof(double));
			auxc =  mxCalloc(columns_x, sizeof(double));
			dist1=0;
			for (i=0; i<columns_x; i++) {
				/* printf("%f\n", *(x+i)); */
				rounding(x+i, dec_point+i, aux);
				dist1 = dist1 + pow(*(x+i) - *(dec_point+i), 2);
			}
			dist2=0;
			for (i=0; i<columns_x; i++) {
				*(auxc+i) = *(x+i) - 0.5;
				rounding(auxc+i, coset2+i, aux);
				*(coset2+i) = *(coset2+i) + 0.5;
				dist2 = dist2 + pow(*(x+i) - *(coset2+i), 2);
			}
			if (dist1>dist2) {
				memcpy(dec_point, coset2, columns_x*sizeof(double));
			}
			break;
		
		case 5:
			/* printf("Gosset lattice E8 \n"); */
            /* Here, the definition of E8 as the union of two cosets of D8 is used */
			arguments[0] = mxCreateDoubleMatrix(1, columns_x, mxREAL);	/* vector to be quantized */
			arguments[1] = mxCreateDoubleMatrix(1, 1, mxREAL);	/* lattice type */
            memcpy(mxGetPr(arguments[0]), x, columns_x*sizeof(double));
            *(mxGetPr(arguments[1])) = 3.0;
			mexCallMATLAB(1, &result_dec, 2, arguments, "lattice_decodingc");   /* decoding with coset [0 0 ... 0] */
            auxc =  mxCalloc(columns_x, sizeof(double));
			for (i=0; i<columns_x; i++) {
                *(auxc+i) = *(x+i) - 0.5;
				/* printf("%f\n", *(mxGetPr(result_dec)+i)); */
			}
            memcpy(mxGetPr(arguments[0]), auxc, columns_x*sizeof(double));
            mexCallMATLAB(1, &result_coset2, 2, arguments, "lattice_decodingc");   /* decoding with coset [1/2 1/2 ... 1/2] */
            dist1 = 0;
            dist2 = 0;                        
            for (i=0; i<columns_x; i++) {
                dist1 = dist1 + pow(*(x+i) - *(mxGetPr(result_dec)+i), 2);
                *(mxGetPr(result_coset2) + i) = *(mxGetPr(result_coset2)+ i) + 0.5; 
                /* printf("%f ", *(mxGetPr(result_coset2) + i)); */
                dist2 = dist2 + pow(*(x+i) - *(mxGetPr(result_coset2) + i), 2); 
            }            
            if (dist1 < dist2) {
                memcpy(dec_point, mxGetPr(result_dec), columns_x*sizeof(double));
            } 
            else {
                memcpy(dec_point, mxGetPr(result_coset2), columns_x*sizeof(double));
            }
            break;
            
        case 6:
            /* printf("Dual of lattice E7 \n"); */
            /* Here, the definition of E7* as the union of 16 cosets of 2Z^7 is used */
            auxc =  mxCalloc(columns_x, sizeof(double));
            w = mxCalloc(1, sizeof(double));
            result_coseti = mxCalloc(columns_x, sizeof(double));
            dist1 = 10000;
            indexmin = 0;
            for (j=0; j<16; j++) {
                /* printf("Decoding with coset %d\n", j+1); */
                dist2 = 0;
                for (i=0; i<columns_x; i++) {                    
                    *(auxc + i) = (*(x+i)/sqrt(0.5) - *(cosets + j*columns_x + i))/2.0;                    
                    rounding(auxc + i, result_coseti + i, w);
                    *(result_coseti + i) = (*(result_coseti + i))*2.0 + *(cosets + j*columns_x + i);
                    dist2 = dist2 + pow(*(x + i)/sqrt(0.5) - *(result_coseti + i), 2);
                }                                
                /* printf("dist2 = %f\n", dist2); */
                if (dist2 <= dist1) {
                    dist1 = dist2;
                    indexmin = j;
                    memcpy(dec_point, result_coseti, columns_x*sizeof(double));
                }
            }
            /* printf("indexmin = %d", indexmin+1); */
            for (i=0; i<columns_x; i++) {
                    *(dec_point + i) = (*(dec_point + i))*sqrt(0.5);
            }
            
	  }

	  
}



