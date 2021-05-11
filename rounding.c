void rounding(double *x, double *f, double *w) {

	if (*x==0) {
		  *f = 0;
		  *w = 1;
	  }
	  else if (*x>0) {
		if ((*x - floor(*x))<=0.5) {
			*f = floor(*x);
			*w = *f + 1;
		}
		else {
			*f = ceil(*x);
			*w = *f - 1;
		}
	  }
	  else {
		if ((ceil(*x) - *x)<=0.5) {
			*f = ceil(*x);
			*w = *f - 1;
		}
		else {
			*f = floor(*x);
			*w = *f + 1;
		}
	  }
	
}