#include <R.h>
#include <Rdefines.h>

#include "entries.h"

static double correlation(double *px,double *py,int count,int rows)
// standard pearson; interleaved mem
{	double x_mean,y_mean;
	x_mean=y_mean=0;
	double *x=px;
	double *y=py;
	for (int i=0; i<count; ++i)
	{	x_mean+=*x;
		y_mean+=*y;
		x+=rows;
		y+=rows;
	}
	x_mean/=count;
	y_mean/=count;

	x=px;
	y=py;

	double cov,sx,sy;
	cov=sx=sy=0;
	for (int i=0; i<count; ++i)
	{	double xt=*x-x_mean;
		double yt=*y-y_mean;
		cov+=xt*yt;
		sx+=xt*xt;
		sy+=yt*yt;
		x+=rows;
		y+=rows;
	}
	double result=cov/sqrt(sx*sy);
	if (result!=result)	// we have NaN -> just return 0 ???
	{	result=0;
	}
	return result;
}

static double correlation_simple(double *px,double *py,int count)
// standard pearson; simple mem layout
{	double x_mean,y_mean;
	x_mean=y_mean=0;
	double *x=px;
	double *y=py;

//printf("x= ");for (int i=0; i<count; ++i) printf("%d ",(int)x[i]); printf("\n");
//printf("y= ");for (int i=0; i<count; ++i) printf("%d ",(int)y[i]); printf("\n");

	for (int i=0; i<count; ++i)
	{
		x_mean+=*(x++);
		y_mean+=*(y++);
	}
	x_mean/=count;
	y_mean/=count;
	x=px;
	y=py;
	double cov,sx,sy;
	cov=sx=sy=0;
	for (int i=0; i<count; ++i)
	{	double xt=*(x++)-x_mean;
		double yt=*(y++)-y_mean;
		cov+=xt*yt;
		sx+=xt*xt;
		sy+=yt*yt;
	}
	double result=cov/sqrt(sx*sy);
	if (result!=result)     // we have NaN -> just return 0 ???
	{       result=0;
	}
	return result;
}

static double correlation_maucher(double *px,double *py,int count)
// maucher et al. 2011 theroem 1; simple mem layout
{	double x_mean,y_mean;
	x_mean=y_mean=0;
	double *x=px;
	double *y=py;
	for (int i=0; i<count; ++i)
	{	x_mean+=*(x++);
		y_mean+=*(y++);
	}
	x_mean/=count;
	y_mean/=count;
	x=px;
	y=py;
	double cov,sx;
	cov=sx=0;
	for (int i=0; i<count; ++i)
	{	double xt=*(x++)-x_mean;
		double yt=*(y++)-y_mean;
		cov+=xt*yt;
		sx+=xt*xt;
	}
	double result=cov/sx;
	if (result!=result)
	{	result=0;
	}
	return (result);
}

static void dump_c_values(double *px,double *py,int count,int rows)
// debugging purpose; prints values to stdout
{	double *x=px;
	double *y=py;
	for (int i=0; i<count; ++i)
	{	printf("%f %f\n",*x,*y);
		x+=rows;
		y+=rows;
	}
}

SEXP InferViaCorrelation(SEXP data_,SEXP rows_,SEXP cols_,SEXP rownames_,SEXP treshold_)
{	int rows=*INTEGER(rows_);
	int cols=*INTEGER(cols_);
	double *data=REAL(data_);
	PROTECT(rownames_=AS_CHARACTER(rownames_));
	double treshold=*REAL(treshold_);
//printf("treshold %f\n",treshold);

	entry_t *results=NULL;

	for (int a=0; a<rows; ++a)
	{	for (int b=0; b<rows; ++b)
		{	double c=correlation(	data+a,
						data+b+rows,
						cols-1,
						rows);
			if (c>=treshold)
			{	results=entry_add(results,CHAR(STRING_ELT(rownames_,b)),CHAR(STRING_ELT(rownames_,a)));
			}
		}
	}
	// build named vector of vector(char) for R
	SEXP result;
	int ce=entry_count(results);
	PROTECT(result=allocVector(VECSXP,ce));
	SEXP list_names;
	PROTECT(list_names = allocVector(STRSXP, ce));

	entry_t *e=entry_first(results);
	for (int i=0; i<ce; ++i)
	{	SET_STRING_ELT(list_names, i,  mkChar(e->s));
		int count_d=dep_count(e->depends);
		int j=0;
		SEXP deps;
		deps=PROTECT(deps=NEW_CHARACTER(count_d));
		dep_t *d=dep_first(e->depends);
		while (d)
		{	SET_STRING_ELT(deps,j,mkChar(d->s));
			++j;
			d=dep_next(d);
		}
		SET_VECTOR_ELT(result,i,deps);
		e=entry_next(e);
	}
	setAttrib(result, R_NamesSymbol, list_names);

	UNPROTECT(3+ce);
	entry_free(results);
	return result;
}

SEXP ivc_multiple(SEXP data_,SEXP count_tr_,SEXP rows_,SEXP cols_,SEXP rownames_,SEXP treshold_,SEXP th_multi_)
// infer from multiple timeseries
{
	int count_tr=*INTEGER(count_tr_);
	int rows=*INTEGER(rows_);
	int cols=*INTEGER(cols_);
	double treshold=*REAL(treshold_);
	double th_multi=*REAL(th_multi_);
	th_multi*=count_tr;
	double *data=REAL(data_);
	PROTECT(rownames_=AS_CHARACTER(rownames_));

	entry_t *results=NULL;

	for (int ts=0; ts<count_tr; ++ts)
	{	double *dts=data+ts*rows*cols;
		for (int a=0; a<rows; ++a)
		{	for (int b=0; b<rows; ++b)
			{	double c=correlation(   dts+a,
							dts+b+rows,
							cols-1,
							rows);
				if (c>=treshold)
				{	results=entry_add(results,CHAR(STRING_ELT(rownames_,b)),CHAR(STRING_ELT(rownames_,a)));
				}
			}
		}
	}

	// return the results
	SEXP result;
	int ce=entry_count(results);
	PROTECT(result=allocVector(VECSXP,ce));
	SEXP list_names;
	PROTECT(list_names = allocVector(STRSXP, ce));

	entry_t *e=entry_first(results);
	for (int i=0; i<ce; ++i)
	{	SET_STRING_ELT(list_names, i,  mkChar(e->s));
		int count_d=dep_count_th(e->depends,th_multi);
		int j=0;
		SEXP deps;
		deps=PROTECT(deps=NEW_CHARACTER(count_d));
		dep_t *d=dep_first(e->depends);
		while (d)
		{	if (d->count>=th_multi)
			{	SET_STRING_ELT(deps,j,mkChar(d->s));
				++j;
			}
			d=dep_next(d);
		}
		SET_VECTOR_ELT(result,i,deps);
		e=entry_next(e);
	}
	setAttrib(result, R_NamesSymbol, list_names);

	UNPROTECT(3+ce);
	entry_free(results);
	return result;

}

int get_index(SEXP *all,int c,const char *what)
{	for (int i=0; i<c; ++i)
	{	if (strcmp(CHAR(STRING_ELT(*all,i)),what)==0)
		{	return i;
		}
	}
	return -1;
}

static void infer_function(int,double*,double*,int,int*,int);
static double	*corr=NULL;		// save last correlation coefficents used in ivc_ts
static int	corr_count_nodes=0;	// dito
static entry_t	*corr_names=NULL;	// the used node names
					// corr and corr_names will not freed if lib exits

static int *last_ordered=NULL;

SEXP ivc_tr_last_corr()
// returns last correlation values
{	SEXP	result;
	PROTECT(result = allocVector(VECSXP, corr_count_nodes));
	double	*ptr=corr;
	for (int i=0; i<corr_count_nodes; ++i)
	{	SEXP	s=PROTECT(allocVector(REALSXP,corr_count_nodes));
		for (int j=0; j<corr_count_nodes; ++j)
		{	REAL(s)[j]=*(ptr++);
		}
		if (i==0)	// add names to first row
		{	SEXP list_names;
			PROTECT(list_names = allocVector(STRSXP, entry_count(corr_names)));
			entry_t *e=entry_first(corr_names);
			for (int j=0; j<entry_count(corr_names); ++j)
			{	SET_STRING_ELT(list_names,j,mkChar(e->s));
				e=entry_next(e);
			}
			setAttrib(s, R_NamesSymbol, list_names);
		}
		SET_VECTOR_ELT(result,i,s);
	}

	UNPROTECT(corr_count_nodes+1+1);

	return result;
}

SEXP ivc_tr_last_ordered()
// returns ordered indices into correlation per row
{	SEXP	result;
	PROTECT(result = allocVector(VECSXP, corr_count_nodes));
	int	*ptr=last_ordered;
	for (int i=0; i<corr_count_nodes; ++i)
	{	SEXP	s=PROTECT(allocVector(INTSXP,corr_count_nodes));
		for (int j=0; j<corr_count_nodes; ++j)
		{	INTEGER(s)[j]=*(ptr++);
		}
		if (i==0)	// add names to first row
		{	SEXP list_names;
			PROTECT(list_names = allocVector(STRSXP, entry_count(corr_names)));
			entry_t *e=entry_first(corr_names);
			for (int j=0; j<entry_count(corr_names); ++j)
			{	SET_STRING_ELT(list_names,j,mkChar(e->s));
				e=entry_next(e);
			}
			setAttrib(s, R_NamesSymbol, list_names);
		}
		SET_VECTOR_ELT(result,i,s);
	}

	UNPROTECT(corr_count_nodes+1+1);

	return result;
}

SEXP ivc_tr_last_corr_ordered()
// retrieves ordered correlation values
{	SEXP	result;
	PROTECT(result = allocVector(VECSXP, corr_count_nodes));
	double	*ptr=corr;
	for (int i=0; i<corr_count_nodes; ++i)
	{	SEXP	s=PROTECT(allocVector(REALSXP,corr_count_nodes));
		int *lo=last_ordered+(i*corr_count_nodes);
		for (int j=0; j<corr_count_nodes; ++j)
		{	REAL(s)[j]=*(ptr+*lo);
			++lo;
		}
		if (i==0)	// add names to first row
		{	SEXP list_names;
			PROTECT(list_names = allocVector(STRSXP, entry_count(corr_names)));
			entry_t *e=entry_first(corr_names);
			for (int j=0; j<entry_count(corr_names); ++j)
			{	SET_STRING_ELT(list_names,j,mkChar(e->s));
				e=entry_next(e);
			}
			setAttrib(s, R_NamesSymbol, list_names);
		}
		SET_VECTOR_ELT(result,i,s);
		ptr+=corr_count_nodes;
	}

	UNPROTECT(corr_count_nodes+1+1);

	return result;
}

SEXP ivc_ts(SEXP data_,SEXP count_tr_,SEXP len_tr_,SEXP rownames_,SEXP treshold_,SEXP alt_)
// infer from state transitions
{	double *data=REAL(data_);
	int count_tr=*INTEGER(count_tr_);
	int len_tr=*INTEGER(len_tr_);
	int count_nodes=len_tr/2;
	double treshold=*REAL(treshold_);
	int alt=*INTEGER(alt_);
	PROTECT(rownames_=AS_CHARACTER(rownames_));
	double *x=(double*)malloc(sizeof(double)*count_tr*count_nodes);
	double *y=(double*)malloc(sizeof(double)*count_tr*count_nodes);

	// reordering data
	for (int i=0; i<count_tr; ++i)
	{	for (int j=0; j<count_nodes; ++j)
		{	x[j*count_tr+i]=data[i*len_tr+j];
			y[j*count_tr+i]=data[i*len_tr+j+count_nodes];
		}
	}

	entry_t *results=NULL;

	double *corr_ptr=NULL;
	corr_count_nodes=count_nodes;
	free(corr);
	corr=(double*)malloc(sizeof(double)*count_nodes*count_nodes);
	corr_ptr=corr;
	entry_free(corr_names);
	corr_names=NULL;

	for (int a=0; a<count_nodes; ++a)
	{	for (int b=0; b<count_nodes; ++b)
		{	double c;
			if (alt==1)
			{	c=correlation_simple(x+a*count_tr,y+b*count_tr,count_tr);
			}
			else
			{	c=fabs(correlation_maucher(x+a*count_tr,y+b*count_tr,count_tr));
			}
			if (c>=treshold)
			{	results=entry_add(results,CHAR(STRING_ELT(rownames_,b)),CHAR(STRING_ELT(rownames_,a)));
			}
			*(corr_ptr++)=c;
		}
//		corr_names=entry_add(corr_names,CHAR(STRING_ELT(rownames_,a)),NULL);
		R_CheckUserInterrupt();
	}

	free(x);
	free(y);

	SEXP result;
	int ce=entry_count(results);
	PROTECT(result=allocVector(VECSXP,ce));
	SEXP list_names;
	PROTECT(list_names = allocVector(STRSXP, ce));

	entry_t *e=entry_first(results);
	for (int i=0; i<ce; ++i)
	{	SET_STRING_ELT(list_names, i,  mkChar(e->s));
		int count_d=dep_count(e->depends);
		int *used=(int*)malloc(count_d*sizeof(int));
		int j=0;
		SEXP deps;
		deps=PROTECT(deps=NEW_CHARACTER(count_d));
		dep_t *d=dep_first(e->depends);
		while (d)
		{	SET_STRING_ELT(deps,j,mkChar(d->s));
			used[j]=get_index(&rownames_,count_nodes,d->s);
			++j;
			d=dep_next(d);
		}
		SET_VECTOR_ELT(result,i,deps);
		// infer function
		// find indices for used
		// infer_function(count_nodes,x,y,count_tr,used,count_d);
		free(used);
		// TODO add infered functions to result
		e=entry_next(e);
	}
	setAttrib(result, R_NamesSymbol, list_names);

	UNPROTECT(3+ce);
	entry_free(results);
	return result;
}

static int cmp_double(const void *a,const void *b)
{	return (double*)a<(double*)b;
}

SEXP ivc_transitions_order(SEXP data_,SEXP count_tr_,SEXP len_tr_,SEXP rownames_,SEXP treshold_,SEXP order_)
// as maucher et.al 2011 Infering Boolean function via higher-order correlations
// algorithm 1
{	double *data=REAL(data_);
	int count_tr=*INTEGER(count_tr_);
	int len_tr=*INTEGER(len_tr_);
	int count_nodes=len_tr/2;
	double treshold=*REAL(treshold_);
	PROTECT(rownames_=AS_CHARACTER(rownames_));
	double *x=(double*)malloc(sizeof(double)*count_tr*count_nodes);
	double *y=(double*)malloc(sizeof(double)*count_tr*count_nodes);
	int order=*INTEGER(order_);
	if (order>count_nodes)	// no higher order possible
	{	order=count_nodes;
	}

	// reordering data
	for (int i=0; i<count_tr; ++i)
	{	for (int j=0; j<count_nodes; ++j)
		{	x[j*count_tr+i]=data[i*len_tr+j];
			y[j*count_tr+i]=data[i*len_tr+j+count_nodes];
		}
	}

	entry_t *results=NULL;

	double *corr_ptr=NULL;
	corr_count_nodes=count_nodes;
	free(corr);
	corr=(double*)malloc(sizeof(double)*count_nodes*count_nodes);
	corr_ptr=corr;
	entry_free(corr_names);
	corr_names=NULL;

	for (int a=0; a<count_nodes; ++a)
	{	for (int b=0; b<count_nodes; ++b)
		{	double c;
			c=fabs(correlation_maucher(x+a*count_tr,y+b*count_tr,count_tr));
			if (c>=treshold)
			{	results=entry_add(results,CHAR(STRING_ELT(rownames_,b)),CHAR(STRING_ELT(rownames_,a)));
			}
			*(corr_ptr++)=c;
		}
		R_CheckUserInterrupt();
	}

	free(last_ordered);
	// init ordering of indices into corr
	last_ordered=(int*)malloc(sizeof(int)*count_nodes*count_nodes);
	for (int i=0; i<count_nodes*count_nodes; ++i)
	{	last_ordered[i]=i%count_nodes;
	}
	if (order>0)
	{	for (int a=0; a<count_nodes; ++a)
		{	int	again;
			do
			{	again=0;
				for (int b=0; b<count_nodes-1; ++b)
				{
					int i=a*count_nodes+b;
					double xl=corr[a*count_nodes+last_ordered[i]];
					double yl=corr[a*count_nodes+last_ordered[i+1]];
					if (xl<yl)
					{	again=1;
						int h=last_ordered[i];
						last_ordered[i]=last_ordered[i+1];
						last_ordered[i+1]=h;
					}
				}
			} while (again);
		}
		// use the order
		for (int o=0; o<order; ++o)
		{	for (int target=0; target<count_nodes; ++target)
			{	int i=target*count_nodes;
				int lo=last_ordered[i+o];
//printf("target=%d last_orderd=%d\n",target,lo);
				for (int j=0; j<count_nodes; ++j)
				{	if (j!=lo)
					{	double c=fabs(correlation_simple(x+target*count_tr,y+j*count_tr,count_tr));
						if (c>=treshold)
						{	results=entry_add(results,CHAR(STRING_ELT(rownames_,target)),CHAR(STRING_ELT(rownames_,j)));
//printf("added	%d %f\n",j,c);
						}
					}
				}
			}
		}

	}

	free(x);
	free(y);

	SEXP result;
	int ce=entry_count(results);
	PROTECT(result=allocVector(VECSXP,ce));
	SEXP list_names;
	PROTECT(list_names = allocVector(STRSXP, ce));

	entry_t *e=entry_first(results);
	for (int i=0; i<ce; ++i)
	{	SET_STRING_ELT(list_names, i,  mkChar(e->s));
		int count_d=dep_count(e->depends);
		int *used=(int*)malloc(count_d*sizeof(int));
		int j=0;
		SEXP deps;
		deps=PROTECT(deps=NEW_CHARACTER(count_d));
		dep_t *d=dep_first(e->depends);
		while (d)
		{	SET_STRING_ELT(deps,j,mkChar(d->s));
			used[j]=get_index(&rownames_,count_nodes,d->s);
			++j;
			d=dep_next(d);
		}
		SET_VECTOR_ELT(result,i,deps);
		// infer function
		// find indices for used
		// infer_function(count_nodes,x,y,count_tr,used,count_d);
		free(used);
		// TODO add infered functions to result
		e=entry_next(e);
	}
	setAttrib(result, R_NamesSymbol, list_names);

	UNPROTECT(3+ce);
	entry_free(results);
	return result;
}

/********************	walsh - hadamard - transformation ********************/

static int log_2(int x)
//	integer log2
{	if (x==0)
	{	return 0;
	}
	int ret=1;
	int n=2;
	while (x>n)
	{	ret++;
		n<<=1;
	}
	return ret;
}

static int pow_2(int x)
//	2 to power of x
{	return 1<<x;
}

void wht(double *x,int n)
//	compute walsh hadamard transformation
//	n must be power of 2
//	normalization is not done here (devide result by n)
{	int m=log2(n);
	double	*h=malloc(sizeof(double)*n);
	int offs=n/2;
	for (int step=0; step<m; ++step)
	{	for (int i=0; i<n; ++i)
		{	if (i%(2*offs)<offs)
			{       h[i]=x[i]+x[i+offs];
			}
			else
			{	h[i]=-x[i]+x[i-offs];
			}
		}
		offs/=2;
		memcpy(x,h,sizeof(double)*n);
	}
	free(h);
}

SEXP wht_r(SEXP x_,SEXP n_)
//	calling wht from R
{	double *x=REAL(x_);
	int n=*INTEGER(n_);
	// do copy; wht() changes it input
	// so input of this function would change in R
	double  *h=malloc(sizeof(double)*n);
	memcpy(h,x,sizeof(double)*n);
	wht(h,n);
	SEXP result;
	PROTECT(result=allocVector(REALSXP,n));
	memcpy(REAL(result),h,sizeof(double)*n);
	UNPROTECT(1);
	return result;
}

/******************** infer function ********************/

static void infer_function(int c_nodes,double *x,double *y,int c_tr,int *used,int c_used)
//	x .. input; y .. output states
//	[ node..0 .. node..n ] * count transistions
{
	if (c_used==0)
	{	return;
	}
	// make sure indices contiained in used are ascending and unique
	for (int i=0; i<c_used-1; ++i)
	{	if (used[i]>=used[i+1])
		{	fprintf(stderr,"ERROR : infer_function : contents of used have wrong order");
			return;
		}
	}
	// reduce to needed
	double *x2=(double*)malloc(c_tr*c_used*sizeof(double));
	double *y2=(double*)malloc(c_tr*c_used*sizeof(double));
	// i .. node; j .. used
	int j=0;
	for (int i=0; i<c_nodes && j<c_used; ++i)
	{	if (i==used[j])
		{	for (int k=0; k<c_tr; ++k)
			{	x2[k*c_used+j]=x[k*c_nodes+i];
				y2[k*c_used+j]=y[k*c_nodes+i];
			}
			++j;
		}
	}
	// find matching operation

	free(x2);
	free(y2);
	// TODO return result
}

