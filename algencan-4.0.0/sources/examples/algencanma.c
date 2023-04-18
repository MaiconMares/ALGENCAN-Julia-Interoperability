#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <limits.h>

typedef struct {
  int counters[5];
} pdata_type;

extern void algencan(void (*evalf)(int *n, double *x, double *f,
				   int *ierr, pdata_type *pdata),
		     void (*evalg)(int *n, double *x, double *g,
				   int *ierr, pdata_type *pdata),
		     void (*evalc)(int *n, double *x, int *m, int *p,
				   double *c, int *ierr, pdata_type *pdata),
		     void (*evalj)(int *n, double *x, int *m, int *p, int *ind,
				   int *jsorted, int *jsta, int *jlen, int *lim,
				   int *jvar, double *jval, int *ierr, pdata_type *pdata),
		     void (*evalhl)(int *n, double *x, int *m, int *p, double *lambda,
				    int *lim, int *inclf, int *hlnnz, int *hlrow, int *hlcol,
				    double *hlval, int *ierr, pdata_type *pdata),
		     int *jnnzmax, int *hlnnzmax, int *n, double *x, int *lind, double *lbnd, int *uind, double *ubnd,
		     int *m, int *p, double *lambda, double *epsfeas, double *epscompl, double *epsopt,
		     int *maxoutit, int *scale, int *rhoauto, double *rhoini, int *extallowed, int *corrin,
		     double *f, double *csupn, double *ssupn, double *nlpsupn, double *bdsvio, int *outiter,
		     int *totiter, int *nwcalls, int *nwtotit, int *ierr, int *istop, pdata_type *pdata);

void evalf (int *n, double *x, double *f, int *ierr, pdata_type *pdata) {
  *f = pow( x[0] + 4.0, 4.0 ) + pow( x[1], 2.0 );
  pdata->counters[0]++;
  printf ("EVALF C: pdata.counters[0] = %d\n", pdata->counters[0]);
}

void evalg (int *n, double *x, double *g, int *ierr, pdata_type *pdata) {
  g[0] = 4.0 * pow( x[0] + 4.0, 3.0 );
  g[1] = 2.0 * x[1];
  pdata->counters[1]++;
}

void evalc(int *n, double *x, int *m, int *p, double *c, int *ierr, pdata_type *pdata) {
  c[0] = pow( x[0], 3.0 ) - x[1] - 1.0;
  pdata->counters[2]++;
}

void evalj(int *n, double *x, int *m, int *p, int *ind,
	   int *jsorted, int *jsta, int *jlen, int *lim, int *jvar,
	   double *jval, int *ierr, pdata_type *pdata) {
  // Only gradients of constraints j such that ind(j) = .true. need
  // to be computed.
    
  if ( ind[0] ) {
    if ( *lim < *n ) {
      *ierr = -94;
      return;
    }

    jsta[0] = 1;
    jlen[0] = *n;

    for (int i = 0; i < *n; i++) jvar[i] = i+1;

    jval[0] = 3.0 * pow( x[0], 2.0 );
    jval[1] = - 1.0;

    // Says whether the variables' indices in jvar (related to this
    // constraint) are in increasing order. In case they are,
    // Algencan takes advantage of this. Implement sorted gradients
    // of constraints if you can do this in a natural (cheap)
    // way. Under no circumnstance use a sorting algorithm. (It is
    // better to set sorted(1) = .false. in this case.)

    jsorted[0] = 1;
  }
  pdata->counters[3]++;
}

void evalhl(int *n, double *x, int *m, int *p, double *lambda, int *lim,
	    int *inclf, int *hlnnz, int *hlrow, int *hlcol, double *hlval,
	    int *ierr, pdata_type *pdata) {
  *hlnnz = 0;

  // If .not. inclf then the Hessian of the objective function must not be included
  if ( *inclf ) {
    if ( *hlnnz + 2 > *lim ) {
      *ierr = -95;
      return;
    }

    hlrow[*hlnnz] = 1;
    hlcol[*hlnnz] = 1;
    hlval[*hlnnz] = 12.0 * pow( x[0] + 4.0, 2.0 );
    *hlnnz = *hlnnz + 1;
    
    hlrow[*hlnnz] = 2;
    hlcol[*hlnnz] = 2;
    hlval[*hlnnz] = 2.0;
    *hlnnz = *hlnnz + 1;
  }

  // Note that entries of the Hessian of the Lagrangian can be
  // repeated. If this is case, them sum of repeated entrances is
  // considered. This feature simplifies the construction of the
  // Hessian of the Lagrangian.
    
  if ( *hlnnz + 1 > *lim ) {
    *ierr = -95;
    return;
  }
    
  hlrow[*hlnnz] = 1;
  hlcol[*hlnnz] = 1;
  hlval[*hlnnz] = lambda[0] * 6.0 * x[0];
  *hlnnz = *hlnnz + 1;

  pdata->counters[4]++;
}
  
int main () {
  // LOCAL SCALARS
  int corrin,extallowed,rhoauto,scale;
  int allocerr,hlnnzmax,i,ierr,inform,istop,jnnzmax,m,maxoutit,
      n,nwcalls,nwtotit,outiter,p,totiter;
  double bdsvio,csupn,epsfeas,epscompl,epsopt,f,nlpsupn,rhoini,ssupn;
  clock_t finish,start;
  pdata_type pdata;
  
  // LOCAL ARRAYS
  int *lind,*uind;
  double *c,*lbnd,*ubnd,*lambda,*x;

  for (i = 0; i < 5; i++) pdata.counters[i] = 0;

  // Number of variables
  n = 2;

  x = malloc(n*sizeof(double));
  lind = malloc(n*sizeof(int));
  lbnd = malloc(n*sizeof(double));
  uind = malloc(n*sizeof(int));
  ubnd = malloc(n*sizeof(double));

  if ( x == NULL || lind == NULL || lbnd == NULL || uind == NULL || ubnd == NULL ) {
    printf("Allocation error.\n");
    return 1;
  }

  // Initial guess and bound constraints
  for (i = 0; i < n; i++) {
    x[i] = 10.0;

    lind[i] = 1;
    lbnd[i] = - 100.0;

    uind[i] = 1;
    ubnd[i] =   100.0;
  }

  // Number equality (m) and inequality (p) constraints.
  m = 0;
  p = 1;

  lambda = malloc((m+p)*sizeof(double));
  c = malloc((m+p)*sizeof(double));

  if ( c == NULL || lambda == NULL) {
    printf("Allocation error.\n");
    return 1;
  }

  // Initial guess for the Lagrange multipliers

  for (i = 0; i < m+p; i++) lambda[i] = 0.0;

  // Number of entries in the JAcobian of the constraints
 
  jnnzmax = n;

  // This should be the number of entries in the Hessian of the
  // Lagrangian. But, in fact, some extra space is need (to store the
  // Hessian of the Augmented Lagrangian, whose size is hard to
  // predict, and/or to store the Jacobian of the KKT system). Thus,
  // declare it as large as possible.
  
  hlnnzmax = INT_MAX / 2;

  // Feasibility, complementarity, and optimality tolerances
  
  epsfeas  = 1.0e-08;
  epscompl = 1.0e-08;
  epsopt   = 1.0e-08;

  // Maximum number of outer iterations
  
  maxoutit = 50;

  // rhoauto means that Algencan will automatically set the initial
  // value of the penalty parameter. If you set rhoauto = .false. then
  // you must set rhoini below with a meaningful value.
  rhoauto = 1;

  if ( !rhoauto ) rhoini = 1.0e-08;

  // scale = .true. means that you allow Algencan to automatically
  // scale the constraints. In any case, the feasibility tolerance
  // (epsfeas) will be always satisfied by the UNSCALED original
  // constraints.
  scale = 0;

  // extallowed = .true. means that you allow Gencan (the active-set
  // method used by Algencan to solve the bound-constrained
  // subproblems) to perform extrapolations. This strategy may use
  // extra evaluations of the objective function and the constraints
  // per iterations; but it uses to provide overal savings. You should
  // test both choices for the problem at hand.
  extallowed = 1;

  // extallowed = .true. means that you allow the inertia of the
  // Jacobian of the KKT system to be corrected during the acceleration
  // process. You should test both choices for the problem at hand.
  corrin = 0;

  start = clock();

  algencan(evalf,evalg,evalc,evalj,evalhl,&jnnzmax,&hlnnzmax,
	   &n,x,lind,lbnd,uind,ubnd,&m,&p,lambda,&epsfeas,&epscompl,
	   &epsopt,&maxoutit,&scale,&rhoauto,&rhoini,&extallowed,
	   &corrin,&f,&csupn,&ssupn,&nlpsupn,&bdsvio,&outiter,&totiter,
	   &nwcalls,&nwtotit,&ierr,&istop,&pdata);

  printf ("\n\nRetorno do Algencan, pdata.counters[0] = %d\n\n", pdata.counters[0]);

  finish = clock();

  printf("\n");
  printf("Number of variables                                   = %d\n",n);
  printf("Number of equality constraints                        = %d\n",m);
  printf("Number of inequality constraints                      = %d\n",p);

  printf("\n");
  printf("(REPORTED BY SOLVER) istop                            = %d\n",istop);
  printf("(REPORTED BY SOLVER) ierr                             = %d\n",ierr);
  printf("(REPORTED BY SOLVER) f                                = %e\n",f);
  printf("(REPORTED BY SOLVER) csupn                            = %e\n",csupn);
  printf("(REPORTED BY SOLVER) ssupn                            = %e\n",ssupn);
  printf("(REPORTED BY SOLVER) nlpsupn                          = %e\n",nlpsupn);
  printf("(REPORTED BY SOLVER) bounds violation                 = %e\n",bdsvio);
  printf("(REPORTED BY SOLVER) Number of outer iterations       = %d\n",outiter);
  printf("(REPORTED BY SOLVER) Number of inner iterations       = %d\n",totiter);
  printf("(REPORTED BY SOLVER) Number of Newton-KKT trials      = %d\n",nwcalls);
  printf("(REPORTED BY SOLVER) Number of Newton-KKT iterations  = %d\n",nwtotit);
  
  printf("\n");
  printf("(COMPUTED BY CALLER) Number of calls to evalf         = %d\n",pdata.counters[0]);
  printf("(COMPUTED BY CALLER) Number of calls to evalg         = %d\n",pdata.counters[1]);
  printf("(COMPUTED BY CALLER) Number of calls to evalc         = %d\n",pdata.counters[2]);
  printf("(COMPUTED BY CALLER) Number of calls to evalj         = %d\n",pdata.counters[3]);
  printf("(COMPUTED BY CALLER) Number of calls to evalhl        = %d\n",pdata.counters[4]);
  printf("(COMPUTED BY CALLER) CPU time in seconds              = %f\n",((float) (finish - start)) / CLOCKS_PER_SEC);

  // *****************************************************************
  // *****************************************************************
  // Just checking ...

  inform = 0;
  
  evalf(&n,x,&f,&inform,&pdata);
  
  if ( inform != 0 ) {
    printf("error when calling evalf in the main file.\n");
    return 1;
  }

  evalc(&n,x,&m,&p,c,&inform,&pdata);

  if ( inform != 0 ) {
    printf("error when calling evalc in the main file.\n");
    return 1;
  }

  csupn = 0.0;
  for (i = 0; i < m; i++)
    if (fabs(c[i]) > csupn) csupn = fabs(c[i]);
  for (i = m; i < m+p; i++)
    if (c[i] > csupn) csupn = c[i];

  bdsvio = 0.0;
  for (i = 0; i < n; i++) {
    if (lind[i])
      if (lbnd[i] - x[i] > bdsvio) bdsvio = lbnd[i] - x[i];

    if (uind[i])
      if (x[i] - ubnd[i] > bdsvio) bdsvio = x[i] - ubnd[i];
  }

  printf("\n");
  printf("(COMPUTED BY CALLER) f                                = %e\n",f);
  printf("(COMPUTED BY CALLER) csupn                            = %e\n",csupn);
  printf("(COMPUTED BY CALLER) bounds violation                 = %e\n\n",bdsvio);

  printf("When a quantity appears as computed by solver and computed by caller, they must coincide.\n");
  printf("(In case they do not coincide, please report it as a bug.)\n");
  // *****************************************************************
  // *****************************************************************
  
  free(lind);
  free(lbnd);
  free(uind);
  free(ubnd);
  free(x);
  free(lambda);
  free(c);

  return 0;
}
