#include <R.h>
#include <Rinternals.h>
#include <string.h>
#include <Rmath.h>
#include <R_ext/Rdynload.h>
#include <R_ext/Applic.h>
#define EPSILON (1e-9)
static double *vector(int n)
{
  double *v;
  v = Calloc(n,double);
  return(v);
}
static int *vector_int(int n)
{
  int *v;
  v = Calloc(n,int);
  return(v);
}
static void free_vector(double *v)
{
  Free(v);
}
static void free_vector_int(int *v)
{
  Free(v);
}
static int sgn(double z)
{
  if (z > 0) return(1);
  else if (z < 0) return(-1);
  else return(0);
}
static double positive(double z)
{
  if (z > 0) return(z); else return(0);
}
static double soft(double lambda, double z)
{
  if (fabs(z) > lambda) return(sgn(z) * (fabs(z) - lambda));
  else return(0);
}
static double fmcp(double lambda, double kappa, double z)
{
  if (fabs(z) > lambda / kappa) return(z); else
  {
    if (z > lambda) return((z - lambda) / (1 - kappa));
    else if (z < -lambda) return ((z + lambda) / (1 - kappa));
    else return(0);
  }
}
static double fscad(double lambda, double kappa, double z)
{
  if (fabs(z) > lambda / kappa) return(z); else
  {
    if (fabs(z) > 2 * lambda) return(((1 / kappa - 1) * z - sgn(z) * lambda / kappa) / (1 / kappa - 2));
    else return(sgn(z) * positive(fabs(z) - lambda));
  }
}

static int equalZero(double num)
{
  if (fabs(num) < EPSILON) return(1); else return(0);
}

static int checkConvergence(double *beta, double *beta_old, double eps,
                            int len)
{
  int i, converged = 1;
  for (i = 0; i <= len; i++)
  {
    if (!equalZero(beta[i]) && !equalZero(beta_old[i]))
    {
      if (fabs((beta[i] - beta_old[i]) / beta_old[i]) > eps)
      {
        converged = 0;
        break;
      }
    }
    else if (!equalZero(beta[i]) && equalZero(beta_old[i]))
    {converged = 0;break;}
    else if (equalZero(beta[i]) && !equalZero(beta_old[i]))
    {converged = 0;break;}
  }
  return(converged);
}

static void cdfit(double *x, double *y, double *init, double lambda,
                  double kappa, char *penalty, double eps, int max, int n,
                  int p, int *iter)
{
  double *beta_old, *beta, *u;
  double  z;
  int i, j, itetime = 0, flag = 1;
  beta = vector(p);
  beta_old = vector(p);
  for (i = 0; i < p; i++)
  {
    beta_old[i] = init[i];
    beta[i] = init[i];
  }
  u = vector(n);
  for (i = 0; i < n; i++)
  {u[i] = y[i];
    for (j = 0; j < p; j++)
      u[i] = u[i] - x[i + j * n] * beta_old[j];
  }
  if (strcmp(penalty, "MCP") == 0)
  {
    while (flag && itetime <= max)
    {
      if (itetime > 0) for (i = 0; i < p; i++) beta_old[i] = beta[i];
      for (i = 0; i < p; i++)
      {
        z = 0;
        for (j = 0; j < n; j++) z = x[i * n + j] * u[j] + z;
        z = z + beta_old[i];
        beta[i] = fmcp(lambda, kappa, z);
        if (!equalZero(beta[i] - beta_old[i]))
          for (j = 0; j < n; j++)
            u[j] = u[j] - (beta[i] - beta_old[i]) * x[i * n + j];
      }
      itetime++;
      flag = !checkConvergence(beta, beta_old, eps, p);
    }
  }
  else if (strcmp(penalty, "LASSO") == 0)
  {
    while (flag && itetime <= max)
    {
      if (itetime > 0) for (i = 0; i < p; i++) beta_old[i] = beta[i];
      for (i = 0; i < p; i++)
      {
        z = 0;
        for (j = 0; j < n; j++) z = x[i * n + j] * u[j] + z;
        z = z + beta_old[i];
        beta[i] = soft(lambda, z);
        if (!equalZero(beta[i] - beta_old[i]))
          for (j = 0; j < n; j++)
            u[j] = u[j] - (beta[i] - beta_old[i]) * x[i * n + j];
      }
      itetime++;
      flag = !checkConvergence(beta, beta_old, eps, p);
    }
  }
  else if (strcmp(penalty, "SCAD") == 0)
  {
    while (flag && itetime <= max)
    {
      if (itetime > 0) for (i = 0;i < p;i++) beta_old[i] = beta[i];
      for (i = 0; i < p; i++)
      {
        z = 0;
        for (j = 0; j < n; j++) z = x[i * n + j] * u[j] + z;
        z = z + beta_old[i];
        beta[i] = fscad(lambda, kappa, z);
        if (!equalZero(beta[i] - beta_old[i]))
          for (j = 0; j < n; j++)
            u[j] = u[j] - (beta[i] - beta_old[i]) * x[i * n + j];
      }
      itetime++;
      flag = !checkConvergence(beta, beta_old, eps, p);
    }
  }
  else if (strcmp(penalty, "adaptive") == 0)
  {
    while (flag && itetime <= max)
    {
      if (itetime > 0) for (i = 0; i < p; i++) beta_old[i] = beta[i];
      for (i = 0; i < p; i++)
      {
        z = 0;
        for (j = 0; j < n; j++) z = x[i * n + j] * u[j] + z;
        z = z + beta_old[i];
        if (beta_old[i] == 0) beta[i] = 0; else 
          beta[i] = soft(lambda / pow(fabs(beta_old[i]), 1 / kappa), z);
        if (!equalZero(beta[i] - beta_old[i]))
          for (j = 0; j < n; j++)
            u[j] = u[j] - (beta[i] - beta_old[i]) * x[i * n + j];
      }
      itetime++;
      flag = !checkConvergence(beta, beta_old, eps, p);
    }
  }
  for (i = 0; i < p; i++) init[i] = beta[i];
  iter[0] = itetime;
  free_vector(beta);
  free_vector(beta_old);
  free_vector(u);
}

SEXP pathSearch(SEXP Rx, SEXP Ry, SEXP Rmethod, SEXP Rlambda, SEXP Rkappa,
                SEXP Rn_lambda, SEXP Rn_kappa, SEXP Reps, SEXP Rmax,
                SEXP Rpenalty)
{
  const char *method, *penalty;
  int  i, j, k, n, p, max, *mark, iter[1], n_lambda, n_kappa;
  double eps;
  double *lambda, *y, *x, *init, *kappa;
  void cdfit(), free_vec();
  SEXP betahat,  betatild, Riter, return_list;
  Rx = coerceVector(Rx, REALSXP);
  Ry = coerceVector(Ry, REALSXP);
  Rlambda = coerceVector(Rlambda, REALSXP);
  Rkappa = coerceVector(Rkappa, REALSXP);
  Rmethod = coerceVector(Rmethod, STRSXP);
  Rn_lambda = coerceVector(Rn_lambda, INTSXP);
  Rn_kappa = coerceVector(Rn_kappa, INTSXP);
  Reps = coerceVector(Reps, REALSXP);
  Rmax = coerceVector(Rmax, INTSXP);
  p = Rf_ncols(Rx);
  n = Rf_nrows(Rx);
  lambda = REAL(Rlambda);
  kappa = REAL(Rkappa);
  method = CHAR(STRING_ELT(Rmethod, 0));
  penalty = CHAR(STRING_ELT(Rpenalty, 0));
  n_lambda = INTEGER(Rn_lambda)[0];
  n_kappa = INTEGER(Rn_kappa)[0];
  eps = REAL(Reps)[0];
  max = INTEGER(Rmax)[0];
  x = REAL(Rx);
  y = REAL(Ry);
  init = vector(p);
  mark = vector_int(p);
  PROTECT(betahat = Rf_allocMatrix(REALSXP, p, n_lambda));
  PROTECT(betatild = Rf_allocMatrix(REALSXP, p, n_lambda * n_kappa));
  PROTECT(Riter = Rf_allocVector(INTSXP, n_lambda * n_kappa));
  PROTECT(return_list = Rf_allocVector(VECSXP, 5));
  for (i = 0;i < p; i++) init[i] = 0;
  /*LASSO path along lambda*/
  for (i = 0; i < n_lambda; i++)
  {
    cdfit(x, y, init, lambda[i], 0.0, "LASSO", eps, max, n, p, iter);
    for (j = 0; j < p; j++) REAL(betahat)[i * p + j] = init[j];
    INTEGER(Riter)[i] = iter[0];
  }
  
  SET_VECTOR_ELT(return_list, 0, betahat);
  /*MCP, SCAD or adaptive LASSO paths along kappa with n_lambda interim points*/
  if (strcmp(penalty, "MCP") == 0)
  {
    if (strcmp(method, "include") == 0)
    {
      for (i = 0; i < n_lambda; i++)
      {
        for (j = 0; j < p; j++) init[j] = REAL(betahat)[i * p + j];
        for (j = 0; j < n_kappa; j++)
        {
          cdfit(x, y, init, lambda[i], kappa[j], "MCP", eps, max, n,
                p, iter);
          for (k = 0; k < p; k++)
            REAL(betatild)[(i + j * n_lambda) * p + k] = init[k];
          INTEGER(Riter)[i + j * n_lambda] = iter[0];
        }
      }
    }
    else if (strcmp(method, "exclude") == 0)
    {
      for (i = 0; i < n_lambda; i++)
      {
        int pnew = 0, t = 0;
        double *xnew, *initnew;
        for (j = 0; j < p; j++)
          if (equalZero(REAL(betahat)[i * p + j])) {init[j] = 0; mark[j] = 0;}
          else {init[j] = REAL(betahat)[i * p + j]; mark[j] = 1; pnew++;}
          if (pnew > 0)
          {
            xnew = vector(n * pnew);
            initnew = vector(pnew);
            for (j = 0; j < p; j++)
              if (mark[j])
              {
                for (k = 0; k < n; k++)
                  xnew[t * n + k] = x[j * n + k];
                t++;
              }
              t = 0;
              for (j = 0; j < p; j++) if (mark[j]) {initnew[t] = init[j]; t++;}
              for (j = 0; j < n_kappa; j++)
              {cdfit(xnew, y, initnew, lambda[i], kappa[j], "MCP", eps,
                     max, n, pnew, iter);
                t = 0;
                for (k = 0; k < p; k++)
                  if (mark[k])
                  {REAL(betatild)[(i + j * n_lambda) * p + k] = initnew[t]; t++;}
                  else REAL(betatild)[(i + j * n_lambda) * p + k] = 0;
                  INTEGER(Riter)[i + j * n_lambda] = iter[0];}
              free(xnew);
              free(initnew);
          }
          else
          {
            for (j = 0; j < n_kappa; j++) {
              for (k = 0; k < p; k++)
                REAL(betatild)[(i + j * n_lambda) * p + k] = 0;
              INTEGER(Riter)[i + j * n_lambda] = 0;}
          }
      }
    }
  }
  if (strcmp(penalty, "SCAD") == 0)
  {
    if (strcmp(method, "include") == 0)
    {
      for (i = 0; i < n_lambda; i++)
      {
        for (j = 0; j < p; j++) init[j] = REAL(betahat)[i * p + j];
        for (j = 0; j < n_kappa; j++)
        {
          cdfit(x, y, init, lambda[i], kappa[j], "SCAD", eps, max, n,
                p, iter);
          for (k = 0; k < p; k++)
            REAL(betatild)[(i + j * n_lambda) * p + k] = init[k];
          INTEGER(Riter)[i + j * n_lambda] = iter[0];
        }
      }
    }
    else if (strcmp(method, "exclude") == 0)
    {
      for (i = 0; i < n_lambda; i++)
      {
        int pnew = 0, t = 0;
        double *xnew, *initnew;
        for (j = 0; j < p; j++)
          if (equalZero(REAL(betahat)[i * p + j])) {init[j] = 0; mark[j] = 0;}
          else {init[j] = REAL(betahat)[i * p + j]; mark[j] = 1; pnew++;}
          if (pnew > 0)
          {
            xnew = vector(n * pnew);
            initnew = vector(pnew);
            for (j = 0; j < p; j++)
              if (mark[j])
              {
                for (k = 0; k < n; k++)
                  xnew[t * n + k] = x[j * n + k];
                t++;
              }
              t = 0;
              for (j = 0; j < p; j++) if (mark[j]) {initnew[t] = init[j]; t++;}
              for (j = 0; j < n_kappa; j++)
              {cdfit(xnew, y, initnew, lambda[i], kappa[j], "SCAD", eps,
                     max, n, pnew, iter);
                t = 0;
                for (k = 0; k < p; k++)
                  if (mark[k])
                  {REAL(betatild)[(i + j * n_lambda) * p + k] = initnew[t]; t++;}
                  else REAL(betatild)[(i + j * n_lambda) * p + k] = 0;
                  INTEGER(Riter)[i + j * n_lambda] = iter[0];}
              free(xnew);
              free(initnew);
          }
          else
          {
            for (j = 0; j < n_kappa; j++) {
              for (k = 0; k < p; k++)
                REAL(betatild)[(i + j * n_lambda) * p + k] = 0;
              INTEGER(Riter)[i + j * n_lambda] = 0;}
          }
      }
    }
  }
  if (strcmp(penalty, "adaptive") == 0)
  {
    for (i = 0; i < n_lambda; i++)
    {
      int pnew = 0, t = 0;
      double *xnew, *initnew;
      for (j = 0; j < p; j++)
        if (equalZero(REAL(betahat)[i * p + j])) {init[j] = 0; mark[j] = 0;}
        else {init[j] = REAL(betahat)[i * p + j]; mark[j] = 1; pnew++;}
        if (pnew > 0)
        {
          xnew = vector(n * pnew);
          initnew = vector(pnew);
          for (j = 0; j < p; j++)
            if (mark[j])
            {
              for (k = 0; k < n; k++)
                xnew[t * n + k] = x[j * n + k];
              t++;
            }
            t = 0;
            for (j = 0; j < p; j++) if (mark[j]) {initnew[t] = init[j]; t++;}
            for (j = 0; j < n_kappa; j++)
            {cdfit(xnew, y, initnew, lambda[i], kappa[j], "adaptive", eps,
                   max, n, pnew, iter);
              t = 0;
              for (k = 0; k < p; k++)
                if (mark[k])
                {REAL(betatild)[(i + j * n_lambda) * p + k] = initnew[t]; t++;}
                else REAL(betatild)[(i + j * n_lambda) * p + k] = 0;
                INTEGER(Riter)[i + j * n_lambda] = iter[0];}
            free(xnew);
            free(initnew);
        }
        else
        {
          for (j = 0; j < n_kappa; j++) {
            for (k = 0; k < p; k++)
              REAL(betatild)[(i + j * n_lambda) * p + k] = 0;
            INTEGER(Riter)[i + j * n_lambda] = 0;}
        }
    }
  }
  SET_VECTOR_ELT(return_list, 1, betatild);
  SET_VECTOR_ELT(return_list, 2, Riter);
  SET_VECTOR_ELT(return_list, 3, Rlambda);
  SET_VECTOR_ELT(return_list, 4, Rkappa);
  UNPROTECT(4);
  free_vector(init);
  free_vector_int(mark);
  return(return_list);
}

static R_CMethodDef DotCEntries[] = {
  {"pathSearch", (DL_FUNC) &pathSearch, 9},
  {NULL}
};

void R_init_pareto(DllInfo *info)
{
  R_registerRoutines(info, DotCEntries, NULL, NULL, NULL);
}
