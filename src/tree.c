#include <R.h>
#include <Rmath.h>
#include <math.h>

void con(int *Nf, int *Ny, double *fuTerm, double *SuTerm, double *con, double *result){
  int i, j, k;
  double matfs1, matfs2;
  for (i = 0; i < Ny[0]; i++) {
    for (j = 0; j < Nf[0]; j++) {
      for (k = 0; k < Nf[0]; k++) {
	matfs1 = fuTerm[i + j * Ny[0]] * SuTerm[i + k * Ny[0]] * con[i];
	matfs2 = fuTerm[i + k * Ny[0]] * SuTerm[i + j * Ny[0]] * con[i];
	if (matfs1 > matfs2) result[0] += matfs1;
	if (matfs1 == matfs2) result[0] += 0.5 * matfs1;
      }
    }
  }
}

void con3(int *Nf, int *Ny, double *fuTerm, double *SuTerm,
	  double *fuTerm2, double *SuTerm2, double *con, double *result){
  int i, j, k;
  double matfs1, matfs3, matfs4;
  for (i = 0; i < Ny[0]; i++) {
    for (j = 0; j < Nf[0]; j++) {
      for (k = 0; k < Nf[0]; k++) {
	matfs1 = fuTerm[i + j * Ny[0]] * SuTerm[i + k * Ny[0]] * con[i];
	// matfs2 = fuTerm[i + k * Ny[0]] * SuTerm[i + j * Ny[0]] * con[i];
	matfs3 = fuTerm2[i + j * Ny[0]] * SuTerm2[i + k * Ny[0]];
	matfs4 = fuTerm2[i + k * Ny[0]] * SuTerm2[i + j * Ny[0]];
	if (matfs3 > matfs4) result[0] += matfs1;
	if (matfs3 == matfs4) result[0] += 0.5 * matfs1;
      }
    }
  }
}

// Used in split
// Variables are:
// Integers: 
// n = length(id)
// cL = length(cutAll)
// m is the maximum node number
// y is an indicator function, indicating which ID satisfy Y <= tau * E
// Ny = sum(Y <= tau * E)
// Doubles:
// minsp and minsp2 are the range for the number of the minimum risk observations required to split.
// X is the pth covariate path, or xlist[[p]].
// ndInd is ndInd in R; indicator matrix tells where to possibly split. 
// cut is cutAll, values to cut the covaraite.
// fmat, Smat: \hat{f} and \hat{S}.
// con: denominator in CON
// f, S, r: constructed by \code{fmat} and \code{Smat}. Needed for computing CON.
// rm: rm, constructed by fm / Sm.
// Nm: number of potential terminals to split (index m in R's for loop).
// Returns CON at each cut point; result is a vector of length equal to length(cutAll).
void cutSearch(int *n, int *cL, int *m, int *y, int *Ny,
	       double *minsp, double *minsp2, double *X, double *ndInd,
	       double *cut, double *fmat, double *Smat, double *con,
	       double *f, double *S, double *r, double *rm, int *Nm, 
	       double *result) {
  int i, j, k;
  double cc, ssL, ssR;
  double *indccL = Calloc(n[0] * n[0], double);
  double *indccR = Calloc(n[0] * n[0], double);
  double *SmatL = Calloc(n[0] * n[0], double);
  double *SmatR = Calloc(n[0] * n[0], double);
  double *rsSmatL = Calloc(n[0], double);
  double *rsSmatR = Calloc(n[0], double);
  double minRS;
  double *fL = Calloc(Ny[0], double);
  double *fR = Calloc(Ny[0], double);
  double *SL = Calloc(Ny[0], double);
  double *SR = Calloc(Ny[0], double);
  double *rL = Calloc(Ny[0], double);
  double *rR = Calloc(Ny[0], double);
  double *c1 = Calloc(Ny[0], double);
  double *c2 = Calloc(Ny[0], double);
  for (i = 0; i < cL[0]; i++) {
    // Reset vectors to zero to start 
    for (j = 0; j < n[0]; j++){
      rsSmatL[j] = 0;
      rsSmatR[j] = 0;
    }
    for (j = 0; j < n[0] * n[0]; j++){
      indccL[j] = 0;
      indccR[j] = 0;
      SmatL[j] = 0;
      SmatR[j] = 0;
    }
    for (j = 0; j < Ny[0]; j++) {
      fL[j] = 0;
      fR[j] = 0;
      SL[j] = 0;
      SR[j] = 0;
      rL[j] = 0;
      rR[j] = 0;
      c1[j] = 0;
      c2[j] = 0;
    }
    // end reset zero
    cc = cut[i];
    for (j = 0; j < n[0] * n[0]; j++) {
      if (ndInd[j] == m[0] && X[j] <= cc) {
	indccL[j] = 1;
	SmatL[j] = Smat[j];
      }
      if (ndInd[j] == m[0] && X[j] > cc) {
	indccR[j] = 1;
	SmatR[j] = Smat[j];
      }
      rsSmatL[j % n[0]] += SmatL[j]; // n by 1
      rsSmatR[j % n[0]] += SmatR[j]; // n by 1
    } 
    ssL = 0.0;
    ssR = 0.0;
    minRS = minsp2[0] / 3;
    for (j = 0; j < Ny[0]; j++) {
      if (rsSmatL[y[j]] < minRS) minRS = 0; 
      if (rsSmatR[y[j]] < minRS) minRS = 0; 
      ssL += indccL[y[j] * (1 + n[0])];
      ssR += indccR[y[j] * (1 + n[0])];
    }
    if (ssL < minsp[0] / 3 || ssR < minsp[0] / 3) {
      result[i] = -1;
    } else {
      for (j = 0; j < n[0]; j++) {
	if (indccL[j * (n[0] + 1)] == 1) {
	  for (k = 0; k < Ny[0]; k++) {
	    fL[k] += fmat[j * n[0] + y[k]] / n[0];
	  }
	}
	if (indccR[j * (n[0] + 1)] == 1) {
	  for (k = 0; k < Ny[0]; k++) {
	    fR[k] += fmat[j * n[0] + y[k]] / n[0];
	  }
	}
      }
      for (j = 0; j < Ny[0]; j++) {
	SL[j] = rsSmatL[y[j]] / n[0];
	SR[j] = rsSmatR[y[j]] / n[0];
      }
      for (j = 0; j < Ny[0]; j++) {
	if (SL[j] == 0) rL[j] = 99999;
	if (SL[j] != 0) rL[j] = fL[j] / SL[j];
	if (SR[j] == 0) rR[j] = 99999;
	if (SR[j] != 0) rR[j] = fR[j] / SR[j];
      }
      if (m[0] == 1) {
	for (j = 0; j < Ny[0]; j++) {
	  result[i] += fabs(fR[j] * SL[j] - fL[j] * SR[j]) * con[j] / 2;
	}
      } else {
	for (j = 0; j < Ny[0] * Nm[0]; j++) {
	  if ((rm[j % Ny[0]] - r[j]) * (r[j] - rL[j % Ny[0]]) > 0)
	    c1[j % Ny[0]] += fabs(f[j] * SL[j % Ny[0]] - S[j] * fL[j % Ny[0]]);
	  if ((rm[j % Ny[0]] - r[j]) * (r[j] - rR[j % Ny[0]]) > 0)
	    c2[j % Ny[0]] += fabs(f[j] * SR[j % Ny[0]] - S[j] * fR[j % Ny[0]]);
	}
	for (j = 0; j < Ny[0]; j++) {
	  result[i] += con[j] * (fabs(fR[j] * SL[j] - fL[j] * SR[j]) / 2 + c1[j] + c2[j]);
	  //result[i] += c1[j];
	}
      }
    } // end else
  } // end i
  Free(indccL);
  Free(indccR);
  Free(SmatL);
  Free(SmatR);
  Free(rsSmatL);
  Free(rsSmatR);
  Free(SL);
  Free(SR);
  Free(fL);
  Free(fR);
  Free(rL);
  Free(rR);
  Free(c1);
  Free(c2);
}

// Used in split; search by \Delta CON
// Variables are the same as those in cutSearch
void cutSearch2(int *n, int *cL, int *m, int *y, int *Ny,
		double *minsp, double *minsp2, double *X, double *ndInd,
		double *cut, double *fmat, double *Smat, double *con,
		double *f, double *S, double *r, double *rm, int *Nm, 
		double *result) {
  int i, j, k;
  double cc, ssL, ssR;
  double *indccL = Calloc(n[0] * n[0], double);
  double *indccR = Calloc(n[0] * n[0], double);
  double *SmatL = Calloc(n[0] * n[0], double);
  double *SmatR = Calloc(n[0] * n[0], double);
  double *rsSmatL = Calloc(n[0], double);
  double *rsSmatR = Calloc(n[0], double);
  double minRS;
  double *fL = Calloc(Ny[0], double);
  double *fR = Calloc(Ny[0], double);
  double *SL = Calloc(Ny[0], double);
  double *SR = Calloc(Ny[0], double);
  double *rL = Calloc(Ny[0], double);
  double *rR = Calloc(Ny[0], double);
  for (i = 0; i < cL[0]; i++) {
    // Reset vectors to zero to start 
    for (j = 0; j < n[0]; j++){
      rsSmatL[j] = 0;
      rsSmatR[j] = 0;
    }
    for (j = 0; j < n[0] * n[0]; j++){
      indccL[j] = 0;
      indccR[j] = 0;
      SmatL[j] = 0;
      SmatR[j] = 0;
    }
    for (j = 0; j < Ny[0]; j++) {
      fL[j] = 0;
      fR[j] = 0;
      SL[j] = 0;
      SR[j] = 0;
    }
    // end reset zero
    cc = cut[i];
    for (j = 0; j < n[0] * n[0]; j++) {
      if (ndInd[j] == m[0] && X[j] <= cc) {
	indccL[j] = 1;
	SmatL[j] = Smat[j];
      }
      if (ndInd[j] == m[0] && X[j] > cc) {
	indccR[j] = 1;
	SmatR[j] = Smat[j];
      }
      rsSmatL[j % n[0]] += SmatL[j]; // n by 1
      rsSmatR[j % n[0]] += SmatR[j]; // n by 1
    } 
    ssL = 0.0;
    ssR = 0.0;
    minRS = minsp2[0] / 3;
    for (j = 0; j < Ny[0]; j++) {
      if (rsSmatL[y[j]] < minRS) minRS = 0; 
      if (rsSmatR[y[j]] < minRS) minRS = 0; 
      ssL += indccL[y[j] * (1 + n[0])];
      ssR += indccR[y[j] * (1 + n[0])];
    }
    if (minRS == 0 && (ssL < minsp[0] / 3 || ssR < minsp[0] / 3)) {
      result[i] = -1;
    } else {
      for (j = 0; j < n[0]; j++) {
	if (indccL[j * (n[0] + 1)] == 1) {
	  for (k = 0; k < Ny[0]; k++) {
	    fL[k] += fmat[j * n[0] + y[k]] / n[0];
	  }
	}
	if (indccR[j * (n[0] + 1)] == 1) {
	  for (k = 0; k < Ny[0]; k++) {
	    fR[k] += fmat[j * n[0] + y[k]] / n[0];
	  }
	}
      }
      for (j = 0; j < Ny[0]; j++) {
	SL[j] = rsSmatL[y[j]] / n[0];
	SR[j] = rsSmatR[y[j]] / n[0];
      }
      for (j = 0; j < Ny[0]; j++) {
	result[i] += fabs(fR[j] * SL[j] - fL[j] * SR[j]) * con[j];
      }
    } // end else
  } // end i
  Free(indccL);
  Free(indccR);
  Free(SmatL);
  Free(SmatR);
  Free(rsSmatL);
  Free(rsSmatR);
  Free(SL);
  Free(SR);
  Free(fL);
  Free(fR);
}

// C function to replace the "giveW" function in R.
void giveWC(int *n, int *lidB2, int *lndTerm,
	    int *ndi, int *idB2, int *ndInd2, int *ndTerm, double *szL2, double *result) {
  int i, j, k;
  for (i = 0; i < lidB2[0]; i++) {
    for (j = 0; j < n[0]; j++) {
      if (ndInd2[i * n[0] + j] == ndi[j]) {
    	for (k = 0; k < lndTerm[0]; k++) {
    	  if (ndTerm[k] == ndi[j] && szL2[k * lndTerm[0] + j] != 0) {
	    result[idB2[i] * n[0] + j] = 1 / szL2[k * lndTerm[0] + j];
    	  }
    	}
      }
    }
  }
}
