/*--------------------------------------------------------------------
  Title: A simulation library for QAOA on sparse, "Grover"-type binary 
  constraint satisfaction problems.
  Author: Aniruddha Bapat
  
  How to cite: 
  -------------------------------------------------------------------*/


#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "grover.h"

// Read a search instance from file and initialize search instance
void readInstance(char *path,searchInstance *S){
  // The instance file should be formatted like a cnf instance file
  // lines starting with 'c' for comments, 'p' for parameters
  // Parameter line format: p n m isGrover isHom
  // Each line of data:     string cost
  // Instance variables
  FILE *cnf;
  int    n; 
  int    m;
  int    isGrover;
  int    isHom;
  
  // Counters, buffers, and other useful constants
  int    i, nstring, ind; 

  // Number of bits cannot exceed buffer length
  int    lenbuff=512;
  char   buff[lenbuff];
  char   trash[lenbuff];
  
  cnf = fopen(path, "r");
  nstring=0;

  while(!feof(cnf)){
    // FIRST CHARACTER DETERMINES IF LINE IS DATA OR NOT
    fscanf(cnf, "%s", buff);
    
    // PARSER FOR COMMENTS, META-DATA, AND TRASH
    if(buff[0]=='\n'){
      continue;
    }
    else if(buff[0]==' '){
      continue;
    }
    else if(buff[0]=='c'){
      fgets(trash, lenbuff, cnf);
      continue;
    }
    else if(buff[0]=='p'){
      fgets(buff, lenbuff, cnf);
      sscanf(buff, "%d %d %d %d\n", &n, &m, &isGrover, &isHom);

      S->n         = n;
      S->m         = m;
      S->isGrover  = isGrover;
      S->isHom     = isHom;
      S->locs      = (int **)malloc(m*sizeof(int *));
      for (i=0;i<m;i++){
	S->locs[i] = (int *)malloc(n*sizeof(int));
      }
      S->dist      = (int *)malloc(m*m*sizeof(int));
      S->cost      = (double *)malloc(m*sizeof(double));
      S->combs     = (double *)malloc((n+1)*sizeof(double));

      fprintf(stderr, "I read a search instance on %d bits with %d marked item bitstrings.\n", n, m);
      continue;
    }

    // PARSER FOR THE MARKED ITEM LOCATIONS AND COSTS    

    // Don't read more than m strings
    if (nstring >= m){
      fprintf(stderr, "Warning: Too many marked strings, or trailing empty lines in file. I only read the first %d string, cost pairs.\n", m);
      break;
    }
    // Grover instance, all energies = -1
    for(i=0;i<n;i++){
      if      (buff[i]=='0') S->locs[nstring][i]=0;
      else if (buff[i]=='1') S->locs[nstring][i]=1;
      else fprintf(stderr, "Error: String %d is not in binary. Please fix it before proceeding!\n", nstring+1);
    }
    (isGrover)? S->cost[nstring] = -1.0 : fscanf(cnf, "%lf", &S->cost[nstring]);
    fgets(trash, lenbuff, cnf);
    nstring++;
    continue;  
  }
  
  // Fill in remaining fields in the instance
  if (isGrover) S->min = -1;
  else {
    ind = 0;
    for(i=0;i<m;i++){
      ind = (S->cost[ind] < S->cost[i])? ind:i;
    }
    S->min=ind;
  }
  measureDist(S);
  combinate(S);
  // Close the file 
  fclose(cnf);
}

// Read a QAOA schedule from file
void readSchedule(char *schd, schedule *U){
  // The schedule file should be formatted like a cnf instance file
  // lines starting with 'c' for comments, 'p' for parameters
  // Parameter line format: p T t isCyclic 
  // Each line of data:     beta-i gamma-i
  // Instance variables
  FILE *cnf;
  double *beta, *gamma;
  int    T, t, isCyclic;
  
  // Counters, buffers, and other useful constants
  int    i, nangle; 
  int    lenbuff=512;
  char   buff[lenbuff];
  char   trash[lenbuff];
  double dbuff[lenbuff];
  
  cnf = fopen(schd, "r");
  nangle=0;

  while(!feof(cnf)){
    // FIRST CHARACTER DETERMINES IF LINE IS DATA OR NOT
    fscanf(cnf, "%s", buff);
    
    // PARSER FOR COMMENTS, META-DATA, AND TRASH
    if(buff[0]=='\n'){
      continue;
    }
    else if(buff[0]==' '){
      continue;
    }
    else if(buff[0]=='c'){
      fgets(trash, lenbuff, cnf);
      continue;
    }
    else if(buff[0]=='p'){
      fgets(buff, lenbuff, cnf);
      sscanf(buff, "%d %d %d\n", &T, &t, &isCyclic);

      U->T         = T;
      U->t         = (isCyclic)? t:1;
      U->isCyclic  = isCyclic;
      // The number of (beta,gamma) pairs L depends on whether the protocol
      // is cyclic or not.
      U->L       = (isCyclic)? t:T;
      U->beta    = (double *)malloc(U->L*sizeof(double));
      U->gamma   = (double *)malloc(U->L*sizeof(double));

      if (isCyclic){
	fprintf(stderr, "I read a cyclic QAOA schedule with %d cycles and %d rounds per cycle.\n"\
		, T, t);
      }
      else fprintf(stderr, "I read an acyclic QAOA schedule with %d rounds.\n", T);
      continue;
    }

    // NOW, WE READ THE ANGLES
    
    // Don't read more than L angles
    if (nangle >= U->L){
      fprintf(stderr, "Warning: Too many angles, or trailing empty lines in file. I only read the first %d pairs of angles.\n", U->L);
      break;
    }
    
    // Read angles in the current line. buff should hold beta right now. Let's first convert
    // it to a double.
    sscanf(buff, "%lf", &U->beta[nangle]);
    fscanf(cnf, "%lf", &U->gamma[nangle]);
    // Scan the rest of the line and trash it
    fgets(trash, lenbuff, cnf);
    nangle++;
  }  
  fclose(cnf);
}

// Measure dist from locs
void measureDist(searchInstance *S) {
  int i,j,a;
  int d = 0;
  int n = S->n;
  int m = S->m;
  for(i=0;i<m;i++){
    for(j=0;j<m;j++){
      for(a=0;a<n;a++){
	d += (S->locs[i][a])^(S->locs[j][a]);
      }
      S->dist[m*i+j] = d;
      d=0;
    }
  }
}

// N choose w
double choose(int n, int w){
  int mini = (w < (n-w))? w:(n-w);
  double prod=1;
  int v;
  if (mini==0){
    return 1;
  }
  for(v=0;v<mini; v++){
    prod *= (double)(n-v)/(double)(mini-v);
  }
  return prod;
}

// A customized multinomial coefficient for this setting
// NOTE: ROW MAJOR FORMAT: row n, column k is n choose k
void multinominate(int n, double *mnoms){
  int i,j;
  for(i=0;i<=n;i++){
    for(j=i;j<=n;j++){
      mnoms[i*(n+1)+j] = 0;
      mnoms[j*(n+1)+i] = choose(j,i);
    }
  }
}

// Store combinatorial terms sqrt(n choose w) in an array
void combinate(searchInstance *S){
  int w;
  int n = S->n;
  for(w=0;w<=n;w++){
    S->combs[w]=sqrt(choose(n,w));
  }
}

// Prepare the initial state
void prepareState(double *state, searchInstance *S){
  int m = S->m;
  int n = S->n;
  int k = m; // Number of displaced sectors the initial state has support on
  int w,i;
  double amp;
  double norm;
  norm = pow(sqrt(2),n);
  for(w=0;w<=n;w++){
    amp = S->combs[w]/(k*norm);
    for(i=0;i<k;i++){
      state[2*((n+1)*i+w)]   = amp;
      state[2*((n+1)*i+w)+1] = 0;
    }
  }
}

// Diffusion for angle beta on the Hamming subspace
// Main idea: diagonalize D (which is tridiagonal), then
// exponentiate and undiagonalize
// Fortran matrices are column major, but this particular matrix
// is symmetric. An intermediate matrix in this routine isn't, however. 
void diffuse(double beta, searchInstance *S, double *Diff){
  int n = S->n;
  int np1 = n+1;
  int w,v;
  int lwork = 3*n;
  char compz = 'I';
  int info;

  double * offdiag   = malloc(n*sizeof(double));
  double * diag      = malloc(np1*sizeof(double));
  double * workspace = malloc(lwork*sizeof(double));
  double * eigvecs   = malloc(np1*np1*sizeof(double));
  double * zeigvecs  = malloc(2*np1*np1*sizeof(double));
  double * DUT       = malloc(2*np1*np1*sizeof(double));

  // Initialize eigvecs
  for(w=0;w<np1; w++){
    for(v=0;v<np1; v++){
      eigvecs[w*np1 + v] = (w==v)?1:0;
    }
  }
  // Specify diagonal and off-diagonal of Diffusion matrix
  for(w=0;w<n;w++){
    offdiag[w] = -sqrt((w+1)*(n-w));
    diag[w]    = n;
    eigvecs[np1*w+w] = 1;
  }
  diag[n]=n;
  eigvecs[np1*n+n] = 1;
  
  // Diagonalize
  dsteqr_(&compz, &np1, diag, offdiag, eigvecs, &np1, workspace, &info);

  // Compute diag*UTranspose, where U is the orthogonal matrix
  // NOTE: FORTRAN MATRICES ARE COLUMN-MAJOR
  // DUT_wv = diag_ww * UT_wv = diag_ww * U_vw
  // Now, flip indices since diag and U are fortran outputs
  for(w=0;w<np1;w++){
    for(v=0;v<np1;v++){
      DUT[2*(np1*v+w)]        = cos(diag[w]*beta)*eigvecs[np1*w + v];
      DUT[2*(np1*v+w)+1]      = -sin(diag[w]*beta)*eigvecs[np1*w + v];
      zeigvecs[2*(np1*v+w)]   = eigvecs[np1*v + w];      
      zeigvecs[2*(np1*v+w)+1] = 0;
    }
  }
  // LAPACK's matrix multiplication needs auxiliary variables 
  double scal1[2], scal2[2];
  scal1[0]=1; scal1[1]=0; scal2[0]=0; scal2[1]=0;
  char tU = 'N';
  char tDUT = 'N';

  // Matrix multiply and output into Diff
  zgemm_(&tU, &tDUT, &np1, &np1, &np1, scal1, zeigvecs, &np1, DUT, &np1, scal2, Diff, &np1);

  free(offdiag);
  free(diag);
  free(workspace);
  free(eigvecs);
  free(zeigvecs);
  free(DUT);
}

// Cost matrix, row major. INITIALIZED TO ALL ZEROS
// Deprecated.
void energize(double gamma, searchInstance *S, double *E){
  int i,j;
  int m = S->m;
  int n = S->n;
  int np1 = n+1;
  double anglej;
  int d;
  // Add transfer terms
  for(i=0; i<m; i++){
    for(j=0; j<m; j++){
      anglej = -(S->cost[j])*gamma;
      d = S->dist[m*i+j];
      E[2*(i*m*np1*np1 + j*np1+d)] = (cos(anglej)-1)/S->combs[d];
      E[2*(m*np1*np1*i + j*np1+d)+1]  = sin(anglej)/S->combs[d];
    }
  }
  // Add identity terms
  for(i=0;i<m*(np1);i++){
    E[2*(m*(np1)*i+i)] += 1;
  }
}

// Take fast matrix powers. ASSUMED COMPLEX in fortran format
// i.e. real value followed by complex value
void mpower(double *M, int T, double *Mpow, int size){
  char tM = 'N';
  double scal1[2], scal2[2];
  scal1[0]=1; scal1[1]=0; scal2[0]=0; scal2[1]=0;
  int i,j, power;

  double *Mcur = malloc(2*size*size*sizeof(double));
  double *Mprod = malloc(2*size*size*sizeof(double));
  double *Mbuf = malloc(2*size*size*sizeof(double));

  // Initialize Mprod to identity (or M if T is odd). This will hold the running product.
  // Initialize Mcur to M. This is the current power M^(2^t)
  // Note: put matrices in column major form for fortran.
  for(i=0;i<size;i++){
    for(j=0;j<size;j++){
      if(T%2){
	Mprod[2*(size*i + j)]   = M[2*(size*j + i)];
	Mprod[2*(size*i + j)+1] = M[2*(size*j + i)+1];
      }
      else   {
	Mprod[2*(size*i + j)]   = (i==j)? 1:0;
	Mprod[2*(size*i + j)+1] = 0;
      }
      Mcur[2*(size*i + j)]      = M[2*(size*j + i)];
      Mcur[2*(size*i + j)+1]    = M[2*(size*j + i)+1];
    }
  }

  // Take matrix powers for the 1s in the binary representation of T
  power = 1;
  while(T>>power){
    zgemm_(&tM, &tM, &size, &size, &size, scal1, Mcur, &size, Mcur, &size, scal2, Mbuf, &size);
    // Copy buffer value into Mcur. This amounts to squaring Mcur.
    for(i=0;i<2*size*size;i++){
      Mcur[i] = Mbuf[i];
    }
    // Now, if T has a digit on the current power of 2, multiply with Mcur
    if((T>>power)&1){
      zgemm_(&tM, &tM, &size, &size, &size, scal1, Mcur, &size, Mprod, &size, scal2, Mbuf, &size);
     // Copy buffer into Mprod
      for(i=0;i<2*size*size;i++){
	Mprod[i] = Mbuf[i];
      }
    }
    // Increment power
    power++;
  }
  
  // Finally, copy Mprod into output, powM, in row-major format
  for(i=0;i<size;i++){
    for(j=0;j<size;j++){
      Mpow[2*(size*i + j)] = Mprod[2*(size*j + i)];
      Mpow[2*(size*i + j)+1] = Mprod[2*(size*j + i)+1];
    }
  }
  free(Mcur);
  free(Mprod);
  free(Mbuf);
}

// Multiple two general matrices (complex)
// Since fortran reads column major, we use essentially pass
// the matrices in reverse, LT * KT = (KL)T
// and then the column-major output can just be read off as KL
void mmmult(double *K, double *L, double *KL, int size){
  int i,j;
  char tM = 'N';
  double scal1[2], scal2[2];
  scal1[0]=1; scal1[1]=0; scal2[0]=0; scal2[1]=0;

  double *KLbuf = malloc(2*size*size*sizeof(double));

  zgemm_(&tM, &tM, &size, &size, &size, scal1, L, &size, K, &size, scal2, KLbuf, &size);

  // Copy buffer into KL
  for(i=0;i<2*size*size;i++){
    KL[i] = KLbuf[i];
  }
  free(KLbuf);
}

// Multiply vector by matrix (complex)
// First, make M column-major by transposing by hand
void mvmult(double *M, double *v, double *Mv, int size){
  int i,j;
  char tM = 'N';
  double scal1[2], scal2[2];
  int inc = 1;
  scal1[0]=1; scal1[1]=0; scal2[0]=0; scal2[1]=0;

  double *Mvbuf = malloc(2*size*sizeof(double));
  double *MT = malloc(2*size*size*sizeof(double));
  // Transpose M to make it column-major
  for(i=0;i<size;i++){
    for(j=0;j<size;j++){
      MT[2*(size*i + j)] = M[2*(size*j + i)];
      MT[2*(size*i + j)+1] = M[2*(size*j + i)+1];
    }
  } 

  zgemv_(&tM, &size, &size, scal1, MT, &size, v, &inc, scal2, Mvbuf, &inc);

  // Copy buffer into Mv
  for(i=0;i<size;i++){
    Mv[2*i] = Mvbuf[2*i];
    Mv[2*i+1] = Mvbuf[2*i+1];
  }
  free(Mvbuf);
  free(MT);
}

// Multiply D, E to find the QAOA1 operator (row-major). 
void qaoa1op(double beta, double gamma, searchInstance *S, double *DE){
  int n = S->n;
  int m = S->m;
  int np1 = n+1;
  double anglej;
  int i,j,w,v;
  int d;

  double *Diff = malloc(2*np1*np1*sizeof(double));

  // Initialize DE to all zeros: time O(m^2*n^2)
  for(i=0;i<m*np1;i++){
    for(j=0;j<m*np1;j++){
      DE[2*(m*np1*i+j)] = 0;
      DE[2*(m*np1*i+j)+1] = 0;
    }
  }
  // It turns out that E is identity plus a very sparse matrix E'
  // So, we write DE = D + D*E'.

  // First, add diffusion blocks to diagonal, D: time O(m*n^2)
  diffuse(beta, S, Diff);
  for(i=0;i<m;i++){
    for(w=0;w<np1;w++){
      for(v=0;v<np1;v++){
	DE[2*((i*np1+w)*m*np1 + (i*np1+v))] = Diff[2*(w*np1 + v)];
	DE[2*((i*np1+w)*m*np1 + (i*np1+v))+1] = Diff[2*(w*np1 + v)+1];
      }
    }
  }
  // Add off-diagonal elements, DE': time O(m^2*n)
  for(i=0;i<m;i++){
    for(j=0;j<m;j++){
      for(w=0;w<np1;w++){
	d = S->dist[m*i+j];
	anglej = -S->cost[j]*gamma;
	DE[2*((i*np1+w)*m*np1 + j*np1 + d)]  += ((cos(anglej)-1)*Diff[2*(0*np1+w)] - sin(anglej)*Diff[2*(0*np1+w)+1])/S->combs[d];
	DE[2*((i*np1+w)*m*np1 + j*np1 + d)+1]+= (sin(anglej)*Diff[2*(0*np1+w)] + (cos(anglej)-1)*Diff[2*(0*np1+w)+1])/S->combs[d];
      }
    }
  }
  free(Diff);
  // Now, the matrix is successfully multiplied and in row-major form.
}

// Perform QAOA1 by preparing D, E and applying to the state
void qaoa1state(double beta, double gamma, searchInstance *S, double *state){
  int i,j, w,v, d, n = S->n, m = S->m;
  int np1 = n+1;
  double anglei;
  double* temp = malloc(2*m*np1*sizeof(double));
  double* Diff = malloc(2*np1*np1*sizeof(double));
  // Action of E matrix. The first two loops give the matrix row index
  // the second two give the vector entry
  for(i=0;i<m;i++){
    for(w=0;w<np1;w++){
      temp[2*(np1*i+w)]   = state[2*(np1*i+w)];
      temp[2*(np1*i+w)+1] = state[2*(np1*i+w)+1];

      // Only in this case, non-trivial multiplication
      if(w==0){
	anglei = -S->cost[i]*gamma;
	for(j=0;j<m;j++){
	  d = S->dist[m*i+j];
	  temp[2*(i*np1+w)]   += ((cos(anglei)-1)*state[2*(np1*j+d)] - sin(anglei)*state[2*(np1*j+d)+1])/(S->combs[d]);
	  temp[2*(i*np1+w)+1] += (sin(anglei)*state[2*(np1*j+d)] + (cos(anglei)-1)*state[2*(np1*j+d)+1])/(S->combs[d]);
	}
      }
    }
  }
  // Action of D matrix
  diffuse(beta, S, Diff);
  for(i=0;i<m;i++){
    for(w=0;w<np1;w++){
      // Reset state
      state[2*(i*np1 + w)]   = 0;
      state[2*(i*np1 + w)+1] = 0;
      // Multiplication loop
      for(v=0;v<np1;v++){
	state[2*(i*np1 + w)]   += Diff[2*(w*np1+v)]*temp[2*(i*np1 + v)] - Diff[2*(w*np1+v)+1]*temp[2*(i*np1+v)+1];
	state[2*(i*np1 + w)+1] += Diff[2*(w*np1+v)]*temp[2*(i*np1 + v)+1] + Diff[2*(w*np1+v)+1]*temp[2*(i*np1+v)];
      }
    }
  }
  free(Diff);
  free(temp);
}


// Perform grover (i.e. same betas, gammas) on a state by applying a matrix power 
void groverT(schedule *U, searchInstance *S, double *state){
  int i,j, n = S->n, m = S->m, L = U->L;
  int np1 = n+1;
  
  // Cyclic protocol, size t iteration repeats T times
  if(U->isCyclic){
    // 1. loop over one cycle, multiple DE's
    double* DEDE = malloc(2*np1*np1*m*m*sizeof(double));
    double* DEcur = malloc(2*np1*np1*m*m*sizeof(double));
    // Initialize both to identity
    for(i=0;i<m*np1; i++){
      for(j=0;j<m*np1; j++){
	DEDE[2*(np1*m*j + i)] = (i==j)?1:0;
	DEDE[2*(np1*m*j + i)+1] = 0;
	DEcur[2*(np1*m*j + i)] = (i==j)?1:0;
	DEcur[2*(np1*m*j + i)+1] = 0;
      }
    }
    if (L==1) {
      qaoa1op(U->beta[0], U->gamma[0], S, DEDE);
    }
    else{
      for(i=0;i<L;i++){
	qaoa1op(U->beta[i], U->gamma[i], S, DEcur);
	mmmult(DEcur, DEDE, DEDE, np1*m); 
      }
    }
    // 2. Exponentiate resulting matrix by T
    mpower(DEDE,U->T, DEDE, np1*m);
    
    // 3. Apply to state
    mvmult(DEDE, state, state, np1*m);
    free(DEDE);
    free(DEcur);
  }
  
  // Not cyclic, run T separate schedule-specified QAOA iterations 
  else           {
    // Check that t is 1
    if (U->t!=1) {
      printf("Acyclic QAOA protocol, running T=%d iterations and assuming t=1.\n", U->T);
    }
    // 1. Apply qaoa1state to state
    for(i=0;i<L;i++){
      //      printf("%lf %lf\n", U->beta[i], U->gamma[i]);
      qaoa1state(U->beta[i], U->gamma[i], S, state);
    }
  }
}

// Norm of the state
double norm(searchInstance *S, double *state, double *mnoms){
  int i,j,w,v,d, k;
  int np1=S->n+1, m=S->m, n=S->n;
  double norm=0;
  for(i=0;i<m;i++){
    for(w=0;w<np1;w++){
      for(j=0;j<m;j++){
	  for(v=0;v<np1;v++){
	    d = S->dist[m*i+j];
	    if ((v-w+d)%2==0&&(v-w+d>=0)&&(w+v-d>=0)){
	      k = (v-w+d)/2;
	      norm+=\
		(state[2*(i*np1+w)]*state[2*(j*np1+v)]+state[2*(i*np1+w)+1]*state[2*(j*np1+v)+1])* \
		mnoms[d*np1+k]*mnoms[(n-d)*np1+v-k]/(S->combs[w]*S->combs[v]);
	    }
	  }
      }
    }
  }
  return norm;
}
// Calculate the success probability of the final state
double successP(searchInstance *S, double *state){
  double amp[2]; amp[0] = 0; amp[1] = 0;
  int d, np1=S->n+1, m=S->m, min=S->min;
  int i, j, lower, upper; // loop variables
  lower = (S->isGrover)?   0:min;
  upper = (S->isGrover)? m-1:min;
  double amp2;
  //double *mnoms = malloc(np1*np1*sizeof(double));
  //multinominate(S->n, mnoms);
  double normalization;
  double success = 0;
  for(i=lower;i<=upper;i++) {
    for(j=0;j<m;j++){
      d = S->dist[m*i+j];
      amp[0] += state[2*(np1*j+d)]/(S->combs[d]); 
      amp[1] += state[2*(np1*j+d)+1]/(S->combs[d]);
    }
    amp2 = amp[0]*amp[0] + amp[1]*amp[1];
    printf("%lf\n",amp2);
    success += amp2;
    amp[0]=0;
    amp[1]=0;
  }
  // TURNS OUT THE STATE STAYS NORMALIZED, SO THIS HAD BETTER BE 1
  //normalization = norm(S, state, mnoms);
  //free(mnoms);
  return success;///normalization;
}
