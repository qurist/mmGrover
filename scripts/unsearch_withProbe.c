#include "grover.h"
#define PI 3.1415926535

int main() {
  searchInstance S;
  schedule       U;
  int i,j, reps;
  char inst[500], schd[500];

  // Read input and prepare setup
  scanf("%s %s", inst, schd);
  readInstance(inst, &S);
  readSchedule(schd, &U);
  //Reset schedule parameters to desired number
  reps = (int)(2*sqrt(pow(2,S.n)/(double)(S.m))/U.T);

  // Prepare state
  double * state = malloc(2*S.m*(S.n+1)*sizeof(double));
  prepareState(state, &S);
  double success;
  for(i=0;i<reps;i++){
    groverT(&U, &S, state);
    success = successP(&S, state);
    printf("\n%d %lf\n", (i+1)*U.T, success);
  }
  free(state);
}
