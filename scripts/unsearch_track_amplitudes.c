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
  reps = (int)(2*sqrt(pow(2,S.n)/(double)(S.m))/U.L);

  // Prepare state
  double * state = malloc(2*S.m*(S.n+1)*sizeof(double));
  prepareState(state, &S);
  for(j=0;j<2*S.m*(S.n+1);j++){
    printf("%lf ", state[j]);
  }
  printf("\n");
  double success;
  for(i=0;i<reps;i++){
    groverT(&U, &S, state);
    for(j=0;j<2*S.m*(S.n+1);j++){
      printf("%lf ", state[j]);
    }
    printf("\n");
  }
  free(state);
}
