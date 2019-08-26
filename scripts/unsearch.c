#include "grover.h"
#define PI 3.1415926535

int main() {
  searchInstance S;
  schedule       U;
  int i,j;
  char inst[500], schd[500];
  scanf("%s %s", inst, schd);
  readInstance(inst, &S);
  readSchedule(schd, &U);
  int np1 = S.n+1;
  printf("n=%d, m=%d, L=%d, T=%d, t=%d, cost=%lf.\n", S.n, S.m, U.L, U.T, U.t, S.cost[0]);
  printf("\n");
  for(i=0;i<S.m;i++){
    for(j=0;j<S.m;j++){
      printf("%d ", S.dist[S.m*i+j]);
    }
    printf("\n");
  }
  // Prepare state
  double * state = malloc(2*S.m*np1*sizeof(double));
  prepareState(state, &S);
  //double success = successP(&S, state);
  //printf("%lf\n", success);
  for(j=0;j<S.m;j++){
    for(i=0;i<np1;i++){
      printf("%5.4f + %5.4f j\n", state[2*(np1*j +i)], state[2*(np1*j+i)+1]);
    }
    printf("\n");
  }
  groverT(&U, &S, state);
  double success = successP(&S, state);
  for(j=0;j<S.m;j++){
    for(i=0;i<np1;i++){
      printf("%5.4f + %5.4f j\n", state[2*(np1*j +i)], state[2*(np1*j+i)+1]);
    }
    printf("\n");
  }
  printf("\n%lf\n", success);
  free(state);
}
