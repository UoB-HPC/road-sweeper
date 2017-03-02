
#ifndef WORK
#define WORK 50
#endif

void compute(void) {
  volatile double x = 0.0;
  for (int i = 0; i < WORK; i++) {
    x += 1.0;
  }
}

