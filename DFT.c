#include <stdio.h>
#include <stdlib.h>
#include <math.h>

/* Functions prototypes for DFT & IDFT */

void DFT(double *xr, double *xi, double *Xr, double *Xi, int N);
void IDFT(double *Xr, double *Xi, double *xr, double *xi, int N);

/* implementation of DFT */
void DFT(double *xr, double *xi, double *Xr, double *Xi, int N)
{
  int k, n;
  for (k = 0; k < N; k++)
  {
    Xr[k] = Xi[k] = 0.0; // initialization
    for (n = 0; n < N; n++)
    {
      double wr = cos(2.0 * M_PI * n * k / N);
      double wi = sin(2.0 * M_PI * n * k / N);
      Xr[k] += xr[n] * wr + xi[n] * wi;  /* should be modified. (hint using wr & wi) */
      Xi[k] += -xr[n] * wi + xi[n] * wr; /* should be modified. (hint using wr & wi) */
    }
  }
  return;
}

/* implementation of IDFT */
void IDFT(double *Xr, double *Xi, double *xr, double *xi, int N)
{
  int k, n;
  for (n = 0; n < N; n++)
  {
    xr[n] = xi[n] = 0.0; // initialization
    for (k = 0; k < N; k++)
    {
      double wr = cos(2.0 * M_PI * n * k / N);
      double wi = sin(2.0 * M_PI * n * k / N);
      xr[n] += Xr[k] * wr - Xi[k] * wi; /* should be modified. (hint using wr & wi) */
      xi[n] += Xr[k] * wi + Xi[k] * wr; /* should be modified. (hint using wr & wi) */
    }
    xr[n] /= N;
    xi[n] /= N;
  }
  return;
}

int main(int argc, char **argv)
{
  /* check the format of input */
  if (argc != 4)
  {
    fprintf(stderr, "Usage: %s DATfile skip[sample] frame_length[sample]\n", argv[0]);
    exit(1);
  }
  FILE *fpDAT;
  int nskip;
  int framelen;
  int i;

  /* check the validity of input */
  if ((fpDAT = fopen(argv[1], "r")) == NULL)
    exit(1);
  if ((nskip = atoi(argv[2])) < 0)
    exit(1);
  if ((framelen = atoi(argv[3])) < 0)
    exit(1);

  fprintf(stderr, "# DATfile = %s\n", argv[1]);
  fprintf(stderr, "# %d samples are skipped.\n", nskip);
  fprintf(stderr, "# 1 frame contains %d sampels.\n", framelen);

  /* memory allocation & initilization */
  /* calloc() puts zero-values for assigned memories. */
  short *sdata = (short *)calloc(framelen, sizeof(short));
  double *xr = (double *)calloc(framelen, sizeof(double));
  double *xi = (double *)calloc(framelen, sizeof(double));
  double *Xr = (double *)calloc(framelen, sizeof(double));
  double *Xi = (double *)calloc(framelen, sizeof(double));
  if (sdata == NULL || xr == NULL || xi == NULL || Xr == NULL || Xi == NULL)
    exit(1);

  fseek(fpDAT, nskip * sizeof(short), SEEK_SET);
  fread(sdata, sizeof(short), framelen, fpDAT);
  fclose(fpDAT);

  /* windowing */
  for (i = 0; i < framelen; i++)
  {
    xr[i] = (0.54 - 0.46 * cos(2 * M_PI * i / (framelen - 1))) * sdata[i];
    xi[i] = 0.0;
  }

/* DFT */
#ifdef STRESSTEST
  for (i = 0; i < 1000; i++)
#endif
    DFT(xr, xi, Xr, Xi, framelen);

/* plot the result */
#ifndef INVERSE
  for (i = 0; i < framelen; i++)
  {
    double spec = log10((Xr[i] * Xr[i] + Xi[i] * Xi[i]) / framelen);
    printf("%f %lf\n", (float)i * 16000 / framelen, spec);
  }

#else
  /* IDFT trial */
  IDFT(Xr, Xi, xr, xi, framelen);

  for (i = 0; i < framelen; i++)
    printf("%d %lf\n", i, xr[i]);
#endif

  return 0;
}
