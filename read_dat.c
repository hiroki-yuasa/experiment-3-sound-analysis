//  音声データにハミング窓をかけるプログラム
//  中身は配布されていたFFTとDFTのプログラムを参考にした

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

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
    if (sdata == NULL || xr == NULL || xi == NULL)
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

    for (i = 0; i < framelen; i++)
    {
        printf("%f %f\n", (float)i / 16000, xr[i]);
    }
}