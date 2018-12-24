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
    int i, tau, j;
    double maxtau = 0, peak = 0;

    /* check the validity of input */
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
    double *r = (double *)calloc(framelen, sizeof(double));
    double *v = (double *)calloc(framelen, sizeof(double));
    double *f = (double *)calloc(53, sizeof(double));
    if (sdata == NULL || xr == NULL || r == NULL || v == NULL || f == NULL)
        exit(1);

    while (nskip < 53 * 256)
    {
        if ((fpDAT = fopen(argv[1], "r")) == NULL)
            exit(1);
        fseek(fpDAT, nskip * sizeof(short), SEEK_SET);
        fread(sdata, sizeof(short), framelen, fpDAT);
        fclose(fpDAT);

        /* windowing */
        for (i = 0; i < framelen; i++)
        {
            xr[i] = (0.54 - 0.46 * cos(2 * M_PI * i / (framelen - 1))) * sdata[i];
        }

        for (tau = 0; tau < framelen; tau++)
        {
            for (i = 0; i < (framelen - tau); i++)
            {
                r[tau] += xr[i] * xr[i + tau];
            }
            v[tau] = r[tau] / r[0];
        }

        for (i = 0; i < framelen; i++)
        {
            // if (nskip == 2 * framelen / 2 || nskip == 30 * framelen / 2)
            // printf("%f %f\n", (float)i / 16000, v[i]);
            if (v[i] > peak && i > 30)
            {
                peak = v[i];
                maxtau = (double)i / 16000;
            }
        }

        f[nskip / 256] = 1 / maxtau;

        maxtau = 0;
        peak = 0;
        nskip += framelen / 2;
        fprintf(stderr, "trial %d\n", nskip / 256);
        free(r);
        r = (double *)calloc(framelen, sizeof(double));
    }

    for (i = 0; i < 53; i++)
        printf("%d %f\n", i, f[i]);
}