#include <stdio.h>
#include <math.h>
#include <time.h>
#include <stdlib.h>
#include <stddef.h>
#include <ctype.h>
//#include "make_tm.h"
//#include "sphere.h"
#include "randn.h"

/* lwseed.c -- Create initial distribution of particles for leeway model
 *
 * This is a seed preprocessor to the leeway model. It reads a simple
 * ASCII input and yields an input file for the leeway model.
 *
 * Eight degrees of freedom are used to define an initial distribution of
 * particles in time and space; two positions (lon,lat), two radii of
 * uncertainty, and a period of time.
 *
 * The particles are spread out in a cone-shaped region starting at position
 * (lon0,lat0) with radius r0 and ending in position
 * (lon1,lat1) with radius r1. The exact positions of the particles are
 * determined by a random normal deviation from a center position along
 * a great circle arc connecting the start and end positions.
 *
 * The radius of uncertainty equals two standard deviations (2*sigma) in a
 * circular normal distribution, i.e., 86% of your particles will on average
 * fall within the radius. Because we want the particle cloud to assume a
 * cone-shape, we let the radius vary linearly from r0 to r1.
 *
 * To compile:
 *
 *  % make lwseed
 *
 * Usage:
 *
 *  % lwseed lwseed.in leeway.in
 *
 * Requires:
 *
 *  - lwseed.in - the input file containing the definition of the initial search
 *     area and the number of the search object
 *  - OBJECTPROP.DAT containing the list of leeway coefficients for
 *     different search objects.
 *
 * oyvind.breivik@met.no
 *
 */


/* Constants */
#define USAGE "Usage: %s lwseed.in leeway.in\n"
#define VER 2.5     /* model version */
#define RW 0        /* right-of-wind placeholder */
#define LW 1        /* left-of-wind placeholder */
#define DW 2        /* downwind placeholder */
#define MAXLEN 200  /* string length */

#define TOL 1.0     /* [m] */
#define C 0.5
#define NMEM 160000    /* default ensemble size */
#define OBJECTFILE "OBJECTPROP.DAT" /* object property data file */
#define COORD "C:/cygwin/home/shchekin/Documents/ChapterMed/HighWindage/Winter2014/coordregion.txt"
// #define COORD "/cygdrive/e/ChapterMed/Winter/coordregion.txt"
#define PJIBE 0.04  /* Hourly probability of jibing */

#define UCSD 0.0 /* East current std dev */
#define VCSD 0.0 /* North current std dev */
#define UWSD 0.6 /* East wind std dev */
#define VWSD 0.6 /* North wind std dev */
#define UWTS 0.0 /* Integral time scale of east wind (not in use yet) */
#define VWTS 0.0 /* Integral time scale of north wind (not in use yet) */
#define UCTS 0.0 /* Integral time scale of east current (not in use yet) */
#define VCTS 0.0 /* Integral time scale of east current (not in use yet) */



 int main(int argc, char *argv[]) {

 extern float decdeg(int deg, float dmin);
 printf("ku-ku");
   float *lon, *lat; /* particle position [deg] */
   float *rdw, *rcw; /* random perturbations downwind & crosswind */
   int *ort;         /* orientation of particle */
   time_t *birth;    /* particle release time [s from start t0] */


 /* Leeway linear regression parameters */

   float lwsd[3], lwa[3], lwb[3]; /* leeway linear regression parameters */


 /* Start and end of seeding */

   int yyyy0, mm0, dd0, hh0, mi0, yyyy1, mm1, dd1, hh1, mi1;


 /* Time parameters */

   time_t simtime; /* System time [UTC] for seeding the random number generator */
   //time_t t0, t1;
   int t0, t1;

 /* Strings */

   char objectname[MAXLEN];
   char line[MAXLEN];
   char datestr[MAXLEN];
   char idtag[MAXLEN];


 /* Various */

   int nmem = NMEM; /* default ensemble size */
   float ver;         /* seeder version */
   float lon0, lat0, lon1, lat1, dmin;
   float arc, darc, r0, r, dr, r1, dir0, dir, x1, x2, clon, clat, a;
   float eps, rmem;
   int k, l, objtypeno, deg;

   int constwind, constcurr; /* logical: true if constant wind or current */
   float uc0, vc0;     /* constant current vector */
   float uw0, vw0;     /* ditto for wind */
   int nostrand;       /* do particles strand? false=0 */
   int outputinterval; /* output timestep [s] */
   int det_cluster; /*maximum distance between points in cluster */
   int NMINPOINTS; /*minimum number of points in cluster */

 /* Wind and current error statistics */

   float ucsd = UCSD, vcsd = VCSD; /* std dev of east & north current comp */
   float uwsd = UWSD, vwsd = VWSD; /* ditto for wind */
   float ucts = UCTS, vcts = VCTS; /* timescale of curr fluctuations (not in use yet) */
   float uwts = UWTS, vwts = VWTS; /* ditto for wind (not in use yet) */

   printf("ku-ku");

 /* File pointers */

   FILE *fin, *fout, *fobject, *fcoord;
 /*  char *windfile = "ECMWFf121014-23.grb";*/
  /* char *windfile = "/home/Elena/Documents/Leeway_gshhs/case_study_wind/20090621-ECMWF-wind-20090623.grb";
   char *currfile = "/home/Elena/Documents/Leeway_gshhs/case_study/20090621-MFS-surf-20090623.grb";*/
 /*  char *windfile = "/home/Elena/Documents/Leeway_gshhs/case_study_ECMWF/case_study/20081030-ECMWF-20081120.grb";*/
 /*  char *windfile = "/home/Elena/Documents/Leeway_gshhs/case_study_ECMWF/case_study/20081030_20081120_ECMWF.grb"; */
 /*  char *windfile = "/home/Elena/Documents/Leeway_gshhs/case_study_ECMWF/case_study/20131004_20131008_ECMWF.grb"; */
  /* char *windfile = "/home/Elena/Documents/Leeway_gshhs/case_study_ECMWF/case_study/20081030_20081120_ECMWF.grb"; */
   char *windfile = "C:/cygwin/home/shchekin/Documents/ChapterMed/HighWindage/Winter2014/mean_Ww2014.txt";
//   char *windfile = "/cygdrive/e/ChapterMed/mean_Ww2013_2016.txt";
 //   char *windfile = "/cygdrive/e/ChapterMed/Winter/mean_Ww2013_2016.txt";
/*   char *currfile = "/home/Elena/Documents/Leeway_gshhs/CURR/20081030_UV_20081120.grb_5m";  */
/*   char *currfile = "/home/Elena/Documents/Leeway_gshhs/case_study_AFS/case_study/20081030_UV_20081120.grb";  */
 /*  char *currfile = "/home/Elena/Documents/Leeway_gshhs/case_study_AFS/case_study/20131004_UV_20131008_30m.grb";*/
  /* char *currfile = "/home/Elena/Documents/Leeway_gshhs/case_study_AFS/case_study/20081030_UV_20081120.grb";  */
   char *currfile = "C:/cygwin/home/shchekin/Documents/ChapterMed/HighWindage/Winter2014/mean_Wc2014_2";
  //  char *currfile = "/cygdrive/e/ChapterMed/mean_Wc2013_2016_2"; // winter current
  //  char *currfile = "/cygdrive/e/ChapterMed/Winter/mean_Wc2013_2016.txt"; //summer current
//  char *currfile = "/home/Elena/Documents/Leeway_gshhs/case_study_AFS/case_study/20081030_UV_20081104.grb";

 /*  char *currfile = ".grb"; */
   char *coastfile = "C:/cygwin/home/shchekin/Documents/ChapterMed/gshhs/gshhs_f.bs";
  // char *coastfile = "/cygdrive/e/ChapterMed/gshhs/gshhs_f.bs";

  //  char *coordfile = "/cygdrive/e/ChapterMed/Winter/coordregion.txt";
    char *coordfile = "C:/cygwin/home/shchekin/Documents/ChapterMed/HighWindage/Winter2014/coordregion.txt";

 /* Size up arrays [nmem] */

   lon =   (float *) malloc(nmem * sizeof(float));
   lat =   (float *) malloc(nmem * sizeof(float));
   ort =   (int *)   malloc(nmem * sizeof(int));
   rdw =   (float *) malloc(nmem * sizeof(float));
   rcw =   (float *) malloc(nmem * sizeof(float));
   birth = (time_t *) malloc(nmem * sizeof(float));


 /* Calculate seeding parameters */

 /* Grab system date & time [UTC] */

   simtime = time(NULL);
   strftime(datestr, MAXLEN, "%Y-%m-%d\t%H:%M:%S", gmtime(&simtime));


/* Seed random number generator with current calendar time (epoch seconds) */

   srand((unsigned) simtime);


 /* Help line */

   if (argc == 1) {
      fprintf(stderr, "Version %f\n", VER);
      fprintf(stderr, USAGE, argv[0]);
      exit(8);
   }


/* Open input and output files */

   if (argc < 4) {
      fprintf(stderr, USAGE, argv[0]);
      exit(1);
   }
   if ((fin = fopen(argv[1],"r")) == NULL) {
      fprintf(stderr,"Could not open file: %s\n", argv[1]);
      exit(1);
   }
   if ((fout = fopen(argv[2],"w")) == NULL) {
      fprintf(stderr,"Could not open file: %s\n", argv[2]);
      exit(1);
   }
   if ((fobject = fopen(OBJECTFILE, "r")) == NULL) {
      fprintf(stderr,"Could not open file: %s\n", OBJECTFILE);
      exit(1);
   }
   if ((fcoord = fopen(argv[3], "r")) == NULL) {
      fprintf(stderr,"Could not open coordinate file: %s\n", argv[3]);
      exit(1);
   }

 /* Read input file fin */

 fgets(line, MAXLEN, fin);    /*  line 1 - ignore */

   fgets(line, MAXLEN, fin);         /*  2 - version */
   sscanf(line, "%g", &ver);

   fgets(line, MAXLEN, fin);         /*  3 - lon0 (startlon) [deg] */
   sscanf(line, "%d", &deg);

   fgets(line, MAXLEN, fin);
   sscanf(line, "%g", &dmin);        /*  4 - decimal minutes */
   lon0 = decdeg(deg, dmin);         /*  convert to decimal [deg] */

   fgets(line, MAXLEN, fin);
   sscanf(line, "%d", &deg);         /*  5 - lat0 (startlat) [deg] */
   fgets(line, MAXLEN, fin);
   sscanf(line, "%g", &dmin);        /*  6 - decimal minutes */
   lat0 = decdeg(deg, dmin);

   fgets(line, MAXLEN, fin);
   sscanf(line, "%d", &deg);         /*  7 - lon1 (endlon) [deg] */
   fgets(line, MAXLEN, fin);
   sscanf(line, "%g", &dmin);        /*  8 - decimal minutes */
   lon1 = decdeg(deg, dmin);

   fgets(line, MAXLEN, fin);
   sscanf(line, "%d", &deg);         /*  9 - lat1 (endlat) [deg] */
   fgets(line, MAXLEN, fin);
   sscanf(line, "%g", &dmin);        /* 10 - decimal minutes */
   lat1 = decdeg(deg, dmin);

   fgets(line, MAXLEN, fin);
   sscanf(line, "%g", &r0);          /* 11 - r0, start radius [km] */
   fgets(line, MAXLEN, fin);
   sscanf(line, "%g", &r1);          /* 12 - r1, end radius [km] */

   fgets(line, MAXLEN, fin);
   sscanf(line, "%s", idtag);        /* 13 - ID tag */

   fgets(line, MAXLEN, fin);
   sscanf(line, "%d", &objtypeno);   /* 14 - object type no */

   fgets(line, MAXLEN, fin);
   sscanf(line, "%d", &yyyy0);       /* 15 - start date - year */

   fgets(line, MAXLEN, fin);
   sscanf(line, "%d", &mm0);         /* 16 - start date - month */

   fgets(line, MAXLEN, fin);
   sscanf(line, "%d", &dd0);         /* 17 - start date - day */

   fgets(line, MAXLEN, fin);
   if (isdigit(line[3]) || isdigit(line[4])) {
      sscanf(line, "%d%d", &hh0, &mi0); /* 18 - start date - hour, minutes */
   }
   else {
      sscanf(line, "%d", &hh0);         /* 18 - start date - hour only */
      mi0=0;
   }

  time_t rawtime;
  struct tm  info;

  time ( &rawtime );
  //info = localtime ( &rawtime );
  info.tm_year = yyyy0-1900;
  info.tm_mon = mm0-1;
  info.tm_mday = dd0;
  info.tm_hour = hh0;
  info.tm_min = mi0;
  info.tm_sec = 1;
  info.tm_isdst = -1;

  //t0 = mktime(&info); /* start time */
  //t0=mktime(&info);
  t0=0;
  //printf("%d\t%s\n",t0,asctime(&info));


   fgets(line, MAXLEN, fin);
   sscanf(line, "%d", &yyyy1);      /* 19 - end date - year */

   fgets(line, MAXLEN, fin);
   sscanf(line, "%d", &mm1);        /* 20 - end date - month */

   fgets(line, MAXLEN, fin);
   sscanf(line, "%d", &dd1);        /* 21 - end date - day */

   fgets(line, MAXLEN, fin);
   if (isdigit(line[3]) || isdigit(line[4])) {
      sscanf(line, "%d%d", &hh1, &mi1); /* 22 - end date - hour, minutes */
   }
   else {
      sscanf(line, "%d", &hh1);         /* 22 - end date - hour only */
      mi1=0;
   }

  time ( &rawtime );
  //info = localtime ( &rawtime );
  info.tm_year = yyyy1-1900;
  info.tm_mon = mm1-1;
  info.tm_mday = dd1;
  info.tm_hour = hh1;
  info.tm_min = mi1;
  info.tm_sec = 1;
   info.tm_isdst = -1;

 //  t1 = mktime(&info); /* end of seeding */
   t1= 24;
   printf("%d\n",t1);
   fgets(line, MAXLEN, fin);
   sscanf(line, "%d%g%g", &constcurr, &uc0, &vc0); /* 23 - const curr */

   fgets(line, MAXLEN, fin);
   sscanf(line, "%d%g%g", &constwind, &uw0, &vw0); /* 24 - const wind */

   fgets(line, MAXLEN, fin);
   sscanf(line, "%d", &nostrand);                  /* 25 - particles do not strand (false=0) */

   fgets(line, MAXLEN, fin);
   sscanf(line, "%d", &outputinterval);            /* 26 - output timestep */

   fgets(line, MAXLEN, fin);
   sscanf(line, "%d", &NMINPOINTS);            /* 27 - minimum number of drifters in cluster */

 /*   * Read leeway drift properties from file fobject.
 *     * Note that this file may be edited and extended at will.
 *       */

 /* Browse to correct object class */
   for (l=0; l<objtypeno; l++) {

      /* Read object name */
      fgets(line,MAXLEN,fobject);       /* line 1 - object name */
      sscanf(line, "%s", objectname);
      fgets(line,MAXLEN,fobject);       /*      2 - ignore */

      /* Read object drift properties.
 *        *  Downwind component:    slope, offset, std dev
 *               *  Right of downwind (+): slope, offset, std dev
 *                      *  Left  of downwind (-): slope, offset, std dev */
      fgets(line,MAXLEN,fobject);       /*      3 - object properties */
      sscanf(line, "%g%g%g%g%g%g%g%g%g", &lwa[DW],&lwb[DW],&lwsd[DW],&lwa[RW],&lwb[RW],&lwsd[RW], &lwa[LW],&lwb[LW],&lwsd[LW]);

   } /* end for l */

   int i, j;
   for (i=0; i<nmem; i++) {
      /* Read object name */
      fscanf(fcoord, "%g\t%g\n", &lon[i],&lat[i]);
      printf("%g\t%g\n", lon[i],lat[i]);

   } /* end for l */

   rmem=(float) nmem;

 /* Seed particles */

   for (k=0; k<nmem; k++)
   {

      /* Set particle orientation */
      ort[k] = (k+1)%2; /* odd to the left, even to the right of downwind */
  /*    ort[k] = 1;  odd to the left, even to the right of downwind */

      /* Generate normal, N(0,1), random perturbations for leeway coeffs.
 *        * Negative downwind slope coefficients must be avoided as particles
 *               * should drift downwind. The problem arises because of high error
 *                      * variances (see e.g. PIW-1).
 *                             */
      do {
         rdw[k] = randn();      /* downwind */
         eps = rdw[k]*lwsd[DW];
      } while (lwa[DW]+eps/20.0<0.0);
      rcw[k] = randn(); /* crosswind */

      /* Compute time of birth of particles [seconds from t0] */
      birth[k] = (time_t) (((float)(t1-t0)*k*3600)/nmem); /* [s] */

   } /* end for k - seed particles */
/* Write seed file, leeway.in */

 fprintf(fout, "# Drift simulation initiated [UTC], output interval [s] & identifier:\n");
   fprintf(fout, "simDate\tsimTime\toutputInterval\tID\n");
   fprintf(fout, "%s\t%7d\t%s\n", datestr, outputinterval, idtag);

   fprintf(fout, "# Model version, no strand option (false=0), coastline:\n");
   fprintf(fout, "modelVersion\tnoStrand\tcoastFile\n");
   fprintf(fout, "%5.2f\t%1d\t%s\n", VER, 0, coastfile);

   fprintf(fout, "# Object class id & name:\n");
   fprintf(fout, "objectClassId\tobjectClassName\n");
   fprintf(fout, "%3d\t%29s\n", objtypeno, objectname);

   fprintf(fout, "# Seeding start time, position & radius:\n");
   fprintf(fout, "startDate\tstartTime\tstartLon\tstartLat\tstartRad\n");
   fprintf(fout, "%4d-%02d-%02d\t%02d:%02d:%02d\t%9.4f\t%8.4f\t%9.3f\n",\
      yyyy0, mm0, dd0, hh0, mi0, 0, lon0, lat0, r0);

   fprintf(fout, "# Seeding end time, position & radius:\n");
   fprintf(fout, "endDate\tendTime\tendLon\tendLat\tendRad\n");
   fprintf(fout, "%4d-%02d-%02d\t%02d:%02d:%02d\t%9.4f\t%8.4f\t%9.3f\n",\
      yyyy1, mm1, dd1, hh1, mi1, 0, lon1, lat1, r1);

   fprintf(fout, "# Total no of seeded particles:\n");
   fprintf(fout, "seedTotal\n");
   fprintf(fout, "%4d\n", nmem);

   fprintf(fout, "# Right leeway coeffs; slope [%%], offset [cm/s], std dev [cm/s]:\n");
   fprintf(fout, "aRight\tbRight\tsdRight\n");
   fprintf(fout, "%7.2f\t%7.2f\t%7.2f\n", lwa[RW], lwb[RW], lwsd[RW]);

   fprintf(fout, "# Left leeway coeffs; slope [%%], offset [cm/s], std dev [cm/s]:\n");
   fprintf(fout, "aLeft\tbLeft\tsdLeft\n");
   fprintf(fout, "%7.2f\t%7.2f\t%7.2f\n", lwa[LW], lwb[LW], lwsd[LW]);

   fprintf(fout, "# Downwind leeway coeffs; slope [%%], offset [cm/s], std dev [cm/s]:\n");
   fprintf(fout, "aDownwind\tbDownwind\tsdDownwind\n");
   fprintf(fout, "%7.2f\t%7.2f\t%7.2f\n", lwa[DW], lwb[DW], lwsd[DW]);

   fprintf(fout, "# Hourly probability of jibing [%%]:\n");
   fprintf(fout, "pJibe\n");
   fprintf(fout, "%5.2f\n", PJIBE*100.0);

   fprintf(fout, "# Wind stats; east & north std dev [m/s]; east & north integral time scale [h]:\n");
   fprintf(fout, "sdEastwind\tsdNorthwind\ttsEastwind\ttsNorthwind\n");
   fprintf(fout, "%7.2f\t%7.2f\t%7.2f\t%7.2f\n", uwsd, vwsd, uwts, vwts);

   fprintf(fout, "# Current stats; east & north std dev [m/s]; east & north integral time scale [h]:\n");
   fprintf(fout, "sdEastcurrent\tsdNorthcurrent\ttsEastcurrent\ttsNorthcurrent\n");
   fprintf(fout, "%7.2f\t%7.2f\t%7.2f\t%7.2f\n", ucsd, vcsd, ucts, vcts);

   fprintf(fout, "# Wind file [Lat-lon GRIB]:\n");
   fprintf(fout, "windFile\n");
   if (constwind==1) {
      fprintf(fout, "constant %7.2f %7.2f\n", uw0, vw0);
   } else {
      fprintf(fout, "%s\n", windfile);
   }

   fprintf(fout, "# Current file [Lat-lon GRIB]:\n");
   fprintf(fout, "currentFile\n");
   if (constcurr==1) {
      fprintf(fout, "constant %7.3f %7.3f\n", uc0, vc0);
   } else {
      fprintf(fout, "%s\n", currfile);
   }

   fprintf(fout, "#Minimum number of drifters in cluster :\n");
   fprintf(fout, "NMINPOINTS\n");
   fprintf(fout, " %d\n",NMINPOINTS);


   /* Write individual particles */
   fprintf(fout, "\n"); /* blank line before particle data */
   fprintf(fout, "# Particle data:\n");
   fprintf(fout, "id\tlon\tlat\torientation\trdw\trcw\tbirth\n");
   for (k=0; k<nmem; k++)
   {
      fprintf(fout, "%4d\t%g\t%g\t%2d\t%g\t%g\t%d\n", k+1,\
            lon[k], lat[k], ort[k], rdw[k], rcw[k], birth[k]);
   } /* end for k */

   fclose(fin);
   fclose(fobject);
   fclose(fout);
   fclose(fcoord);

  return 0;

}
float decdeg(int deg, float dmin)
{
   int sign;

   sign = (deg >= 0 && dmin >= 0.0? 1 : -1);
   return sign*(fabs((float) deg) + fabs(dmin/60.0));
}


