#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "./minIni/dev/minIni.h"

#include <math.h>
#ifndef M_PI
  #define M_PI      3.1415926535897932
#endif

#define sizearray(a)  (sizeof(a) / sizeof((a)[0]))

#define GETSEC(s, section, file) \
        ini_getsection(s, section, sizearray(section), inifile)

typedef struct {
  const char *config;
  const char *output;
  int help;
} ArgParse;

static void printUsage( const char *program_name );
static int argumentsParse( int argc, char *argv[], ArgParse *argparse );
static void getParams( const char *inifile, int *n, int *m, size_t *size );
static void allocator( size_t size, double **e, double **r, double **v );
static void getProblems( const char *inifile, double *e, double *r, double *v );
static inline double dist2BC( const double e, const double E );
static double Kepler( const double e, const double M );
static double rhs( const double z, const double e, const double time );
static void rk4( double *zPos, double *zVel, const double e, const double time, double dt );



int main( int argc, char **argv ) {
  int n, m;
  size_t size;
  double *e, *r, *v;
  double t, h;
  FILE *f;
  ArgParse argparse;

  if ( argumentsParse(argc, argv, &argparse) != 0 ) {
    printUsage( argv[0] );
    return 1;
  }
  if ( argparse.help ) {
    printUsage( argv[0] );
    return 0;
  }
  fprintf( stdout, "using Config: %s\n", argparse.config );
  fprintf( stdout, "writing Output: %s\n", argparse.output );

  getParams( argparse.config, &n, &m, &size );
  allocator( size, &e, &r, &v );
  getProblems( argparse.config, e, r, v );

  h = 2.0 * M_PI / m;
  f = fopen( argparse.output, "w" );
  if ( f ) {
    /* fprintf( f, "%d\t%d\n", n, m ); */
    for ( int k = 0; k < size; ++k ) {
      t = 0.0;
      /* fprintf( f, "%lf\t%lf\t%lf\n", e[k], r[k], v[k] ); */
      for ( int j = 0; j < n; ++j ) {
        for ( int i = 0; i < m; ++i ) {
          rk4( &r[k], &v[k], e[k], t, h );
          t += h;
        }
        fprintf( f, "%lf\t%lf\n", r[k], v[k] );
      }
    }
    fprintf( stdout, "success: written to %s\n", argparse.output );
    fclose( f );
  } else {
    fprintf( stderr, "error: couldn't open file\n" );
    return 1;
  }

  free( e ); free( r ); free( v );
  return 0;
}

/* distance to binary system's barycenter 
 * from one of main bodies */
static inline double 
dist2BC( const double e, const double E )
{
  return 1.0 - e * cos( E );
}
/* solves Kepler's equation via simple iterations */
static double
Kepler( const double e, const double M )
{
  double temp1, temp2, diff;
  int i;

  temp1 = M;
  i = 0;
  for ( ;; ) {
    temp2 = e * sin( temp1 ) + M;
    diff = temp2 - temp1; 
    temp1 = temp2;
    if ( 1e-15 > fabs(diff) )
      break;
    if ( 100 < ++i ) {
      fprintf( stdout, "warning: Eccentric anomaly might be inaccurate\n" );
      break;
    }
  }
  return temp1;
}

static double
rhs( const double z, const double e, const double time )
{
  double n, M, E, rho;
  /* mean motion */
  n = 1.0;
  /* mean anomaly */
  M = n * time;
  M = fmod( M, 2*M_PI );
  /* eccentric anomaly */
  E = Kepler( e, M );
  /* distance to barycenter */
  rho = dist2BC( e, E );
  /* equation of motions */
  return - z / pow( rho * rho + z * z, 3.0 / 2 );
}

/* Runge-Kutta 4th order solver */
static void 
rk4( double *zPos, double *zVel, const double e, const double time, double dt ) 
{
  double k1r, k2r, k3r, k4r, tempr;
  double k1v, k2v, k3v, k4v, tempv;
  double t;

  t = time;
  tempv = *zVel;
  tempr = *zPos; 
  k1r = tempv;
  k1v = rhs( tempr, e, t );

  t = time + dt/2;
  tempv = *zVel + dt/2 * k1v;
  tempr = *zPos + dt/2 * k1r; 
  k2r = tempv;
  k2v = rhs( tempr, e, t );

  t = time + dt/2;
  tempv = *zVel + dt/2 * k2v;
  tempr = *zPos + dt/2 * k2r; 
  k3r = tempv;
  k3v = rhs( tempr, e, t );

  t = time + dt;
  tempv = *zVel + dt * k3v;
  tempr = *zPos + dt * k3r; 
  k4r = tempv;
  k4v = rhs( tempr, e, t );

  *zVel += dt / 6 * ( k1v + 2 * k2v + 2 * k3v + k4v );
  *zPos += dt / 6 * ( k1r + 2 * k2r + 2 * k3r + k4r );
}

static void
getParams( const char *inifile, int *n, int *m, size_t *size )
{
  char str[100];
  char section[50];
  int s = 0, k = 0;
  long x, y;
  for (s = 0; GETSEC(s, section, inifile) > 0; s++) {
    if ( strcmp(section, "integration") == 0 ) {
      int ni = 0, mi = 0;
      for (k = 0; ini_getkey(section, k, str, sizearray(str), inifile) > 0; k++) {
        if ( strcmp(str, "n") == 0 ) ni += 1;
        if ( strcmp(str, "m") == 0 ) mi += 1;
      } 
      if ( k != 2 || ni > 1 || mi > 1 ) {
        fprintf( stderr, "error\n" );
        exit( 1 );
      }
      *n = ini_getl(section, "n", *n, inifile );
      *m = ini_getl(section, "m", *m, inifile );
    }
  } 
  *size = ( size_t ) ( s - 1 );
}

static void
allocator( size_t size, double **e, double **r, double **v )
{
  *e = malloc( size * sizeof(double) );
  *r = malloc( size * sizeof(double) );
  *v = malloc( size * sizeof(double) );
  if ( !e || !r || !v ) {
    fprintf( stderr, "error: memory allocation failed\n" );
    free( e ); free( r ); free( v );
  }
}

static void
getProblems( const char *inifile, double *e, double *v, double *r )
{
  char str[100];
  char section[50];
  int iter = 0;
  int k = 0;
  int s;
  for (s = 0; GETSEC(s, section, inifile) > 0; s++) {
    if ( strncmp(section, "problem", 7) == 0 ) {
      int ei = 0, vi = 0, ri = 0;
      for (k = 0; ini_getkey(section, k, str, sizearray(str), inifile) > 0; k++) {
        if ( strcmp(str, "e") == 0 ) ei += 1;
        if ( strcmp(str, "v") == 0 ) vi += 1;
        if ( strcmp(str, "r") == 0 ) ri += 1;
      } 
      if ( k != 3 || ei > 1 || vi > 1 || ri > 1 ) {
        fprintf( stderr, "error\n" );
        exit( 1 );
      }
      e[iter] = ini_getf(section, "e", e[iter], inifile );
      v[iter] = ini_getf(section, "v", v[iter], inifile );
      r[iter] = ini_getf(section, "r", r[iter], inifile );
      iter++;
    }
  }
}

static void 
printUsage( const char *programName ) {
  fprintf( stderr, "usage: %s [-c config] [-o output] [-h]\n", programName );
  fprintf( stderr, "options:\n" );
  fprintf( stderr, "  -c <file>   Path to configuration file (default: ./config.ini)\n" );
  fprintf( stderr, "  -o <file>   Path to output file (default: ./output.dat)\n" );
  fprintf( stderr, "  -h          Show this help message\n" );
}

static int 
argumentsParse( int argc, char **argv, ArgParse *argparse ) {
  argparse->config = "./config.ini";
  argparse->output = "./output.dat";
  argparse->help = 0;

  for ( int i = 1; i < argc; i++ ) {
    if ( strcmp(argv[i], "-h") == 0 || strcmp(argv[i], "--help") == 0 ) {
      argparse->help = 1;
      return 0;
    }
    else if ( strcmp(argv[i], "-c") == 0 ) {
      if ( i + 1 >= argc ) {
        fprintf( stderr, "Error: -c requires a file path argument.\n" );
        return 1;
      }
      argparse->config = argv[++i]; 
    }
    
    else if ( strcmp(argv[i], "-o" ) == 0) {
      if ( i + 1 >= argc ) {
        fprintf( stderr, "Error: -o requires a file path argument.\n" );
        return 1;
      }
      argparse->output = argv[++i];
    }
    
    else {
      fprintf( stderr, "Error: Unknown argument '%s'\n", argv[i] );
      return 1;
    }
  }

  return 0;
}
