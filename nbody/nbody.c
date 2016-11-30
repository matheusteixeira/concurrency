#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <pthread.h>

/*
 * pRNG based on http://www.cs.wm.edu/~va/software/park/park.html
 */
#define MODULUS    2147483647
#define MULTIPLIER 48271
#define DEFAULT    123456789

static long seed = DEFAULT;
int n_threads;
double dt, dt_old; /* Alterado de static para global */
double max_f;
pthread_mutex_t global_var_mutex;

double Random(void)
/* ----------------------------------------------------------------
 * Random returns a pseudo-random real number uniformly distributed
 * between 0.0 and 1.0.
 * ----------------------------------------------------------------
 */
{
  const long Q = MODULUS / MULTIPLIER;
  const long R = MODULUS % MULTIPLIER;
        long t;

  t = MULTIPLIER * (seed % Q) - R * (seed / Q);
  if (t > 0)
    seed = t;
  else
    seed = t + MODULUS;
  return ((double) seed / MODULUS);
}

/*
 * End of the pRNG algorithm
 */

typedef struct {
    double x, y, z;
    double mass;
} Particle;

typedef struct {
    double xold, yold, zold;
    double fx, fy, fz;
} ParticleV;

struct particles_data{
  Particle* myparticles;
  Particle* others;
  ParticleV* pv;
  int npart, id;
  int start, end;
};

typedef struct { int start, end; } batches;

void InitParticles( Particle[], ParticleV [], int );
double ComputeForces( void *thread_data );
double ComputeNewPos( void *thread_data );


int main(int argc, char **argv)
{
    Particle  * particles;   /* Particles */
    ParticleV * pv;          /* Particle velocity */
    int         npart, i, j;
    int         cnt;         /* number of times in loop */
    double      sim_t;       /* Simulation time */

    pthread_t threads[n_threads];
    pthread_mutex_init(&global_var_mutex, NULL);

    int tmp;
    if(argc != 4){
  		printf("Wrong number of parameters.\nUsage: nbody num_bodies timesteps num_of_threads\n");
    	exit(1);
	  }

  	npart = atoi(argv[1]);
  	cnt = atoi(argv[2]);
    n_threads = atoi(argv[3]);

    dt = 0.001;
  	dt_old = 0.001;

    /* Allocate memory for particles */
    particles = (Particle *) malloc(sizeof(Particle)*npart);
    pv = (ParticleV *) malloc(sizeof(ParticleV)*npart);

   struct particles_data *data;
   data = malloc(sizeof(struct particles_data));

    /* Generate the initial values */
    InitParticles( particles, pv, npart);

    while (cnt--) {
      //compute forces
      for(i = 0; i < n_threads; i++) {
        data->npart = npart;
        data->myparticles = particles;
        data->others = particles;
        data->pv = pv;
        data->start = i*(npart/n_threads);
        data->end = i*(npart/n_threads) + (npart/n_threads);
        pthread_create(&threads[i], NULL, ComputeForces, (void*)data);
      }

      for(i = 0; i < n_threads; i++){
        pthread_join(threads[i], NULL);
      }
      /* Once we have the forces, we compute the changes in position */
      for(i = 0; i < n_threads; i++) {
        data->npart = npart;
        data->myparticles = particles;
        data->others = particles;
        data->pv = pv;
        data->start = i*(npart/n_threads);
        data->end = i*(npart/n_threads) + (npart/n_threads);
        pthread_create(&threads[i], NULL, ComputeNewPos, (void*)data);
      }

      for(i = 0; i < n_threads; i++){
        pthread_join(threads[i], NULL);
      }
    }
    for (i=0; i<npart; i++)
      fprintf(stdout,"%.5lf %.5lf\n", particles[i].x, particles[i].y);

    pthread_mutex_destroy(&global_var_mutex);
    pthread_exit(NULL);
    free(data);
    return 0;
}

void InitParticles(Particle particles[], ParticleV pv[], int npart)
{
    int i;
    for (i=0; i<npart; i++) {
      particles[i].x	  = Random();
      particles[i].y	  = Random();
      particles[i].z	  = Random();
      particles[i].mass = 1.0;
      pv[i].xold	  = particles[i].x;
      pv[i].yold	  = particles[i].y;
      pv[i].zold	  = particles[i].z;
      pv[i].fx	  = 0;
      pv[i].fy	  = 0;
      pv[i].fz	  = 0;
    }
}

double ComputeForces( void *thread_data )
{
  struct particles_data *received_data = thread_data;

  int i;
  max_f = 0.0;

  int npart = received_data->npart;
  Particle* myparticles = received_data->myparticles;
  Particle* others = received_data->others;
  ParticleV* pv = received_data->pv;

  int start = received_data->start;
  int end = received_data->end;

  for (i = start; i <  end; i++) {
    int j;
    double xi, yi, mi, rx, ry, mj, r, fx, fy, rmin;
    rmin = 100.0;
    xi   = myparticles[i].x;
    yi   = myparticles[i].y;
    fx   = 0.0;
    fy   = 0.0;
    for (j= 0; j < npart; j++) {
      rx = xi - others[j].x;
      ry = yi - others[j].y;
      mj = others[j].mass;
      r  = rx * rx + ry * ry;
      /* ignore overlap and same particle */
      if (r == 0.0) continue;
      if (r < rmin) rmin = r;
      r  = r * sqrt(r);
      fx -= mj * rx / r;
      fy -= mj * ry / r;
    }
    pv[i].fx += fx;
    pv[i].fy += fy;
    fx = sqrt(fx*fx + fy*fy)/rmin;
    pthread_mutex_lock (&global_var_mutex);
    if (fx > max_f) max_f = fx;
    pthread_mutex_unlock (&global_var_mutex);
    pthread_exit((void*) 0);
  }
  return max_f;
}

double ComputeNewPos( void *thread_data )
{
  struct particles_data *received_data = thread_data;

  Particle* particles = received_data->myparticles;
  ParticleV* pv = received_data->pv;

  int start = received_data->start;
  int end = received_data->end;

  int i;
  double a0, a1, a2;
  double dt_new;

  a0	 = 2.0 / (dt * (dt + dt_old));
  a2	 = 2.0 / (dt_old * (dt + dt_old));
  a1	 = -(a0 + a2);

  for (i=start; i < end; i++) {
    double xi, yi;
    xi	           = particles[i].x;
    yi	           = particles[i].y;
    particles[i].x = (pv[i].fx - a1 * xi - a2 * pv[i].xold) / a0;
    particles[i].y = (pv[i].fy - a1 * yi - a2 * pv[i].yold) / a0;
    pv[i].xold     = xi;
    pv[i].yold     = yi;
    pv[i].fx       = 0;
    pv[i].fy       = 0;
  }
  dt_new = 1.0/sqrt(max_f);
  /* Set a minimum: */
  if (dt_new < 1.0e-6) dt_new = 1.0e-6;
  /* Modify time step */
  pthread_mutex_lock (&global_var_mutex);
  if (dt_new < dt) {
    dt_old = dt;
    dt     = dt_new;
  }
  else if (dt_new > 4.0 * dt) {
      dt_old = dt;
    dt    *= 2.0;
  }
  pthread_mutex_unlock (&global_var_mutex);
  return dt_old;
}
