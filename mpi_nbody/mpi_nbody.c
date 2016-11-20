#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <mpi.h>

/*
 * pRNG based on http://www.cs.wm.edu/~va/software/park/park.html
 */
#define MODULUS    2147483647
#define MULTIPLIER 48271
#define DEFAULT    123456789

static long seed = DEFAULT;
double dt, dt_old; /* Alterado de static para global */
int *buffer, *global_buffer;
int batch_size, batch_start, batch_end;
int rank, size;
double global_max_f;

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

void InitParticles( Particle[], ParticleV [], int );
double ComputeForces( Particle [], Particle [], ParticleV [], int , int, int);
double ComputeNewPos( Particle [], ParticleV [], int, double);

int main(int argc, char **argv)
{
    Particle  * particles;   /* Particles */
    ParticleV * pv;          /* Particle velocity */
    int         npart, i, j;
    int         cnt;         /* number of times in loop */
    double      sim_t;       /* Simulation time */
    int tmp;

    /*function void *calloc(size_t nitems, size_t size) allocates the requested
    memory and returns a pointer to it. The difference in malloc and calloc is
    that malloc does not set the memory to zero where as calloc sets allocated
    memory to zero.
    https://www.tutorialspoint.com/c_standard_library/c_function_calloc.htm
    */
    buffer = calloc(npart, sizeof(int));
    global_buffer = calloc(npart, sizeof(int));

    if(argc != 3){
		    printf("Wrong number of parameters.\nUsage: nbody num_bodies timesteps\n");
		    exit(1);
	  }

    npart = atoi(argv[1]);
  	cnt = atoi(argv[2]);
  	dt = 0.001;
  	dt_old = 0.001;

    /* Allocate memory for particles */
    particles = (Particle *) malloc(sizeof(Particle)*npart);
    pv = (ParticleV *) malloc(sizeof(ParticleV)*npart);

    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    batch_size = (npart/size);
    int ns_remainder = npart % size;

      if (rank == 0) {
        for (int i = 1; i < size; i++) {
            /*
            MPI_RECv Params:
            buf: ponteiro para o dado a ser recebido
            count: num de elementos a serem recebidos
            dtype: tipo dos elementos
            src: rank do remetente
            tag: um num para classificar msgs (de quem vai receber)
            comm: comunicador
            status: infos sobre a mensagem recebida
            */
            MPI_Recv(global_buffer, npart, MPI_INT, i, 0, MPI_COMM_WORLD, NULL);
            /*
            Note that MPI_STATUS_IGNORE is not a special type of MPI_STATUS object;
            rather, it is a special value for the argument.
            In C one would expect it to be NULL, not the address of a special MPI_STATUS.
            http://mpi-forum.org/docs/mpi-2.0/mpi-20-html/node47.htm
            */

            batch_start = batch_size * i;
            batch_end = batch_start + batch_size - 1;

            if (i == size - 1) {
              //avoid skipping elements when there's a odd size
              batch_end += ns_remainder;
            }

            for (int j = batch_start; j < batch_end; j++) {
              buffer[j] = global_buffer[j];
            }
        }
      } else {
        /*
        MPI_Send Params:
        buf: ponteiro para o dado a ser enviado
        count: num de elementos a serem enviados
        dtype: tipo de elementos do buffer
        dest: rank do destinatário
        tag: num para classificar mensagens
        comm: MPI Communicator
        */
        MPI_Send(buffer, npart, MPI_INT, 0, 0, MPI_COMM_WORLD);
      }

      /*
      Broadcasts a message from the process with rank 0 to all other processes of the group.
      */
      MPI_Bcast(buffer, npart, MPI_INT, 0, MPI_COMM_WORLD);

    /* Generate the initial values */
    InitParticles( particles, pv, npart);
    sim_t = 0.0;

    while (cnt--) {
      double max_f;

      /* Compute forces (2D only) */
      max_f = ComputeForces( particles, particles, pv, npart, batch_start, batch_end );
      /* Compute single max_f for each process and than reduce it to
      global_max_f using MPI_Allreduce */
      printf("BEFORE REDUCE\n");
      printf("Max Force: %f\n", max_f);

      MPI_Allreduce(&max_f, &global_max_f, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

      printf("AFTER REDUCE\n");
      printf("GLOBAL Max Force: %f\n", global_max_f);

      MPI_Barrier(MPI_COMM_WORLD);

      /* Once we have the forces, we compute the changes in position */
      sim_t += ComputeNewPos( particles, pv, npart, global_max_f);
    }
    for (i=0; i<npart; i++)
      fprintf(stdout,"%.5lf %.5lf\n", particles[i].x, particles[i].y);


    MPI_Finalize();
    return 0;
}

void InitParticles( Particle particles[], ParticleV pv[], int npart )
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

double ComputeForces( Particle myparticles[], Particle others[], ParticleV pv[], int npart, int start, int end )
{
  double max_f;
  int i;
  max_f = 0.0;
  for (i = start; i < end; i++) {
    int j;
    double xi, yi, mi, rx, ry, mj, r, fx, fy, rmin;
    rmin = 100.0;
    xi   = myparticles[i].x;
    yi   = myparticles[i].y;
    fx   = 0.0;
    fy   = 0.0;
    for (j=0; j < npart; j++) {
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
    if (fx > max_f) max_f = fx;
  }
  return max_f;
}

double ComputeNewPos( Particle particles[], ParticleV pv[], int npart, double max_f)
{
  int i;
  double a0, a1, a2;
  double dt_new;
  a0	 = 2.0 / (dt * (dt + dt_old));
  a2	 = 2.0 / (dt_old * (dt + dt_old));
  a1	 = -(a0 + a2);
  for (i=0; i<npart; i++) {
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
  if (dt_new < dt) {
    dt_old = dt;
    dt     = dt_new;
  }
  else if (dt_new > 4.0 * dt) {
    dt_old = dt;
    dt    *= 2.0;
  }
  return dt_old;
}
