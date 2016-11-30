#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>

unsigned int compute_pi(unsigned int, unsigned int, int, int);

int main(int argc, char **argv){
  unsigned int pontos;
  unsigned int i, size, rank;
  unsigned int global_points;
  unsigned int pontos_no_circulo;

  if(argc != 2){
    printf("Uso:\n");
    printf("\t%s <numero de pontos a serem sorteados>\n", argv[0]);
    return 1;
  }

  pontos = atoi(argv[1]);

  MPI_Init(&argc, &argv);
  MPI_Comm_size(MPI_COMM_WORLD, &size);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);


  // retorna quantos pontos sorteados cairam dentro do circulo
  // aqui estamos considerando uma semente para o gerador de
  // numeros aleatorios fixa = 0

  int batch_size = (pontos/size);
  int batch_start = (rank*batch_size);
  int batch_end = batch_start + batch_size;

  pontos_no_circulo = compute_pi(0, pontos, batch_start, batch_end);
  MPI_Allreduce(&pontos_no_circulo, &global_points, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);

  // calcula a aproximacao de Pi baseado nos pontos sorteados
  if(rank == 0){
    printf("Pi = %.040f\n", ((double)global_points/(double)pontos)*4);
  }

  MPI_Finalize();
  return 0;
}

unsigned int compute_pi(unsigned int seed, unsigned int pontos, int start, int end){
  unsigned int i;
  unsigned int pontos_no_circulo;
  double x, y;

  pontos_no_circulo = 0;
  srand(seed);

  for(i = start; i < end; i++){
  	// sorteia um ponto: coordenadas x e y dentro do quadrado
  	// consideramos que R = 1, então x e y pertencem ao intervalo [0; 1]
    x = (double)rand()/(double)(RAND_MAX);
    y = (double)rand()/(double)(RAND_MAX);

    // verifica se o ponto sorteado encontra-se dentro do circulo
    // um ponto (x, y) esta dentro do circulo se: x^2 + y^2 < R^2
    // nesse caso, consideramos R = 1
    if( (x*x + y*y) < 1 ){
      pontos_no_circulo++;
    }
  }
  return pontos_no_circulo;
}
