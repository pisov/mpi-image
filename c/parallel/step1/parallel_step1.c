#include <stdio.h>
#include <mpi.h>

int main (int argc, char **argv)
{
   void datread(char *, void *, int , int );
   void pgmwrite(char *, void *, int , int );

   int i,j,k,mloc;
   int m = 600;
   int n = 450;
   int numit = 50000;

   int size, rank, root;
   int scount[16], displs[16];
   float im[m+2][n+2], old[m+2][n+2], new[m+2][n+2];
   float buf[m][n];


   MPI_Init(&argc, &argv);
   
   MPI_Comm_size(MPI_COMM_WORLD, &size);
   MPI_Comm_rank(MPI_COMM_WORLD, &rank);

   root = 0;

   for (i = 0; i < size; i++) {
     scount[i] = (i+1)*m/size - i*m/size; 
     scount[i] *= n;
     displs[i] = (i*m/size)*n;
   }

   mloc = scount[rank]/n;
   for (i = 0; i < mloc+2; i++)
     for (j = 0; j < n+2; j++)
       {
	  im[i][j] = 0.0;
	  old[i][j] = 0.0;
       }

   if (rank == root) {
     datread("edge600x450.dat", buf, m, n);
   }

   MPI_Scatterv(buf, scount, displs, MPI_FLOAT, buf, scount[rank], MPI_FLOAT, root, MPI_COMM_WORLD);

   for (i = 1; i <= mloc; i++)
     for (j = 1; j <= n; j++)
       {
	  im[i][j] = buf[i-1][j-1];
	  old[i][j] = buf[i-1][j-1];
       }

   for (k = 1; k <= numit; k++)
     {
	for (i = 1; i <= mloc; i++)
	  for (j = 1; j <= n; j++)
	    new[i][j] = 0.25*(old[i-1][j]+old[i+1][j]+old[i][j-1]+old[i][j+1]-
			      im[i][j]);
	for (i = 1; i <= mloc; i++)
	  for (j = 1; j <= n; j++)
	    old[i][j] = new[i][j];
	if (!(k % 1000)) printf("%d iterations done\n", k);
     }

   for (i = 1; i <= mloc; i++)
     for (j = 1; j <= n; j++)
       buf[i-1][j-1] = old[i][j];

   MPI_Gatherv(buf, scount[rank], MPI_FLOAT, buf, scount, displs, MPI_FLOAT, root, MPI_COMM_WORLD);

   if (rank == root) {
     pgmwrite("image.pgm", buf, m, n);
   }

   MPI_Finalize();    
   return 0;
}
