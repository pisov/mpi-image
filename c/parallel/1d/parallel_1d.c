#include <stdio.h>
#include <mpi.h>

int main (int argc, char **argv)
{
   void datread(char *, void *, int , int );
   void pgmwrite(char *, void *, int , int );

   int i,j,k;
   int m = 239;
   int n = 432;
   int numit = 5000;

   int size, rank;
   int up, down;
   int dims[1],periods[1];
   MPI_Comm MPI_COMM_OneD;
   int srow, erow, mcounts;
   int sendcnts[16], displs[16];
// Declare row type variable name;
   MPI_Datatype MPI_ONE_ROW;

 
   float im[m+2][n+2], old[m+2][n+2], new[m+2][n+2];
   float buf[m][n];

   MPI_Init(&argc, &argv);
   MPI_Comm_size(MPI_COMM_WORLD,&size);
   dims[0] = size;
   periods[0] = 0;
   MPI_Cart_create(MPI_COMM_WORLD,1,dims,periods,0,&MPI_COMM_OneD);
   MPI_Comm_rank(MPI_COMM_OneD,&rank);
   MPI_Cart_shift(MPI_COMM_OneD,0,1,&up,&down);

// Define new type
   MPI_Type_contiguous(n+2,MPI_FLOAT,&MPI_ONE_ROW);
   MPI_Type_commit(&MPI_ONE_ROW);

   srow = (m * rank) / size + 1;
   erow = (m * (rank + 1)) / size;
   mcounts = erow - srow + 1;

   for (i = 0; i < m+2; i++)
     for (j = 0; j < n+2; j++)
       {
	  im[i][j] = 0.0;
	  old[i][j] = 0.0;
       }

   if (rank == 0) datread("edge239x432.dat", buf, m, n);

   for (i = 0; i < size; i++) {
      displs[i] = ((m * i) / size) * n;
      sendcnts[i] = ((m * (i + 1)) / size - (m * i) / size) * n;
      //printf(" scnt[%d] = %d dspl[%d] = %d ",i,displs[i],i,sendcnts[i]);
   }

   if (rank != 0)
      MPI_Scatterv(buf, sendcnts, displs, MPI_FLOAT, buf, mcounts * n, MPI_FLOAT, 0, MPI_COMM_OneD);
   else
      MPI_Scatterv(buf, sendcnts, displs, MPI_FLOAT, MPI_IN_PLACE, mcounts * n, MPI_FLOAT, 0, MPI_COMM_OneD);



   for (i = 1; i <= mcounts; i++)
     for (j = 1; j <= n; j++)
       {
	  im[i][j] = buf[i-1][j-1];
	  old[i][j] = buf[i-1][j-1];
       }

   for (k = 1; k <= numit; k++)
     {
/*
	MPI_Sendrecv(&old[1][0]                  ,n+2,MPI_FLOAT,up  ,0,
                     &old[mcounts+1][0],n+2,MPI_FLOAT,down,0,
                     MPI_COMM_OneD,MPI_STATUS_IGNORE);
	MPI_Sendrecv(&old[mcounts][0],n+2,MPI_FLOAT,down,0,
                     &old[0][0]        ,n+2,MPI_FLOAT,up  ,0,
                     MPI_COMM_OneD,MPI_STATUS_IGNORE);
*/

        MPI_Sendrecv(&old[1][0]      ,1,MPI_ONE_ROW,  up,0,&old[mcounts+1][0],1,MPI_ONE_ROW,down,0,MPI_COMM_OneD,MPI_STATUS_IGNORE);
        MPI_Sendrecv(&old[mcounts][0],1,MPI_ONE_ROW,down,0,&old[0][0],        1,MPI_ONE_ROW,up  ,0,MPI_COMM_OneD,MPI_STATUS_IGNORE);

	for (i = 1; i <= mcounts; i++)
	  for (j = 1; j <= n; j++)
	    new[i][j] = 0.25*(old[i-1][j]+old[i+1][j]+old[i][j-1]+old[i][j+1]-
			      im[i][j]);
	for (i = 1; i <= mcounts; i++)
	  for (j = 1; j <= n; j++)
	    old[i][j] = new[i][j];
	if (rank == 0)
		if (!(k % 10000)) printf("%d iterations done\n", k);
     }

   for (i = 1; i <= mcounts; i++)
     for (j = 1; j <= n; j++)
       buf[i-1][j-1] = old[i][j];

   if (rank != 0)
      MPI_Gatherv(buf, mcounts * n, MPI_FLOAT, buf, sendcnts, displs, MPI_FLOAT, 0, MPI_COMM_OneD);
   else
      MPI_Gatherv(MPI_IN_PLACE, mcounts * n, MPI_FLOAT, buf, sendcnts, displs, MPI_FLOAT, 0, MPI_COMM_OneD);


   if (rank == 0) pgmwrite("image.pgm", buf, m, n);

   MPI_Finalize();
    
   return 0;
}
