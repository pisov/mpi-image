#include <stdio.h>
#include <mpi.h>

int main (int argc, char **argv)
{
   void datread(char *, void *, int , int );
   void pgmwrite(char *, void *, int , int );

   int i,j,k;
   int m = 600;
   int n = 450;
   int numit = 100000;

   int size, rank;
   int up, down;
   int left, right;
   int dims[2],periods[2], crd[2];
   MPI_Comm MPI_COMM_TwoD;
   int srow, erow, mcounts;
   int scol, ecol, ncounts;
   int sendcnts[32], displs[32];
// Declare row type variable name;
   MPI_Datatype MPI_ONE_ROW, MPI_ONE_COL, MPI_BLOCK2, MPI_BLOCK;
 
   float im[m+2][n+2], old[m+2][n+2], new[m+2][n+2];
   float buf[m][n];

   double time;

   MPI_Init(&argc, &argv);
   MPI_Comm_size(MPI_COMM_WORLD,&size);
// Get optimal grid
   dims[0] = 0;
   dims[1] = 0;
   MPI_Dims_create(size,2,dims);
   periods[0] = 0;
   periods[1] = 0;
   MPI_Cart_create(MPI_COMM_WORLD,2,dims,periods,0,&MPI_COMM_TwoD);
   MPI_Comm_rank(MPI_COMM_TwoD,&rank);
// Get directions in up/down
   MPI_Cart_shift(MPI_COMM_TwoD,0,1,&up,&down);
// Get directions in left/right
   MPI_Cart_shift(MPI_COMM_TwoD,1,1,&left,&right);

// Calculate stripes
   mcounts = m / dims[0];
   ncounts = n / dims[1];

   if (rank == 0) {
       printf("mcount: %d ncount %d\n",mcounts,ncounts);
   }
   
// Define new type ROW
   MPI_Type_contiguous(ncounts+2,MPI_FLOAT,&MPI_ONE_ROW);
   MPI_Type_commit(&MPI_ONE_ROW);
// Define new type COLUMN
   MPI_Type_vector(mcounts+2,1,n+2,MPI_FLOAT,&MPI_ONE_COL);
   MPI_Type_commit(&MPI_ONE_COL);
// Define BLOCK type
   MPI_Type_vector(mcounts,ncounts,n,MPI_FLOAT,&MPI_BLOCK2);
   MPI_Type_create_resized(MPI_BLOCK2,0,sizeof(float),&MPI_BLOCK);
   MPI_Type_commit(&MPI_BLOCK);

   for (i = 0; i < m+2; i++)
     for (j = 0; j < n+2; j++)
       {
	  im[i][j] = 0.0;
	  old[i][j] = 0.0;
       }

   if (rank == 0) datread("edge600x450.dat", buf, m, n);

   for (k = 0; k < size; k++) {
       MPI_Cart_coords(MPI_COMM_TwoD, k, 2, crd);
       i = crd[0];
       j = crd[1];
       displs[k] = i * n * mcounts + j * ncounts;
       sendcnts[k] = 1; 
   }

//   if (rank == 0)
//      for (i=0; i < dims[0]*dims[1];i++)
//         printf("%d\t",displs[i]);
//   printf("\n");

   if (rank != 0)
      MPI_Scatterv(buf, sendcnts, displs, MPI_BLOCK, buf, 1, MPI_BLOCK, 0, MPI_COMM_TwoD);
   else
      MPI_Scatterv(buf, sendcnts, displs, MPI_BLOCK, MPI_IN_PLACE, 1, MPI_BLOCK, 0, MPI_COMM_TwoD);



   for (i = 1; i <= mcounts; i++)
     for (j = 1; j <= ncounts; j++)
       {
	  im[i][j] = buf[i-1][j-1];
	  old[i][j] = buf[i-1][j-1];
       }

    time = MPI_Wtime();
   for (k = 1; k <= numit; k++)
     {

        MPI_Sendrecv(&old[1][0]      ,1,MPI_ONE_ROW,  up,0,&old[mcounts+1][0],1,MPI_ONE_ROW,down,0,MPI_COMM_TwoD,MPI_STATUS_IGNORE);
        MPI_Sendrecv(&old[mcounts][0],1,MPI_ONE_ROW,down,0,&old[0][0],        1,MPI_ONE_ROW,up  ,0,MPI_COMM_TwoD,MPI_STATUS_IGNORE);

        MPI_Sendrecv(&old[0][1]      ,1,MPI_ONE_COL,left ,0,&old[0][ncounts+1],1,MPI_ONE_COL,right,0,MPI_COMM_TwoD,MPI_STATUS_IGNORE);
        MPI_Sendrecv(&old[0][ncounts],1,MPI_ONE_COL,right,0,&old[0][0]        ,1,MPI_ONE_COL,left ,0,MPI_COMM_TwoD,MPI_STATUS_IGNORE);

	for (i = 1; i <= mcounts; i++)
	  for (j = 1; j <= ncounts; j++)
	    new[i][j] = 0.25*(old[i-1][j]+old[i+1][j]+old[i][j-1]+old[i][j+1]-
			      im[i][j]);
	for (i = 1; i <= mcounts; i++)
	  for (j = 1; j <= ncounts; j++)
	    old[i][j] = new[i][j];
	if (!(k % 1000) && rank == 0) printf("%d iterations done\n", k);
     }
    time = MPI_Wtime() - time;

   for (i = 1; i <= mcounts; i++)
     for (j = 1; j <= ncounts; j++)
       buf[i-1][j-1] = old[i][j];

   if (rank != 0)
      MPI_Gatherv(buf, 1, MPI_BLOCK, buf, sendcnts, displs, MPI_BLOCK, 0, MPI_COMM_TwoD);
   else
      MPI_Gatherv(MPI_IN_PLACE, 1, MPI_BLOCK, buf, sendcnts, displs, MPI_BLOCK, 0, MPI_COMM_TwoD);


   if (rank == 0) { 
       pgmwrite("image.pgm", buf, m, n);
       printf("Execution time: %.2f s\n",time);
    }

   MPI_Finalize();
    
   return 0;
}
