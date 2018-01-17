#include <stdio.h>

int main (int argc, char **argv)
{
   void datread(char *, void *, int , int );
   void pgmwrite(char *, void *, int , int );

   int i,j,k;
   int m = 600;
   int n = 450;
   int numit = 50000;
  
   float im[m+2][n+2], old[m+2][n+2], new[m+2][n+2];
   float buf[m][n];

   for (i = 0; i < m+2; i++)
     for (j = 0; j < n+2; j++)
       {
	  im[i][j] = 0.0;
	  old[i][j] = 0.0;
       }

   datread("edge600x450.dat", buf, m, n);

   for (i = 1; i <= m; i++)
     for (j = 1; j <= n; j++)
       {
	  im[i][j] = buf[i-1][j-1];
	  old[i][j] = buf[i-1][j-1];
       }

   for (k = 1; k <= numit; k++)
     {
	for (i = 1; i <= m; i++)
	  for (j = 1; j <= n; j++)
	    new[i][j] = 0.25*(old[i-1][j]+old[i+1][j]+old[i][j-1]+old[i][j+1]-
			      im[i][j]);
	for (i = 1; i <= m; i++)
	  for (j = 1; j <= n; j++)
	    old[i][j] = new[i][j];
	if (!(k % 1000)) printf("%d iterations done\n", k);
     }

   for (i = 1; i <= m; i++)
     for (j = 1; j <= n; j++)
       buf[i-1][j-1] = old[i][j];

   pgmwrite("image.pgm", buf, m, n);
    
   return 0;
}
