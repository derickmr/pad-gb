/*
  This program is an adaptation of the Mandelbrot program
  from andrejbauer project on github, see:
  https://gist.github.com/andrejbauer/7919569
 
  It uses andrejbauer algorithm as a basis, but processing is done in parallel using pthread library, in a master-slave model.
 
  Compile the program with:
  gcc -pthread -o mandelbrot mandelbrot.c
  
  Usage:
  ./mandelbrot <xmin> <xmax> <ymin> <ymax> <maxiter> <xres> <out.ppm> <numThreads>
  
  Example:
  ./mandelbrot 0.27085 0.27100 0.004640 0.004810 1000 1024 pic.ppm 10
 
  The resulting file can be opened in any software that handles .ppm extensions (e.g. Gimp).
  
*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <stdint.h>
#include <pthread.h>
#include <time.h>
#define COLOR_SIZE 6

unsigned char *colorsToBeWrittenOnFile;

typedef struct{
    int xres;
    double xmin;
    double ymin;
    double ymax;
    uint16_t maxiter;
    int counter;
    int threadStart;
    int threadEnd;
    double dx;
    double dy;
}thread_arg, *ptr_thread_arg;


void *calculate_mandelbrot(void *arg){
    
    ptr_thread_arg targ = (ptr_thread_arg)arg;
    
    int xres = targ->xres;
    double xmin = targ->xmin;
    double ymin = targ->ymin;
    double ymax = targ->ymax;
    uint16_t maxiter = targ->maxiter;
    int counter = targ->counter;
    int threadStart = targ->threadStart;
    int threadEnd = targ->threadEnd;
    double dx = targ->dx;
    double dy = targ->dy;
    
    printf ("xres: %d \n xmin: %f \n ymin: %f \n counter: %d \n threadStart: %d \n threadEnd: %d \n dx: %f \n dy: %f \n\n", xres, xmin, ymin, counter, threadStart, threadEnd, dx, dy);
    
    double x, y; /* Coordinates of the current point in the complex plane. */
    double u, v; /* Coordinates of the iterated point. */
    int i,j; /* Pixel counters */
    int k; /* Iteration counter */
    
        
    for (j = threadStart; j < threadEnd; j++) {
      y = ymax - j * dy;
      for(i = 0; i < xres; i++) {
      
        double u = 0.0;
        double v= 0.0;
        double u2 = u * u;
        double v2 = v*v;
        x = xmin + i * dx;
        /* iterate the point */
        for (k = 1; k < maxiter && (u2 + v2 < 4.0); k++) {
              v = 2 * u * v + y;
              u = u2 - v2 + x;
              u2 = u * u;
              v2 = v * v;
        };
        /* compute pixel color*/
        if (k >= maxiter) {
          /* interior */
            int colorCounter;
            for (colorCounter = 0; colorCounter < COLOR_SIZE; colorCounter++){
                colorsToBeWrittenOnFile[counter++] = 0;
            }
          
        }
        else {
          /* exterior */
            colorsToBeWrittenOnFile[counter++] = k >> 8;
            colorsToBeWrittenOnFile[counter++] = k & 255;
            colorsToBeWrittenOnFile[counter++] = k >> 8;
            colorsToBeWrittenOnFile[counter++] = k & 255;
            colorsToBeWrittenOnFile[counter++] = k >> 8;
            colorsToBeWrittenOnFile[counter++] = k & 255;

        }
      }
    }
    return NULL;
  }
    

int main(int argc, char* argv[])
{
  /* Parse the command line arguments. */
  if (argc != 9) {
    printf("Usage:   %s <xmin> <xmax> <ymin> <ymax> <maxiter> <xres> <out.ppm> <numthreads>\n", argv[0]);
    printf("Example: %s 0.27085 0.27100 0.004640 0.004810 1000 1024 pic.ppm 10\n", argv[0]);
    exit(EXIT_FAILURE);
  }
    
    clock_t begin = clock();

    
  /* The window in the plane. */
    double xmin = atof(argv[1]);
    double xmax = atof(argv[2]);
    double ymin = atof(argv[3]);
    double ymax = atof(argv[4]);
    int numThreads = atoi(argv[8]);
    
    pthread_t threads[numThreads];


  /* Maximum number of iterations, at most 65535. */
   uint16_t maxiter = (unsigned short)atoi(argv[5]);

  /* Image size, width is given, height is computed. */
  int xres = atoi(argv[6]);
  int yres = (xres*(ymax-ymin))/(xmax-xmin);
    
  int arraySize = yres * xres * COLOR_SIZE;
    
  colorsToBeWrittenOnFile = (unsigned char *)malloc(arraySize * sizeof(unsigned char));
    
  printf ("array size: %d \n", arraySize);

  /* The output file name */
  const char* filename = argv[7];

  /* Open the file and write the header. */
  FILE * fp = fopen(filename,"wb");
  char *comment="# Mandelbrot set";/* comment should start with # */

  /*write ASCII header to the file*/
  fprintf(fp,
          "P6\n# Mandelbrot, xmin=%lf, xmax=%lf, ymin=%lf, ymax=%lf, maxiter=%d\n%d\n%d\n%d\n",
          xmin, xmax, ymin, ymax, maxiter, xres, yres, (maxiter < 256 ? 256 : maxiter));

    /* Precompute pixel width and height. */
    double dx=(xmax-xmin)/xres;
    double dy=(ymax-ymin)/yres;
        
    thread_arg arguments[numThreads];
 
    int i, j;
    
    printf ("test 1\n");
    
    //Initializing threads struct
    for (i = 0; i < numThreads; i++){
        arguments[i].xres = xres;
        arguments[i].xmin = xmin;
        arguments[i].ymin = ymin;
        arguments[i].ymax = ymax;
        arguments[i].maxiter = maxiter;
        arguments[i].dx = dx;
        arguments[i].dy = dy;
        arguments[i].counter = (arraySize/numThreads) * i;
        arguments[i].threadStart = (yres/numThreads) * i;
        arguments[i].threadEnd = (yres/numThreads) * (i+1);
    }
    
    printf ("test 2\n");
    
    arguments[numThreads-1].threadEnd += yres%numThreads;

    //Computing slaves
    for (i = 1; i < numThreads; i++){
        pthread_create(&(threads[i]), NULL, calculate_mandelbrot, &(arguments[i]));
    }
    
    printf ("test 3\n");
    
    //Master
    calculate_mandelbrot(&(arguments[0]));
    
    for (i = 1; i < numThreads; i++){
        pthread_join(threads[i], NULL);
    }
    
    printf ("test 4\n");
        
    //Writing result to file
    unsigned char color[COLOR_SIZE];
    
    for (i = 0; i < arraySize; ){
        for (j = 0; j < COLOR_SIZE; j++){
                color[j] = colorsToBeWrittenOnFile[i++];
        }
        fwrite(color, COLOR_SIZE, 1, fp);
    }
    
    printf ("test 5\n");
        
  fclose(fp);
  free(colorsToBeWrittenOnFile);
     
    clock_t end = clock();

    double time_spent = (double)(end - begin) / CLOCKS_PER_SEC;
    
    printf ("Execution time: %f \n", time_spent);

  return 0;
}
