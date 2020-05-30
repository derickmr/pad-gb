/*
  This program is an adaptation of the Mandelbrot program
  from the Programming Rosetta Stone, see
  http://rosettacode.org/wiki/Mandelbrot_set
  Compile the program with:
  gcc -o mandelbrot -O4 mandelbrot.c
  Usage:
 
  ./mandelbrot <xmin> <xmax> <ymin> <ymax> <maxiter> <xres> <out.ppm>
  Example:
  ./mandelbrot 0.27085 0.27100 0.004640 0.004810 1000 1024 pic.ppm
  The interior of Mandelbrot set is black, the levels are gray.
  If you have very many levels, the picture is likely going to be quite
  dark. You can postprocess it to fix the palette. For instance,
  with ImageMagick you can do (assuming the picture was saved to pic.ppm):
  convert -normalize pic.ppm pic.png
  The resulting pic.png is still gray, but the levels will be nicer. You
  can also add colors, for instance:
  convert -negate -normalize -fill blue -tint 100 pic.ppm pic.png
  See http://www.imagemagick.org/Usage/color_mods/ for what ImageMagick
  can do. It can do a lot.
*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <stdint.h>
#include <pthread.h>

unsigned char *colorsToBeWrittenOnFile;

typedef struct{
    int yres;
    int xres;
    double xmin;
    double xmax;
    double ymin;
    double ymax;
    uint16_t maxiter;
    int arrayCounter;
    int threadStart;
    int threadEnd;
}thread_arg, *ptr_thread_arg;


void *calculate_mandelbrot(void *arg){
    
    ptr_thread_arg targ = (ptr_thread_arg)arg;
    
    /* Precompute pixel width and height. */
     double dx=(targ->xmax-targ->xmin)/targ->xres;
     double dy=(targ->ymax-targ->ymin)/targ->yres;
    
    double x, y; /* Coordinates of the current point in the complex plane. */
    double u, v; /* Coordinates of the iterated point. */
    int i,j; /* Pixel counters */
    int k; /* Iteration counter */
    
    int counter = targ->arrayCounter;
        
    for (j = targ->threadStart; j < targ->threadEnd; j++) {
      y = targ->ymax - j * dy;
      for(i = 0; i < targ->xres; i++) {
      
        double u = 0.0;
        double v= 0.0;
        double u2 = u * u;
        double v2 = v*v;
        x = targ->xmin + i * dx;
        /* iterate the point */
        for (k = 1; k < targ->maxiter && (u2 + v2 < 4.0); k++) {
              v = 2 * u * v + y;
              u = u2 - v2 + x;
              u2 = u * u;
              v2 = v * v;
        };
        /* compute  pixel color and write it to file */
        if (k >= targ->maxiter) {
          /* interior */
            int colorCounter;
            for (colorCounter = 0; colorCounter < 6; colorCounter++){
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
  if (argc != 8) {
    printf("Usage:   %s <xmin> <xmax> <ymin> <ymax> <maxiter> <xres> <out.ppm>\n", argv[0]);
    printf("Example: %s 0.27085 0.27100 0.004640 0.004810 1000 1024 pic.ppm\n", argv[0]);
    exit(EXIT_FAILURE);
  }
    
    pthread_t threads[4];

  /* The window in the plane. */
    double xmin = atof(argv[1]);
    double xmax = atof(argv[2]);
    double ymin = atof(argv[3]);
    double ymax = atof(argv[4]);

  /* Maximum number of iterations, at most 65535. */
   uint16_t maxiter = (unsigned short)atoi(argv[5]);

  /* Image size, width is given, height is computed. */
  int xres = atoi(argv[6]);
  int yres = (xres*(ymax-ymin))/(xmax-xmin);
    
  int arraySize = yres * xres * 6;
    
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

    
    int i, j;
    
    int numThreads = 4;
    
    thread_arg arguments[numThreads];
 
    for (i = 0; i < numThreads; i++){
        arguments[i].yres = yres;
        arguments[i].xres = xres;
        arguments[i].xmin = xmin;
        arguments[i].xmax = xmax;
        arguments[i].ymin = ymin;
        arguments[i].ymax = ymax;
        arguments[i].maxiter = maxiter;
        arguments[i].threadCounter = (arraySize/numThreads) * i;
        arguments[i].threadStart = (yres/numThreads) * i;
        arguments[i].threadEnd = (yres/numThreads) * (i+1);
    }
    
    arguments[numThreads-1] += yres&numThreads;

    for (i = 0; i < numThreads; i++){
        pthread_create(&(threads[i]), NULL, calculate_mandelbrot, &(arguments[i]));
    }
    
    for (i = 0; i < numThreads; i++){
        pthread_join(threads[i], NULL);
    }
        
    unsigned char color[6];
    
    for (i = 0; i < arraySize; ){
        for (j = 0; j < 6; j++){
                color[j] = colorsToBeWrittenOnFile[i++];
        }
        fwrite(color, 6, 1, fp);
    }
        
  fclose(fp);
  return 0;
}
