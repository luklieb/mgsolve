#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <errno.h>
#include <sys/time.h>
#include <math.h>
#include <vector>

#define GHOST false

// converts a timeval struct to seconds
double timevalToDouble( struct timeval *t ){
    return (double)(t->tv_sec) + (((double)(t->tv_usec))/1000000.0);
}

double residuum(int nx, int ny, std::vector<double> &grid, std::vector<double> &f_x_y, double h){
    double residuum = 0;
    double sum = 0;
    for(int i=1; i<ny-1; i++){
        for(int k=1; k<nx-1; k++){
                double temp = ((-1.0f/h)*(grid[i*nx+(k-1)]+grid[i*nx+(k+1)]) - (1.0f/h)*(grid[(i-1)*nx+k]+grid[(i+1)*nx+k]) + ((2.0f/h+2.0f/h)*grid[i*nx+k]) - f_x_y[i*nx+k]);
                sum += temp*temp;
        }
    }
    residuum = sqrt((1.0f/(nx*ny))*sum);
    return residuum;
}

void residual( int l, int nx, int ny, std::vector<double> &grid, std::vector<double> &f_x_y, std::vector<double> &res, double h){
    for(int i=1; i<ny-1; i++){
        for(int k=1; k<nx-1; k++){
                res[i*nx+k] = ((-1.0/h)*(grid[i*nx+(k-1)]+grid[i*nx+(k+1)]) - (1.0/h)*(grid[(i-1)*nx+k]+grid[(i+1)*nx+k]) + ((2.0/h+2.0/h)*grid[i*nx+k]) - f_x_y[i*nx+k]);
        }
    }
}

// restrict a grid (from) with the full weighting stencile to another grid (to). From has size nx[l] and to has size nx[l-1].
void coarsening( int l, std::vector<double>& from, std::vector<double>& to, std::vector<double>& nx, std::vector<double>& ny ){
    
  for( int i=0; i<nx[l-1]; i++ ){
    for( int j=0; j<ny[l-1]; j++ ){
      to[i*nx[l-1]+j] = (   from[(2*i+1)*nx[l]+2*j-1] + 2*from[(2*i+1)*nx[l]+2*j] +   from[(2*i+1)*nx[l]+2*j+1]
			+ 2*from[2*i*nx[l]+2*j-1]     + 4*from[2*i*nx[l]+2*j]     + 2*from[2*i*nx[l]+2*j+1]
			+   from[(2*i-1)*nx[l]+2*j-1] + 2*from[(2*i-1)*nx[l]+2*j] +   from[(2*i-1)*nx[l]+2*j+1] ) /16.0;
    }
  }
}


void Red_Black_Gauss(int nx, int ny, std::vector<double> &grid, std::vector<double> &f_x_y, double h, int numIterations){
    
     for(int iterations=0; iterations<numIterations; iterations++){
		for(int m=1; m<ny-1; m++){
			int q=1;
			if(m%2 == 0){
				 q++;
			}
			for(; q<nx-1; q=q+2){
				 grid[m*nx+q] = (1.0/4.0) * (h*h*f_x_y[m*nx+q] + (grid[m*nx+(q-1)]+grid[m*nx+(q+1)]+grid[(m-1)*nx+q]+grid[(m+1)*nx+q]));
			}
		}
		for(int m=1; m<ny-1; m++){
			int q=1;
			if(m%2 != 0){
			      q++;
			}
			for(; q<nx-1; q=q+2){
				grid[m*nx+q] = (1.0/4.0) * (h*h*f_x_y[m*nx+q] + (grid[m*nx+(q-1)]+grid[m*nx+(q+1)]+grid[(m-1)*nx+q]+grid[(m+1)*nx+q]));
			}
		}
	}
}

void multigrid( int l, std::vector<std::vector<double>>& grid, std::vector<std::vector<double>>& f, std::vector<double>& nx,
		std::vector<double>& ny, std::vector<double>& h, std::vector<std::vector<double>>& res, int v1=2, int v2=1, int gamma=1 ){
  
  // l = 0 entspricht Level 1 
  l--;
  
  //Presmoothing
  Red_Black_Gauss( nx[l], ny[l], grid[l], f[l], h[l], v1 );
  
  // Residuum
  residual( l, nx[l], ny[l], grid[l], f[l], res[l], h[l] );
  
  // restrict residual
  coarsening( l, res[l], f[l-1], nx, ny );
  
  if( l == 1 ){
    Red_Black_Gauss( nx[l-1], ny[l-1], grid[l-1], f[l-1], h[l-1], 10 );
  }else{
    for( int i=0; i<gamma; i++ ){
      multigrid( l, grid, f, nx, ny, h, res, v1, v2, gamma );
    }
    // interpolation
    std::vector<double> c( nx[l]*ny[l], 0.0 );
    interpolate( c2h to c );
    for( int i=0; i<nx[l-1]; i++ ){
      for( int j=0; j<ny[l-1]; j++ ){
	grid[l][i*nx[l]+j] += c[i*nx[l]+j];
      }
    }
  }
  
  //Postsmothing
  Red_Black_Gauss( nx[l], ny[l], grid[l], f[l], h[l], v2 );
}

int main(int argc, char **argv){

	// Ueberpruefung, ob Eingabeparamter passen
    if(argc != 3){
		fprintf(stderr, "Usage: ./mgsolve l n\n");
		exit(EXIT_SUCCESS);
	}
	
    // Einlesen der Eingabeparameter l, n
    int l, n;
    l = atoi(argv[1]);
    n = atoi(argv[2]);
    fprintf(stderr, "%d/n", n);
	
	// nx and ny are the total points in x and y direction
    std::vector<int> nx( l, 0 );
    std::vector<int> ny( l, 0 );
    std::vector<double>  h( l, 0 );
    for( int i=l-1; i>=0; i-- ){
      nx[i] = (int)(pow(2,i)+1);
      ny[i] = (int)(pow(2,i)+1);
      h[i] = 1.0/(nx[i]-1);
    }

	// Bei Neumann Randwerte werden in x-Richtung ghost-Layers benoetigt
    if(GHOST){
      for( int i=0; i<l; i++)
        nx[i]=nx[i]+2;
    }
	
	// Anlegen von f_x_y:
    std::vector<std::vector<double>> f_x_y(l);
    for( int i=l-1; i>=0; i--){
	f_x_y[i] = std::vector<double>(nx[i]*ny[i],0.0);
    }
    
    // Speicher fuer das residual allokieren
    std::vector<std::vector<double>> res(l);
    for( int i=l-1; i>=0; i--){
      res[i] = std::vector<double>(nx[i]*ny[i], 0.0);
    }
  
    // Speicher fuer das Gitter allokieren:
    std::vector<std::vector<double>> grid(l);
    
    // Initialisierung des Gitters
    if(GHOST){
        for( int i=l-1; i>=0; i--){
	  grid[i] = std::vector<double>(nx[i]*ny[i], 2.0);
	}
	
	// Dirichlet RDB setzten aber nicht auf den Ghost-Layers 
	for( int j=l-1; j>=0; j--){
	  for(int i=0; i<nx[j]-2; i++){
	      grid[j][i+1] = i*h[j]*(1-i*h[j]);
	      grid[j][nx[j]*(ny[j]-1)+i+1] = i*h[j]*(1-i*h[j]);
	  }
	  for(int i=0; i<ny[j]; i++){
	      grid[j][i*nx[j]] = -h[j];
	      grid[j][i*nx[j]+(nx[j]-1)] = -h[j];
	  }
	}
    }
    else{
      for( int i=l-1; i>=0; i--){
	grid[i] = std::vector<double>(nx[i]*ny[i],0.0);
      }
      for( int j=l-1; j>=0; j--){
        for(int i=0; i<nx[j]; i++){
            grid[j][nx[j]*(ny[j]-1)+i] = sin(M_PI*i*h[j])*sinh(M_PI);
        }
      }
    }

	// TIMER RED_BLACK_GAUSS
    struct timeval t1;
    gettimeofday(&t1, NULL);

    Red_Black_Gauss(nx[l-1], ny[l-1], grid[l-1], f_x_y[l-1], h[l-1], 1000);

    struct timeval t2;
    gettimeofday(&t2, NULL);
    fprintf(stdout, "Timer Red_Black_Gauss: %lf\n", timevalToDouble(&t2)-timevalToDouble(&t1));

    double resid = residuum(nx[l-1], ny[l-1], grid[l-1], f_x_y[l-1], h[l-1]);

    // Oeffnen der Ausgabedatei
	FILE *out = fopen("./solution.txt", "w");
	if(out == NULL){
		perror("fopen");
		exit(EXIT_FAILURE);
    }

    // Ausgabe fuer solution.txt
    fprintf(out, "# x y u(x,y)\n");
    for(int j=0; j<ny[l-1]; j++){
        for(int t=0; t<nx[l-1]; t++){
            double f = (double)t/(double)(nx[l-1]-1);
            double q = (double)j/(double)(ny[l-1]-1);

            fprintf(out, "%.5lf %.5lf %.8lf\n", f, q, grid[l-1][j*ny[l-1]+t]);
        }
	fprintf(out, "\n");
    }
    fclose(out);
    fprintf(stdout, "Residuum: %lf\n", resid);
}
