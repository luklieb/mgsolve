#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <errno.h>
#include <sys/time.h>
#include <math.h>
#include <vector>

#define GHOST true

// converts a timeval struct to seconds
double timevalToDouble( struct timeval *t ){
    return (double)(t->tv_sec) + (((double)(t->tv_usec))/1000000.0);
}

double residuum(int nx, int ny, std::vector<double> &grid, std::vector<double> &f_x_y, double h){
    double residuum = 0;
    double sum = 0;
    for(int i=1; i<ny-1; i++){
        for(int k=1; k<nx-1; k++){
                double temp = ((grid[i*nx+(k-1)]+grid[i*nx+(k+1)]) + (grid[(i-1)*nx+k]+grid[(i+1)*nx+k]) - (4*grid[i*nx+k]))/(h*h) + f_x_y[i*nx+k];
                sum += temp*temp;
        }
    }
    residuum = sqrt((1.0/((nx-1)*(ny-1)))*sum);
    return residuum;
}

void residual( int nx, int ny, std::vector<double> &grid, std::vector<double> &f_x_y, std::vector<double> &res, double h){
    for(int i=1; i<ny-1; i++){
        for(int k=1; k<nx-1; k++){
                res[i*nx+k] = ((grid[i*nx+(k-1)]+grid[i*nx+(k+1)]) + (grid[(i-1)*nx+k]+grid[(i+1)*nx+k]) - (4*grid[i*nx+k]))/(h*h) + f_x_y[i*nx+k];
        }
    }
}

// restrict a grid (from) with the full weighting stencile to another grid (to). From has size nx[l] and to has size nx[l-1].
void coarsening( int l, std::vector<double>& from, std::vector<double>& to, std::vector<int>& nx, std::vector<int>& ny ){
    
  for( int i=1; i<nx[l-1]-1; i++ ){
    for( int j=1; j<ny[l-1]-1; j++ ){
      to[i*nx[l-1]+j] = (   from[(2*i+1)*nx[l]+2*j-1] + 2*from[(2*i+1)*nx[l]+2*j] +   from[(2*i+1)*nx[l]+2*j+1]
			+ 2*from[ 2*i   *nx[l]+2*j-1] + 4*from[ 2*i   *nx[l]+2*j] + 2*from[ 2*i*   nx[l]+2*j+1]
			+   from[(2*i-1)*nx[l]+2*j-1] + 2*from[(2*i-1)*nx[l]+2*j] +   from[(2*i-1)*nx[l]+2*j+1] ) / 16.0;
    }
  }
}

void interpolation( int l, std::vector<double>& from, std::vector<double>& to, std::vector<int>& nx, std::vector<int>& ny ){
    for( int i=1; i<nx[l]-1; i++ ){
      for( int j=1; j<ny[l]-1; j++ ){
	if( i%2 == 0 && j%2 == 0 ){
	  // wert uebernehmen
	  to[i*nx[l]+j] = from[(i/2)*nx[l-1]+(j/2)];
	}else if( i%2 != 0 && j%2 != 0 ){
	  // kreuz
	  to[i*nx[l]+j] = (double)( from[(i/2)*nx[l-1]+(j/2)] + from[((i/2)+1)*nx[l-1]+(j/2)] + from[(i/2)*nx[l-1]+(j/2)+1] + from[((i/2)+1)*nx[l-1]+(j/2)+1] ) / 4.0;
	}else if( i%2 == 0 && j%2 != 0 ){
	  // vertikal
	  to[i*nx[l]+j] = (double)(from[(i/2)*nx[l-1]+(j/2)] + from[(i/2)*nx[l-1]+(j/2)+1]) / 2.0;
	}else{
	  // horizontal 
	  to[i*nx[l]+j] = (double)(from[(i/2)*nx[l-1]+(j/2)] + from[((i/2)+1)*nx[l-1]+(j/2)]) / 2.0;
	}
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
	
  if(GHOST){
    for(int i=0; i<ny; i++){
      grid[i*nx] = grid[i*nx+1] - h;
      grid[i*nx+(nx-1)] = grid[i*nx+(nx-1)-1] - h;
    }
  }
}

void multigrid( int l, std::vector<std::vector<double>>& grid, std::vector<std::vector<double>>& f, std::vector<int>& nx,
        std::vector<int>& ny, std::vector<double>& h, std::vector<std::vector<double>>& res, int v1=2, int v2=1, int gamma=1){
  //Presmoothing
  Red_Black_Gauss( nx[l], ny[l], grid[l], f[l], h[l], v1 );
    
  // Residuum
  residual( nx[l], ny[l], grid[l], f[l], res[l], h[l] );
  
  // restrict residual
  coarsening( l, res[l], f[l-1], nx, ny );

  if( l <= 1 ){
    // solve
    Red_Black_Gauss( nx[l-1], ny[l-1], grid[l-1], f[l-1], h[l-1], 1 );
  
  }else{
    for( int i=1; i<nx[l-1]-1; i++ ){
      for( int j=1; j<ny[l-1]-1; j++ ){
	grid[l-1][i*nx[l-1]+j] = 0.0;
      }
    }
    for( int i=0; i<gamma; i++ ){
      multigrid( l-1, grid, f, nx, ny, h, res, v1, v2, gamma );
    }
  
  
  // interpolation
    std::vector<double> c( nx[l]*ny[l], 0.0 );
    interpolation( l, grid[l-1], c, nx, ny );
   
    // corretion
    for( int i=1; i<nx[l]-1; i++ ){
      for( int j=1; j<ny[l]-1; j++ ){
	grid[l][i*nx[l]+j] += c[i*nx[l]+j];
      }
    }
  }
  

  
  //Postsmothing
  Red_Black_Gauss( nx[l], ny[l], grid[l], f[l], h[l], v2 );
}

void testInterpolation(){
	std::vector<double> to( { 0.0, 0.0, 0.0, 0.0, 0.0,
					  0.0, 0.0, 0.0, 0.0, 0.0,
					  0.0, 0.0, 0.0, 0.0, 0.0,
					  0.0, 0.0, 0.0, 0.0, 0.0, 
					  0.0, 0.0, 0.0, 0.0, 0.0 } );
	std::vector<double> from( { 0.0, 0.0, 0.0,
				    0.0, 25.0, 0.0,
				    0.0, 0.0, 0.0 } );

	std::vector<int> nx( {3, 5} );
	std::vector<int> ny( {3, 5} );
	interpolation( 1, from, to, nx, ny );
	
	for( int i=0; i<5; i++ ){
		for( int j=0; j<5; j++ ){
			fprintf( stderr, "%lf, ", to[i*5+j]);
		}
		fprintf( stderr, "\n" );
	}
}

void testCoarsening(){
	std::vector<double> from( { 1.0, 1.0, 1.0, 1.0, 1.0,
					  2.0, 2.0, 2.0, 2.0, 2.0,
					  3.0, 3.0, 3.0, 3.0, 3.0,
					  4.0, 4.0, 4.0, 4.0, 4.0, 
					  5.0, 5.0, 5.0, 5.0, 5.0 } );
	std::vector<double> to( { 0.0, 0.0, 0.0,
				    0.0, 0.0, 0.0,
				    0.0, 0.0, 0.0 } );

	std::vector<int> nx( {3, 5} );
	std::vector<int> ny( {3, 5} );
	coarsening( 1, from, to, nx, ny );
	
	for( int i=0; i<3; i++ ){
		for( int j=0; j<3; j++ ){
			fprintf( stderr, "%lf, ", to[i*3+j]);
		}
		fprintf( stderr, "\n" );
	}
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
	
	// nx and ny are the total points in x and y direction
    std::vector<int> nx( l, 0 );
    std::vector<int> ny( l, 0 );
    std::vector<double>  h( l, 0 );
    for( int i=l-1; i>=0; i-- ){
      nx[i] = (int)(pow(2,i+1)+1);
      ny[i] = (int)(pow(2,i+1)+1);
      h[i] = 1.0/(nx[i]-1);
    }
    double convergence = 0.0;
    double old_residuum = 0.0;
    double new_residuum = 0.0;

	// Bei Neumann Randwerte werden in x-Richtung ghost-Layers benoetigt
    if(GHOST){
      for( int i=0; i<l; i++)
        nx[i]=nx[i]+2;
    }
	
	// Anlegen von f_x_y:
    std::vector<std::vector<double>> f_x_y(l);
    for( int i=l-1; i>=0; i--){
      if(GHOST){
	f_x_y[i] = std::vector<double>(nx[i]*ny[i],2.0);
      }else{
	f_x_y[i] = std::vector<double>(nx[i]*ny[i],0.0);
      }
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
	  grid[i] = std::vector<double>(nx[i]*ny[i], 0.0);
	}
	
	// Dirichlet RDB setzten aber nicht auf den Ghost-Layers 
	for( int j=l-1; j>=l-1; j--){
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
        for(int i=0; i<nx[l-1]; i++){
            grid[l-1][nx[l-1]*(ny[l-1]-1)+i] = sin(M_PI*i*h[l-1])*sinh(M_PI);
        
	}
    }
    
	// TIMER RED_BLACK_GAUSS
    struct timeval t1;
    gettimeofday(&t1, NULL);



    for(int j=0; j<n; j++){
        multigrid( l-1, grid, f_x_y, nx, ny, h, res, n);
        new_residuum = residuum(nx[l-1], ny[l-1], grid[l-1], f_x_y[l-1], h[l-1]);
        fprintf(stdout, "L2 Norm: %lf\n", new_residuum);
	if(j>0){
            convergence = new_residuum / old_residuum;
	    fprintf(stdout, "Convergence: %lf\n", convergence);
        }
        old_residuum = new_residuum;
    }


    struct timeval t2;
    gettimeofday(&t2, NULL);
    fprintf(stdout, "Timer Multigrid: %lf\n", timevalToDouble(&t2)-timevalToDouble(&t1));

    // Oeffnen der Ausgabedatei
    FILE *out = fopen("./solution.txt", "w");
	if(out == NULL){
		perror("fopen");
		exit(EXIT_FAILURE);
    }

  double errorSum = 0.0;
    // Ausgabe fuer solution.txt
    fprintf(out, "# x y u(x,y)\n");
    for(int j=0; j<ny[l-1]; j++){
        for(int t=0; t<nx[l-1]; t++){
            double x = (double)t/(double)(nx[l-1]-1);
            double y = (double)j/(double)(ny[l-1]-1);
	    double temp = grid[l-1][j*ny[l-1]+t] - sin(x*M_PI)*sinh(y*M_PI);
	    errorSum += temp * temp;
	    
            fprintf(out, "%.5lf %.5lf %.8lf\n", x, y, grid[l-1][j*ny[l-1]+t]);
        }
	fprintf(out, "\n");
    }
    // Ausgabe fuer error.txt
    fprintf(out, "\n" );
    fclose(out);

    
    double error = sqrt((1.0/((nx[l-1]-1)*(ny[l-1]-1)))*errorSum);
    fprintf( stdout, "Fehler zur korrekten Lsg: %lf\n", error );
    
}
