#include "functions.h"

void input(double &del_x, double &del_t){
	
	std::cout<<"Enter the value of delta x: "<<std::endl;
	std::cin>>del_x;
	std::cout<<"Enter the value of delta t: "<<std::endl;
	std::cin>>del_t;
	
}

void write_file(Vector V, const char *filename){
	std::ofstream plot;
	plot.open(filename, std::ios::app);
	plot<<V;
	plot.close();
}

void dufort(double dx, double dt, int n){
	//Remove the output file if it already exists
	remove("dufort.plot");

	//Factors used for DuFort-Frenkel Method
	const double ALPHA = 2*D*dt/dx/dx;
	double a, b;
	a = (1-ALPHA)/(1+ALPHA); //a and b are constants used in DeFort Frenkel equation
	b = ALPHA/(1+ALPHA);
	
	//Integers to be used in for loops
	int i, j, k;
	
	//Total time is 0.5hrs
	double tot_time = 0.5;
	
	//Total time steps
	int TIME_STEPS = tot_time/dt;

	//Temporary temperature vectors for DuFort-Frenkel method
	Vector t(n), temp_n(n), temp_n_minus(n);
	
	//For loop for time steps start here	
	for (i=0; i<TIME_STEPS; i++){
		
		//For loop for updating temperature at each node
		for (j=1; j<n-1; j++){
			t[j] = (a*temp_n_minus[j]) + (b*(temp_n[j+1]+temp_n[j-1]));
		}
		
		//For loop to update temperature after one timestep
		for (k=1; k<n-1; k++){
			temp_n_minus[k] = temp_n[k];
			temp_n[k] = t[k];
		}
		
		//Write to output file
		if (i == 9 || i==19 || i==29 || i==39 || i==49)
		write_file(t, "dufort.plot");
	}
	write_file(t, "dufort.plot");
	std::cout<<"\nOut put has been stored in file named 'DuFort.plot' \n";
}

void richardson(double dx, double dt, int n){
	//Remove the output file if it already exists
	remove("richardson.plot");

	//Factors used for DuFort-Frenkel Method
	const double ALPHA = 2*D*dt/dx/dx;
	
	//Integers to be used in for loops
	int i, j, k;
	
	//Total time is 0.5hrs
	double tot_time = 0.5;
	
	//Total time steps
	int TIME_STEPS = tot_time/dt;

	//Temporary temperature vectors for DuFort-Frenkel method
	Vector t(n), t_n(n), t_n_minus(n);
	
	//For loop for time steps start here	
	for (i=0; i<TIME_STEPS; i++){
		
		//For loop for updating temperature at each node
		for (j=1; j<n-1; j++){
			t[j] = t_n_minus[j] + (ALPHA*(t_n[j+1]-(2*t_n[j])+t_n[j-1]));
		}
		
		//For loop to update temperature after one timestep
		for (k=1; k<n-1; k++){
			t_n_minus[k] = t_n[k];
			t_n[k] = t[k];
		}
		
		//Write to output file
		if (i == 9 || i==19 || i==29 || i==39 || i==49)
		write_file(t, "richardson.plot");
	}
	std::cout<<"\nOut put has been stored in file named 'richardson.plot' \n";
}

void laasonen(double dx, double dt, double n){
	//Remove the output file if it already exists
	remove("Laasonen.plot");

	//Factors used for DuFort-Frenkel Method
	const double ALPHA = D*dt/dx/dx;
	
	//Total time is 0.5hrs
	double tot_time = 5;
	
	//Total time steps
	int TIME_STEPS = tot_time/dt;
	
	Matrix smatrix(n,n);
	
	int i, j;
	smatrix[0][0] = 1;
	smatrix[n-1][n-1] = 1;
	
	for (i=1; i<n-1; i++){
		for(j=1; j<n-1; j++){
			if(i==j){
				smatrix[i][j] = 1 + (2*ALPHA);
				smatrix[i][j-1] = -ALPHA;
				smatrix[i][j+1] = -ALPHA;
			}
			else smatrix[i][j] == 0;
		}
	}
	
	Matrix p(n,n);
	
	Vector x(n), b(n);
	
	Matrix l(n,n), u(n,n);
	
	//reorder(smatrix, n, p);
	
	//pa = p*smatrix;
	// factorise
	lu_fact(smatrix,l,u,n);
		
	for(i=0; i<TIME_STEPS; i++){
		lu_solve(l,u,b,n,x);
		b = x;
		
		//Write to output file
		//if (i == 9 || i==19 || i==29 || i==39 || i==49)
		//if (i == 3 || i==7 || i==11 || i==15 || i==19)
		//if (i == 39 || i==79 || i==119 || i==159 || i==199)
		//write_file(b, "Laasonen.plot");
	}
	write_file(b, "Laasonen.plot");

	std::cout<<"\nOut put has been stored in file named 'richardson.plot' \n";
	
}

void crank(double dx, double dt, double n){
	//Remove the output file if it already exists
	remove("CrankNicklson.plot");

	//Factors used for DuFort-Frenkel Method
	const double ALPHA = D*dt/dx/dx/2;
	
	//Total time is 0.5hrs
	double tot_time = 0.5;
	
	//Total time steps
	int TIME_STEPS = tot_time/dt;
	
	Matrix smatrix(n,n);
	
	int i, j;
	smatrix[0][0] = 1;
	smatrix[n-1][n-1] = 1;
	
	for (i=1; i<n-1; i++){
		for(j=1; j<n-1; j++){
			if(i==j){
				smatrix[i][j] = 1 + (2*ALPHA);
				smatrix[i][j-1] = -ALPHA;
				smatrix[i][j+1] = -ALPHA;
			}
			else smatrix[i][j] == 0;
		}
	}
	
	Matrix p(n,n);
	
	Vector x(n), b(n);
	
	Matrix l(n,n), u(n,n);
	
	//reorder(smatrix, n, p);
	
	//pa = p*smatrix;
	// factorise
	lu_fact(smatrix,l,u,n);
	Vector temp(n);
		
	for(i=0; i<TIME_STEPS; i++){
		for (j=1; j<n-1; j++){
			temp[j] = (ALPHA*b[j+1]) + ((1-(2*ALPHA))*b[j]) +(ALPHA*b[j-1]);
		}
		lu_solve(l,u,temp,n,x);
		b = x;
		
		//Write to output file
		if (i == 9 || i==19 || i==29 || i==39 || i==49)
		write_file(b, "CrankNicklson.plot");
	}
	
	std::cout<<"\nOut put has been stored in file named 'richardson.plot' \n";
}

void analytical(double n){
	remove("analytical.plot");
	
	int i,j, k, m;
	double sum, time=0.1;
	
	Vector t(n), dx(n);
	
	for(i=0; i<n; i++){
		dx[i] = i*0.05;
	}
	
	for (i=0; i<5; i++){
		for(j=0; j<n; j++){
			sum = 0;
			for(m=1; m<10; m++){
				sum = sum + (exp(-D*(m*3.14)*(m*3.14)*time)*((1-pow(-1,m))/m/3.14)*sin(m*3.14*dx[j]));
			}
			t[j] = 300 + (2*(-200)*sum);
		}
		write_file(t, "analytical.plot");
		time += 0.1;
	}
	
}
