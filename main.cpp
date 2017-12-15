#include "functions.h"

using namespace std;

typedef double real;
typedef int whole;

int main(){
	try{
		//To choose the method for solving the model
		int method;
		real delx, delt;
		whole n;
		

		
		do{
			do{
				fflush(stdin);
				cout<<"\nEnter 1 for DeFort-Frenkel method ";
				cout<<"\nEnter 2 for Richardson method ";
				cout<<"\nEnter 3 for Laasonen Simple Implicit method ";
				cout<<"\nEnter 4 for Crank Nickolson method ";
				cout<<"\nEnter 5 for Analytical solution";
				cout<<"\nEnter 0 to exit \n";
				cin>>method;
				
				
			}
			while(method < 0 || method > 5);
			

			switch (method){
				case 1:
				
					//Taking input from user
					input(delx, delt);
				
					//Total nodes
					n = 1/delx+1;
				
					dufort(delx, delt, n);
					break;
				
				case 2:
				
					input(delx, delt);
				
					n = 1/delx+1;
				
					richardson(delx, delt, n);
					break;
				
				case 3:
					
					input(delx, delt);
					
					n = 1/delx+1;
					
					laasonen(delx, delt, n);
					break;
					
				case 4:
					
					input(delx, delt);
					
					n = 1/delx+1;
					
					crank(delx, delt, n);
					break;
					
				case 5:
					input(delx, delt);
					
					n = 1/delx + 1;
					
					analytical(n);
					break;
			}
		}
		while(method);	
		
	}
	catch (exception &e){
		cerr<<"Your error is: "<<endl;
		cerr<<"Type: "<<typeid(e).name()<<endl;
		cerr<<"What: "<<e.what()<<endl;
	}
	

	
}
