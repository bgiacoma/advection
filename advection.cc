#include <fstream>
#include <iostream>
#include <math.h>
using namespace std;

double sqr(double x) {return x*x;}

int main(){


  //Ask the user for the domain extension
  cout << "xmax (e.g., 10)? \n";
  double x2;
  cin >> x2;

  //Ask the user how many points to use
  cout << "Number of Points J (e.g., 101)? \n";
  int J;
  cin >> J;

  //Ask the user the value of the courant factor
  cout << "Courant factor cf (e.g., 0.5)? \n";
  double cf;
  cin >> cf;

  //Ask the user at which time we end the simulation
  cout << "Final time tend (e.g., 20)? \n";
  double tend;
  cin >> tend;

  double x[J+2], u[J+2], unpo[J+2], unmo[J+2];
  //note x[0] and x[J+1] are ghost points used for boundary conditions
  const double x1=0.0, x0=5.0;
  double dx = (x2-x1)/(J-1);
  double dt = cf*dx;
  
  double norm=0; //initialize the norm to be zero

  ofstream outs; //output stream
  outs.open("u.dat");
  ofstream norm_out;
  norm_out.open("l2norm.dat");
  norm_out << "Time \t L2norm \n";


  //ask the user to select the initial data
  cout<<"Select initial data \n";
  cout<<"1) Gaussian \n";
  cout<<"2) Step Function \n";
  int initial;
  cin >> initial;

  //ask the user to select the boundary conditions
  cout<<"Select boundary conditions \n";
  cout<<"1) Periodic \n";
  cout<<"2) Outflow \n";
  int boundary;
  cin >> boundary;

  
  //initial conditions
  outs <<"#Time=0" << "\n";
  switch(initial){
  case 1: //Gaussian
    for (int i=0;i<=J+1;i++){
      x[i] = x1 + (i-1)*dx; //spatial grid
      u[i] = exp(-1.0*sqr(x[i]-x0));
      unmo[i] = exp(-1.0*sqr(x[i]+dt-x0));//used by Leapfrog
      unpo[i] = 0.0;
      norm+=sqr(u[i]);
      //saving data
      // (we exclude ghost points from the output)
      if (i>0 && i<J+1) outs << x[i] << "\t" << u[i] << "\n";
    }
    break;
  case 2: //Step Function
    for (int i=0;i<=J+1;i++){
      x[i] = x1 + (i-1)*dx; //spatial grid
      if(x[i]>=4.0 && x[i]<=6.0) 
	u[i] = 1.0;
      else
	u[i] = 0.0;
      if(x[i]>=4.0-dt && x[i]<=6.0-dt) 
	unmo[i] = 1.0;
      else
	unmo[i] = 0.0;
      unpo[i] = 0.0;
      norm+=sqr(u[i]);
      //saving data
      // (we exclude ghost points from the output)
      if (i>0 && i<J+1) outs << x[i] << "\t" << u[i] << "\n";
    }
    break;
  default :
    cerr << "BAD CHOICE IN INTIAL DATA\n";
    return -1;
    break;
  }
  outs << "\n";
  outs << "\n";
  norm=sqrt(norm/J);
  norm_out << 0 << "\t" << norm << "\n";
  norm = 0.0; //reinitialize the norm

  double t=0.0;
  int count=0;

  //ask the user to select the algorithm
  cout<<"Select algorithm \n";
  cout<<"1) FTCS \n";
  cout<<"2) Lax-Friedrichs \n";
  cout<<"3) Lax-Wendroff \n";
  cout<<"4) Leapfrog \n";
  int alg;
  cin >> alg;
  

  if(boundary==1){
    //apply periodic boundary conditions
    u[0]=u[J];
    u[J+1]=u[1];
    unmo[0]=unmo[J];
    unmo[J+1]=unmo[1];}
  else if(boundary==2){
    //apply outflow boundary conditions
    u[0]=u[1];
    u[J+1]=u[J];
    unmo[0]=unmo[1];
    unmo[J+1]=unmo[J];}
  else{
    cerr << "Erroin the choice of boundary conditions";
    return 0;
  }
    

  //evolution
  while (t<tend){    
    cout << "Time=" << t << "\n";
    switch(alg){
    case 1 : //FTCS
      for (int i=1;i<=J;i++){
	unpo[i] = u[i] - dt/(2.0*dx)*(u[i+1]-u[i-1]);
      }
      break;
    case 2 : //Lax-Friedrichs
      for (int i=1;i<=J;i++){
	unpo[i] = 0.5*(u[i+1]+u[i-1]) - dt/(2.0*dx)*(u[i+1]-u[i-1]);
      }
      break;
    case 3 : //Lax-Wendroff
      for (int i=1;i<=J;i++){
	unpo[i] = u[i] - dt/(2.0*dx)*(u[i+1]-u[i-1]) + 
	  dt*dt/(2.0*dx*dx)*(u[i+1]-2.0*u[i]+u[i-1]);
      }
      break;
    case 4 : //Leapfrog
      for (int i=1;i<=J;i++){
	unpo[i] = unmo[i] - dt/dx*(u[i+1]-u[i-1]);
      }
      break;
    default :
      cerr << "BAD CHOICE \n";
      t=tend+1;//stop
      break;
    }

    count++;//increase counter
    t+=dt;//update time

    if(boundary==1){
      //apply periodic boundary conditions
      unpo[0]=unpo[J];
      unpo[J+1]=unpo[1];}
    else if(boundary==2){
      //apply outflow boundary conditions
      unpo[0]=unpo[1];
      unpo[J+1]=unpo[J];
      // The following use a Taylor expansion at first order
      // unpo[0]=2.0*unpo[1]-unpo[2];
      // unpo[J+1]=2.0*unpo[J]-unpo[J-1];
    }
    else{
      cerr << "Erroin the choice of boundary conditions";
      return 0;
    }

    //save data
    if (count==1){
      count=0;
      norm=0;
      outs << "#Time="<< t << "\n";
      for (int i=1;i<=J;i++){
	norm+=sqr(unpo[i]);
	outs << x[i] << "\t" << unpo[i] << "\n";
      }
      outs << "\n";
      outs << "\n";
      norm=sqrt(norm/J);
      norm_out << t << "\t" << norm << "\n";      
    }

    //initialize for new iteration
    for (int i=0; i<=J+1; i++){
      unmo[i]=u[i];
      u[i]=unpo[i];
    }
  }//end of while loop

  //Closing files and terminating the program
  outs.close();
  norm_out.close();
  return 0;

}
