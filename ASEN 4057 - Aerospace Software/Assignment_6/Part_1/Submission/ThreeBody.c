#include <stdio.h>
#include <stdlib.h>
#include <math.h>


struct Initals {
    double mM;
    double mS;
    double mE;
    double rM;
    double rE;
    double G;

    double dES0;
    double xS0;
    double yS0;
    double vSx0;
    double vSy0;

    double dEM0;
    double vM;
    double xM0;
    double yM0;
    double vMx0;
    double vMy0;

    //double data[][8];

};
// -----------------------------------------------------------------------------
// Function Prototypes


int eulerFunc(double *data, int time, int dt,double* usr_x,int* ptr);

void eomFunc(double* stateVars, double *stateDers);

void GridSearch(struct Initals OG,double* usr_x);

double OptimizationFunction(double delta_V_Sx, double delta_V_Sy,struct Initals OG, double* usr_x, double * data, int time, int dt, int* ptr);

void saveFile(double* data, int index, int time, int dt);

// -----------------------------------------------------------------------------
// Main
int main(int argc, char const *argv[]){

// Input argument formating
double usr_x[4] = {atof(argv[0]), atof(argv[1]), atof(argv[2])}; // For actual useage
//double usr_x[3] = {1, 100000, 0.1}; // For Testing/Debugging
// Iniitial Conditions
struct Initals C; // Initial Variables


// Population of Structure values
C.mM = 7.34767309E22;
C.mE = 5.97219E24;
C.mS = 28833;
C.rM = 1737100;
C.rE = 6371000;
C.G = 6.674E-11;

C.dES0 = 340000000;
C.xS0 = C.dES0 * cos(50* M_PI / 180);
C.yS0 = C.dES0 * sin(50* M_PI / 180);
C.vSx0 = 1000 *  cos(50* M_PI / 180);
C.vSy0 = 1000 *  sin(50* M_PI / 180);

C.dEM0 = 384403000;
C.vM = sqrt((C.G*pow(C.mE,2)) / ((C.mE + C.mM)*C.dEM0));
C.xM0 = C.dEM0 * cos(42.5* M_PI / 180);
C.yM0 = C.dEM0 * sin(42.5* M_PI / 180);
C.vMx0 = -C.vM *  sin(42.5* M_PI / 180);
C.vMy0 = C.vM *  cos(42.5* M_PI / 180);

//data[0] = C.xS0;
//data[1] = C.yS0;
//data[2] = C.vSx0;
//data[3] = C.vSy0;

//data[4] = C.xM0;
//data[5] = C.yM0;
//data[6] = C.vMx0;
//data[7] = C.vMy0;

GridSearch(C,usr_x);

//printf("Dad?");

}
// -----------------------------------------------------------------------------
// Grid Search Function


void GridSearch(struct Initals OG,double* usr_x){
    /*Initial guess, arbitrary but the closer it is, the
    faster the grid search is function works */
    double delta_V_S_01 = 1000;
    double delta_V_S_02 = 1000;
    // Initializing a few Variables
    double delta_V_X, delta_V_Y,v,minX,minY;

    int time = 1E6; // Time Scale for integration
    int dt = 10; // Time scale resolution
    int index;
    double* data = (double*)malloc((8*(time/dt) +8) * sizeof(*data)); // data array for position and velocites of spacecraft and moon {xS,yS,vxS,vyS,xM,yM,vxM,vyM}

    // Grid search to further narrow down the initial guess (step size of 10)
    double max_delta_V_S = 100; // max delta V as per constraints


if (usr_x[0] == 1){
  /* Upper hemisphere */
  /* Upper hemisphere */
  for (delta_V_X = -100; delta_V_X <= 100; delta_V_X+=10){
            for (delta_V_Y = 0; delta_V_Y <= 100; delta_V_Y +=10){
          if (sqrt(pow(delta_V_X,2) + pow(delta_V_Y,2)) < max_delta_V_S){
              v = OptimizationFunction(delta_V_X, delta_V_Y, OG ,usr_x,data,time,dt,&index);
              if (v < 100){
                  max_delta_V_S = v;
              }
              if (v < sqrt(pow(delta_V_S_01,2)+pow(delta_V_S_02,2))){
                  delta_V_S_01 = delta_V_X;
                  delta_V_S_02 = delta_V_Y;
              }
          }
      }
  }
  // test
  /* Lower hemisphere */
  for (delta_V_X = -100; delta_V_X <= 100; delta_V_X+=10){
      for (delta_V_Y = 0; delta_V_Y >= -100; delta_V_Y -= 10){
          if (sqrt(pow(delta_V_X,2) + pow(delta_V_Y,2)) < max_delta_V_S){
              v = OptimizationFunction(delta_V_X, delta_V_Y, OG ,usr_x, data,time, dt, &index);
              if (v < 100){
                  max_delta_V_S = v;
              }
              if (v < sqrt(pow(delta_V_S_01,2)+pow(delta_V_S_02,2))){
                  delta_V_S_01 = delta_V_X;
                  delta_V_S_02 = delta_V_Y;
              }
          }
      }
  }

  minX =   delta_V_S_01;
  minY =   delta_V_S_02;

  for (double i = minX-1;i <= minX+1;i+=usr_x[2]){
    for (double j = minY-1; j<= minY+1;j+=usr_x[2]){
      v = OptimizationFunction(i, j, OG ,usr_x, data,time, dt, &index);
      if (v < sqrt(pow(delta_V_S_01,2)+pow(delta_V_S_02,2))){
          delta_V_S_01 = i;
          delta_V_S_02 = j;
      }
    }
  }
}
else if(usr_x[0] == 2){

  int time_stp = 0;
  int min_time = time;


  /* Upper hemisphere */
  for (delta_V_X = -100; delta_V_X <= 100; delta_V_X+=10){
            for (delta_V_Y = 0; delta_V_Y <= 100; delta_V_Y +=10){
          if (sqrt(pow(delta_V_X,2) + pow(delta_V_Y,2)) < 100){
              time_stp = OptimizationFunction(delta_V_X, delta_V_Y, OG ,usr_x,data,time,dt,&index);
              if (time_stp < min_time){
                  min_time = time_stp;
                  delta_V_S_01 = delta_V_X;
                  delta_V_S_02 = delta_V_Y;
              }
          }
      }
  }
  // test
  /* Lower hemisphere */


  for (delta_V_X = -100; delta_V_X <= 100; delta_V_X+=10){
      for (delta_V_Y = 0; delta_V_Y >= -100; delta_V_Y -= 10){
          if (sqrt(pow(delta_V_X,2) + pow(delta_V_Y,2)) < 100){
              time_stp = OptimizationFunction(delta_V_X, delta_V_Y, OG ,usr_x, data,time, dt, &index);
              if (time_stp < min_time){
                  min_time = time_stp;
                  delta_V_S_01 = delta_V_X;
                  delta_V_S_02 = delta_V_Y;

              }
          }
      }
  }

  minX =   delta_V_S_01;
  minY =   delta_V_S_02;

  for (double i = minX -1; i <= minX+1; i += usr_x[2]){
            for (double j = minY-1; j <= minY+1; j += usr_x[2]){
          if (sqrt(pow(i,2) + pow(j,2)) < 100){
              time_stp = OptimizationFunction(i, j, OG ,usr_x,data,time,dt,&index);
              if (time_stp < min_time){
                  min_time = time_stp;
                  delta_V_S_01 = i;
                  delta_V_S_02 = j;
              }
          }
      }
  }
}
else{printf("Haha you thought! you silly goose, its gotta be a 1 or 2!");}


v = OptimizationFunction(delta_V_S_01, delta_V_S_02, OG ,usr_x, data,time, dt, &index);

printf("Smallest delta V required to return to Earth : %f\n",sqrt(pow(delta_V_S_01,2)+pow(delta_V_S_02,2)));
printf("dv_01 = %f\n",delta_V_S_01);
printf("dv_02 = %f\n",delta_V_S_02);
saveFile(data,index, time, dt);


}

// -----------------------------------------------------------------------------

// Optimzation Function
double OptimizationFunction(double delta_V_Sx, double delta_V_Sy,struct Initals OG, double* usr_x, double * data, int time, int dt, int* ptr){
  double result;
  struct Initals D = OG;


  //RESET INITIAL CONDITIONS WITH ADDED DELTA V
  //IC = [xS_0, yS_0, vSx_0 + delta_V_S(1), vSy_0+ delta_V_S(2), xM_0, yM_0, vMx_0, vMy_0, xE_0, yE_0, vEx_0, vEy_0]';%vector of initial conditions
  D.vSx0 += delta_V_Sx;
  D.vSy0 += delta_V_Sy;
  //Use eulers function to integrate eomFunc with new delta_V_S
  //INTEGRATE HERE

  data[0] = D.xS0;
  data[1] = D.yS0;
  data[2] = D.vSx0;
  data[3] = D.vSy0;

  data[4] = D.xM0;
  data[5] = D.yM0;
  data[6] = D.vMx0;
  data[7] = D.vMy0;

  // if simulation results in the spaceship returning back to earth, return
  // the total delta V that resulted in that successful simulation, otherwise
  // return very high value
  if (eulerFunc(data,time,dt,usr_x,ptr)==1){
      //return magnitude of added delta_V_Sx and delta_V_Sy
      if (usr_x[0] == 1){
      result = sqrt(pow(delta_V_Sx,2) + pow(delta_V_Sy,2));
      }
      else if (usr_x[0] == 2){
        result = *ptr/8 -8;
      }
    }
  else
    if (usr_x[0] == 1){result=100000;}
    else if(usr_x[0] == 2){result = 2E6;}

return result;
}

// -----------------------------------------------------------------------------
// Euler Fnction
int eulerFunc(double *data, int time, int dt ,double* usr_x,int*ptr){
// Allocation of function variables
  // Array of state variable derivatives
  double stateDers[8];
  // Array of state varaibles
  double stateVars[8];
  // Positioning and distance variables
  double xS,yS,xM,yM,dES,dEM,dMS;

  int timespan = 8*time/dt;


  for (int i = 8; i< timespan; i=i+8){
        // Deifnes values of state varaibles for derivites in equation of motion function
        for (int j=0; j<8; j++){
          stateVars[j] = data[i+j-8];
        }

        // Calculates the derivites for the current state variables
        eomFunc(stateVars,stateDers);

        // Preforms Euler method of intengration
        data[i] = stateVars[0]  + dt * stateDers[0];
        data[i+1] = stateVars[1]  + dt * stateDers[1];
        data[i+2] = stateVars[2]  + dt * stateDers[2];
        data[i+3] = stateVars[3]  + dt * stateDers[3];

        data[i+4] = stateVars[4]  + dt * stateDers[4];
        data[i+5] = stateVars[5]  + dt * stateDers[5];
        data[i+6] = stateVars[6]  + dt * stateDers[6];
        data[i+7] = stateVars[7]  + dt * stateDers[7];

      //  printf("data = {%e,%e,%e,%e,%e,%e,%e,%e,}\n\n",data[i],data[i+1],data[i+2],data[i+3],data[i+4],data[i+5],data[i+6],data[i+7]);
      // Collects spacecraft and Moon positions to determine status
        xS =  data[i];
        yS =  data[i+1];
        xM =  data[i+4];
        yM =  data[i+5];
        dES = sqrt(pow(xS,2)+pow(yS,2));
        dMS = sqrt(pow(xS-xM,2)+pow(yS-yM,2));
        dEM = sqrt(pow(xM,2) + pow(yM,2));

        //printf("dMS = %f\n",dMS);

      // Checks if termination conditions have been meet

        // If spacecraft is at Earth surface, teminate with 1 status (Earth rendezvous)
        if(dES <= 6371000){*ptr = i;return 1;}
        // If spacecraft is at user specified distance from Moon surface, terminate with 2 status (Lunar distance threshold)
        else if ( dMS <= (1737100+usr_x[1])){*ptr = i;return 2;}
        // If spacecraftis at a distance from Earth greater than 2x the Moon-Earth distance, terminate with 3 status (Lost in Space)
        else if (dES >= (2*dEM)){*ptr = i;return 3;}

        }
        // If spacecraft has not hit a termination condition, exit with 0 status (Time period exhausted)
      return 0;}


      void eomFunc(double *P_state, double *P_stateDers){
      //set up constants
              double G = 6.674E-11; // Gravitational constant
              double mM = 7.34767309E22; /*mass of moon in kg*/
              double mE = 5.97219E24;  /*mass of Earth kg*/
              double mS = 28833; /*mass of spacecraft in kg*/
              //double rM = 1737100; /*radius of moon in m*/
              //double rE = 6371000; /*radius of the Earth in m*/


      //Extract from input vec
              double xS = P_state[0];
              double yS = P_state[1];
              double vSx = P_state[2];
              double vSy = P_state[3];
              double xM = P_state[4];
              double yM = P_state[5];
              double vMx = P_state[6];
              double vMy = P_state[7];


      // Force Equations //

      //Force from Earth to Moon (F_EM)
              double d_EM = sqrt(pow((xM),2) + pow((yM),2));


              double F_EMx=(G*mE*mM*(-xM))/(pow(d_EM,3)); //force on Moon by Earth in the x-direction
              double F_EMy=(G*mE*mM*(-yM))/(pow(d_EM,3)); // force on Moon by Earth in the y-direction

      //Force from Satellite to Moon (F_SM)
              double d_SM = sqrt(pow((xS-xM),2) + pow((yS-yM),2));

              double F_MSx=(G*mS*mM*(xM-xS))/(pow(d_SM,3)); // force on Satellite by Moon in the x-direction
              double F_MSy=(G*mS*mM*(yM-yS))/(pow(d_SM,3)); // force on Satellite by Moon in the y-direction


      //Force from Satellite to Earth (F_ES)
              double d_ES = sqrt(pow((xS),2) + pow((yS),2));

              double F_ESx=(G*mS*mE*(-xS))/(pow(d_ES,3)); // force on Satellite by Moon in the x-direction
              double F_ESy=(G*mS*mE*(-yS))/(pow(d_ES,3)); // force on Satellite by Moon in the y-direction

      // Acceleration Equations //

      //Newton's second law and acceleration
              double aSx=(F_MSx+F_ESx)/mS;
              double aSy=(F_MSy+F_ESy)/mS;
              double aMx=(F_EMx-F_MSx)/mM;
              double aMy=(F_EMy-F_MSy)/mM;

      // Output derivatives
              P_stateDers[0] = vSx;
              P_stateDers[1] = vSy;
              P_stateDers[2] = aSx;
              P_stateDers[3] = aSy;
              P_stateDers[4] = vMx;
              P_stateDers[5] = vMy;
              P_stateDers[6] = aMx;
              P_stateDers[7] = aMy;
              //printf("vSx = %e\n",vSx);
             // double output_derivative = [vSx, vSy, aSx, aSy, vMx, vMy, aMx, aMy, vEx, vEy, aEx, aEy];
      return;
}

// -----------------------------------------------------------------------------
// saveFile function

void saveFile(double* data, int index, int time, int dt){

  //printf("Index = %d\n",index);
  FILE *file;

  file = fopen("ThreeBody_data.csv","w");
  // File formating
  fprintf(file,"xS  yS  vSx   vSy   xM   yM   vMx   vMy\n");
  // Writes the data into the csv file

    for (int i =0; i < index; i=i+8){
      //printf("Im Here\n");
      fprintf(file,"%e   %e   %e   %e   %e   %e   %e   %e\n", data[i],data[i+1],data[i+2],data[i+3],data[i+4],data[i+5],data[i+6],data[i+7]);
    }


  //closes the file
  fclose(file);
  free(data);
}
