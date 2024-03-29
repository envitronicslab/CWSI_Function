#include <stdio.h>
#include <math.h>       /* pow */ /* exp */ /* abs */ 
#include "CWSIFunction.h"



int main(int argc, char * argv[])
{
    /* Input data to the CWSI model */
    models_input_data data;
    data.Ta = 25.5;    
    data.Sgl = 800.02;
    data.RH = 60.60;
    data.U = 2.3;
    data.b2 = 3.0;
    data.b0 = 0.1;
    data.lat = 42.3;    
    data.Elev = 2500;    
    double Tc {24.1};
    double DOY {180};

    /* Put data into the model */
    CWSIModel cwsi_model(data); 

    /* Calculate CWSI: */
    double cwsi = cwsi_model.CWSI(Tc);

    printf("*********** CWSI Model Test **********\n");
    
    /* Calculated output: */
    printf("CWSI (cal) = %f\n", cwsi);
    /* Expected output: */
    double CWSI {0.216767};
    printf("CWSI (exp) = %f\n", CWSI);
    /* Check if the two values are approximately equal */
    if (abs(CWSI - cwsi) < 0.001) {
        printf("\nTest passed!\n\n");
    } else {
        printf("\nTest did NOT pass!\n\n");
    }

    /* Calculate CWSI for a range of canopy temperatures*/
    printf("Num    Tc     CWSI\n");
    for (int i = 0; i < 20; i++)
    {
        Tc += 0.1;
        cwsi = cwsi_model.CWSI(Tc);
        if (i < 10) 
        {
            printf("%d      %.1f   %.2f\n", i, Tc, cwsi);
        } else {
            printf("%d     %.1f   %.2f\n", i, Tc, cwsi);
        }
    }

    printf("************** End **************\n");

    return 0;
}