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
    data.DOY = 180;  
    double Tc = 24.1;    

    /* Put data into the model */
    CWSIModel models(data); 
    
    /* Calculate actual transpiration */ 
    double t_a = models.T_Actual_F(Tc);
    /* Calculate potential transpiration */ 
    double t_p = models.Tp_F();       
    printf("*********** Transpiration Models Test **********\n");
    printf("Actual Transpiration = %.1f (mm/day)\n", t_a);
    printf("Potential Transpiration = %.1f (mm/day)\n", t_p);
    printf("\n");

    /* Calculate CWSI: */
    double cwsi = models.CWSI(Tc);

    printf("***********  CWSI Model Test  **********\n");
    
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
    printf("Num\t\tTc\t\tCWSI\n");
    for (int i = 0; i < 20; i++)
    {
        Tc += 0.1;
        cwsi = models.CWSI(Tc);
        if (i < 10) 
        {
            printf("%d\t\t%.1f\t\t%.2f\n", i, Tc, cwsi);
        } else {
            printf("%d\t\t%.1f\t\t%.2f\n", i, Tc, cwsi);
        }
    }

    printf("**************    End    **************\n");

    return 0;
}