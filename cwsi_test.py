from cwsi_function import CWSIModel, ModelsInputData


def main():
    """
    Main function to test the CWSI  and transpiration models.
    """

    # Input data to the CWSI model
    data = ModelsInputData()
    data.Ta = 25.5
    data.Sgl = 800.02
    data.RH = 60.60
    data.U = 2.3
    data.b2 = 3.0
    data.b0 = 0.1
    data.lat = 42.3
    data.Elev = 2500
    data.DOY = 180    

    # Put data into the CWSI model
    models = CWSIModel(data)

    # Canopy temperature
    Tc = 24.1  

    # Calculate actual transpiration
    t_a = models.T_Actual_F(Tc)
    # Calculate potential transpiration
    t_p = models.Tp_F()       
    print("*********** Transpiration Models Test **********")
    print("Actual Transpiration = {:.1f} (mm/day)".format(t_a))
    print("Potential Transpiration = {:.1f} (mm/day)".format(t_p))
    print()
 
    # Calculate CWSI
    cwsi = models.CWSI(Tc) 
    print("**************** CWSI Model Test  ****************")
    print("CWSI (cal) = {:.6f}".format(cwsi))
    expected_cwsi = 0.216767
    print("CWSI (exp) = {:.6f}".format(expected_cwsi))

    # Check if the two values are approximately equal
    if abs(expected_cwsi - cwsi) < 0.001:
        print("\nTest passed!")

        # Calculate CWSI for a range of canopy temperatures
        print("\nNum\t\t Tc\t\t CWSI")
        for i in range(20):
            Tc += 0.1
            cwsi = models.CWSI(Tc)
            print(f"{i}\t\t {Tc:.1f}\t\t {cwsi:.2f}")        
    else:
        print("\nTest did NOT pass!")

    print("******************   End   *********************")
    
if __name__ == "__main__":
    main()