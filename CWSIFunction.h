#ifndef MYCWSIMODEL_
#define MYCWSIMODEL_



// Models input data
struct models_input_data {
    /* Variables loaded w/ default values */
    double Ta {25.5};  // Air Temperature (deg C) -> measured
    // double Tc;   // Canopy Temperature (deg C) -> measured
    double RH {60.6};  // Relative Humidity (%) -> measured
    double Sgl {800.02}; // Incoming short wave solar radiation (W/m^2) -> estimated
    double U {2.3};   // Wind Speed (m/s) -> measured
    double m {1};   // Calibration coeff.
    double b1 {0};  // Calibration coeff.
    double b2 {3.0}; // Calibration Coefficient -> default value = 8.00, 4.00 for Daylight
    double b0 {0.1}; // Calibration Coefficient -> default value = 0.00
    double lat {42.3}; // Latitude (degree)
    double Elev {2500.0}; // Elevation (m)
    int DOY {180}; // Day of Year (Julian Day)
};


class CWSIModel {
public:
    CWSIModel(models_input_data &data): data_(data) {}
    ~CWSIModel() {};

    // Function to calculate the crop water stress index (CWSI)
    double CWSI(double Tc);

    // Potential transpiration function
    double Tp_F();

    // Non-Linearized actual transpiration function
    double T_Actual_F(const double &Tc);  

private:
    // Main microclimate data
    models_input_data data_; 

    // Model constants
    static constexpr double GT_MAX {99};    // gT_Max, mmol/m2/s
    static constexpr double  GT_MIN {0.0};         //gT_Min, 0.00; mmol/m2/s        
    static constexpr double  DTL_MIN {-20};  // dTl_Min, deg C        
    static constexpr double  ET_MAX {99};  // ET_Max, mm
    static constexpr double  T_ACTUAL_MIN {0.001};  // T_Actual_Min, mm        
    static constexpr double  CD {0.036};  // Cd, = 0.72 * 0.05  // To Leaf Width = 5 cm
    static constexpr double  WIND_F {1.0};  // Wind_f, 0.3 // WS(In-Canopy)/WS(Roza) at the middel of distance between the center of caopy and free air        
    static constexpr double  ALBEDO_LEAF {0.09};  // Albedo_leaf
    static constexpr double  ALBEDO_LEAF_G {0.1};  // Albedo_leaf_g, Albedo of Ground = Thermal Absortptivity of Ground
    static constexpr double  TRANS {0.06};  // Trans
    static constexpr double  ALPHA_LW {0.95};  // alpha_lw, Longwave absorptivity of apple leaf        
    static constexpr double  CP_AIR {29.17};  // Cp_Air, Heat capacity of air J/mol/C
    static constexpr double  Y {44000};  // y, J/mol
    static constexpr double  SIGMA {5.671e-8};  // sigma, Stefan-Boltzmann constant (W/m^2/K^4)        
    static constexpr double  FE {1.0};     // Fe, 1.00; Fe is the view factor between the entire surface (only one side) of the leaf and the complete sphere of the view
    // static constexpr double PI {3.14159265359};     // Already defined in math lib

    // Sub Functions
    double Rnet(double Tc);
    double dTu_F();
    double dTp_F();
    double gT_F(); // Total Condcutance Model (Potential)
    double Q_F();
    double VPD_F(double Tc);
    double Delta_F();
    double D_F();
    double Qn(double Tc);
    double La_F();
    double emiss_cloudy_F();
    double RaPot();
    double ea_F();
    double gH_Leaf(double Tc, double dT);
    double gV_Leaf(double Tc, double dT);
    double Kv_F();
    double gHa_Free(double Tc);
    double gHv_Free(double Tc);

}; // class CWSIModel
#endif // MYCWSIMODEL_

