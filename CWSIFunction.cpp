//@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
// 
// Equations to calculate theoretical CWSI, actual transpiration (Ta), and potential transpiration (Tp) 
// Note: models developed and calibrated for apple tree canopies
// 
// Written by Yasin Osroosh, Ph.D.
// Email: yosroosh@gmail.com
// 
//@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
//@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

#include <math.h>       /* pow */ /* exp */ /* abs */ 
#include "CWSIFunction.h"


//@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
//@ Crop Water Stress Index (Model)
//@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
// Function to calculate the crop water stress index (CWSI)
double CWSIModel::CWSI(double Tc) 
{
    double dT{ 0.0 };        //Tc-Ta for measured conditions (deg C)
    double dTu{ 0.0 };       //Tc-Ta of the upper limit (nontranspiring crop), hot, maximum stress (deg C)
    double dTl{ 0.0 };       //Tc-Ta of the lower limit (crop canopy transpiration not limited by availabel moisture), cold, nonstressed (deg C)
    double resultCWSI{ 0.0 };
    
    dT = Tc - data_.Ta;
    dTu = dTu_F();
    dTl = dTp_F();
    
    resultCWSI = (dT - dTl) / (dTu - dTl);
  
    if (dT <= dTl) {    	
        resultCWSI = 0.0001;
    }
    if (dT >= dTu) {    	
        resultCWSI = 1;
    }
    if (dTl == dTu) {
        resultCWSI = 0.0001;
    }
    
    if (resultCWSI < 0.0001) {
		resultCWSI = 0.0001;
	}
	
	return resultCWSI;
}

//@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
//@ Potential Transpiration
//@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
// Potential transpiration function
double CWSIModel::Tp_F()
{
    double ga{ 0.0 }; 
    double Qnc{ 0.0 }; 
    double emiss_cloudy{ 0.0 }; 
    double DeltaT{ 0.0 }; 
    double dT{ 0.0 }; 
    double T{ 0.0 }; 
    double d{ 0.0 }; 
    double VPD{ 0.0 }; 
    double resultTp_F{ 0.0 };
    double Tc{ 0.0 };

    d = D_F();
    VPD = data_.m * d + data_.b1;
    
    dT = dTp_F();
    Tc = dT + data_.Ta;
    
    emiss_cloudy = emiss_cloudy_F();

    ga = gH_Leaf(0, 0); // mol/m^2/s
    
    Qnc = Qn(Tc) -FE * ALPHA_LW * (SIGMA) * (pow ((Tc + 273), 4));
    
    resultTp_F = (Qnc - ga * CP_AIR * (dT)) * (0.018 * 24 * 3600) / Y; // mm/day
    
    //If Tp_F < 0 Then Tp_F = T_ACTUAL_MIN
    //If Tp_F < ETc_Min Then Tp_F = T_ACTUAL_MIN
    //If Tp_F > ET_Max Then Tp_F = T_Actual_Max
    
    return resultTp_F;
}

//@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
//@ Actual Transpiration
//@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
// Non-Linearized actual transpiration function
double CWSIModel::T_Actual_F(const double &Tc) 
{
    double ga{ 0.0 }; 
    double Qnc{ 0.0 }; 
    double DeltaT{ 0.0 }; 
    double dTu{ 0.0 }; 
    double resultT_Actual_F{ 0.0 };
    
    dTu = Tc - data_.Ta;
    ga = gH_Leaf(Tc, dTu); // mol/m^2/s
        
    Qnc = Qn(Tc) -FE * ALPHA_LW * (SIGMA) * (pow((Tc + 273), 4));  // Non-Linearized function
    resultT_Actual_F = (Qnc - CP_AIR * ga * dTu) * (0.018 * 24 * 3600) / Y; // mm/day
    
    if (resultT_Actual_F < 0) {
		resultT_Actual_F = T_ACTUAL_MIN;
	}
    if (resultT_Actual_F < T_ACTUAL_MIN) {
		resultT_Actual_F = T_ACTUAL_MIN;
	}
    //If T_Actual_F > ET_Max Then T_Actual_F = ET_Max
	
	return resultT_Actual_F; 	
}

//@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
//@ Sub Functions
//@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
double CWSIModel::Rnet(double Tc) 
{
    double ga{ 0.0 }; 
    double Qnc{ 0.0 }; 
    double emiss_cloudy{ 0.0 }; 
    double DeltaT{ 0.0 }; 
    double dTu{ 0.0 }; 
    double resultRnet{ 0.0 };
        
    resultRnet = Qn(Tc) -FE * ALPHA_LW * (SIGMA) * (pow((Tc + 273), 4.0));  // Non-Linearized function
    //If resultRnet < 0 Then resultRnet = 0
    
    return resultRnet;
}

double CWSIModel::dTu_F() 
{
    // gs = 0
    double ga{ 0.0 }; 
    double Q{ 0.0 }; 
    double emiss_cloudy{ 0.0 };
	double resultdTu_F{ 0.0 };

    emiss_cloudy = emiss_cloudy_F();
        
    ga = gH_Leaf(0, 0);  // mol/m^2/s
    Q = Q_F();
    
    resultdTu_F = (Q) / (ga * CP_AIR);
    
    return resultdTu_F;
}

double CWSIModel::dTp_F() 
{
    double Pa{ 0.0 }; 
    double VPD{ 0.0 };        // Vapor pressure deficit (kPa)
    double ga{ 0.0 };         // Aerodynamic conductance to heat
    double gv{ 0.0 };         // Aerodynamic conductance to vapor
    double gs{ 0.0 }; 
    double gT{ 0.0 }; 
    double dTu{ 0.0 }; 
    double es_Ta{ 0.0 };         // Saturation vapor pressure (kPa)
    double s{ 0.0 }; 
    double emiss_cloudy{ 0.0 }; 
    double d{ 0.0 }; 
    double dT{ 0.0 }; 
    double Q{ 0.0 }; 
    double n{ 0.0 }; 
    double resultdTp_F{ 0.0 };
    
    emiss_cloudy = emiss_cloudy_F();
    Pa = 101.3 * pow(((293 - 0.0065 * data_.Elev) / 293), 5.26);
        
    d = D_F();
    s = Delta_F() / Pa; // 1/C '
        
    ga = gH_Leaf(0, 0);  // mol/m^2/s
    gT = gT_F(); // mol/m^2/s
            
    n = (3 * ALPHA_LW - 4 *FE) * emiss_cloudy * SIGMA * (pow((data_.Ta + 273.15), 3));
    Q = Q_F();
    
    resultdTp_F = (Q - gT * Y * d / Pa) / (ga * CP_AIR - n + Y * gT * s);

    if (resultdTp_F < DTL_MIN) {
		resultdTp_F = DTL_MIN;
	}
	
	return resultdTp_F;
} 

// Total Condcutance Model (Potential)
double CWSIModel::gT_F() 
{
    double Pa{ 0.0 }; 
    double gT{ 0.0 }; 
    double gs{ 0.0 }; 
    double dTl{ 0.0 }; 
    double Tavg{ 0.0 }; 
    double d{ 0.0 }; 
    double Q{ 0.0 };         // Saturation vapor pressure (kPa)
    double s{ 0.0 }; 
    double Lc{ 0.0 }; 
    double la{ 0.0 }; 
    double emiss_cloudy{ 0.0 }; 
    double SC{ 0.0 }; 
    double VPD{ 0.0 }; 
    double gv{ 0.0 }; 
    double ga{ 0.0 }; 
    double Qnc{ 0.0 }; 
    double dT{ 0.0 }; 
    double T{ 0.0 };
    double n{ 0.0 }; 
    double resultgT_F{ 0.0 };
    double Tc{ 0.0 };
    
    emiss_cloudy = emiss_cloudy_F();
    d = D_F();
    VPD = data_.m * d + data_.b1;
    
    dT = ((data_.m - 1) * d + data_.b1) / Delta_F();
    Tc = dT + data_.Ta;

    emiss_cloudy = emiss_cloudy_F();
    Pa = 101.3 * pow(((293 - 0.0065 * data_.Elev) / 293), 5.26);
    
    // YO: check if Tc needs to be 0 or calculated
    ga = 0; //gH_Leaf(Tc, 0);  // mol/m^2/s
      
    n = (3 * ALPHA_LW - 4 *FE) * emiss_cloudy * SIGMA * (pow((data_.Ta + 273.15), 3));
    Q = Q_F();
    
    resultgT_F = ((Q + (n - ga * CP_AIR) * dT) / (Y * VPD / Pa));
        
    //Exit Function
    if (resultgT_F <= 0) {
        if (-0.5 < resultgT_F) {
            resultgT_F = GT_MIN;
        } else {
        	resultgT_F = GT_MAX;
		}
    } else {
    	resultgT_F = data_.b2 * resultgT_F + data_.b0;
	}
    
    if (resultgT_F > GT_MAX) {
		resultgT_F = GT_MAX;
	}	
	
	return resultgT_F;
}

double CWSIModel::Q_F() 
{
    // Ra: Total shortwave irradiance measured by pyranometer
    double Absorb{ 0.0 };  // Absorptivity
    double St1{ 0.0 };   // Transmitted radiation
    double la{ 0.0 }; 
    double resultQ_F{ 0.0 };
          
    Absorb = (1 - ALBEDO_LEAF - TRANS);
    la = La_F();
    St1 = TRANS * data_.Sgl; // Transmitted radiation
    
    resultQ_F = 0.25 * (Absorb * data_.Sgl + Absorb * St1 + 4 * (ALPHA_LW -FE) * la);
    
    return resultQ_F;
}

double CWSIModel::VPD_F(double Tc) 
{
    double es_Tc{ 0.0 };         // Saturation vapor pressure (kPa)
    double es_Ta{ 0.0 };         // Saturation vapor pressure (kPa)
    double ea_Tc{ 0.0 };         // Actual vapor pressure (kPa)
    double ea_Ta{ 0.0 };         // Actual vapor pressure (kPa)
    double resultVPD_F{ 0.0 };

    es_Tc = 0.6108 * exp(17.27 * Tc / (Tc + 237.3));
    es_Ta = 0.6108 * exp(17.27 * data_.Ta / (data_.Ta + 237.3));
    ea_Ta = es_Ta * data_.RH * 0.01;
    resultVPD_F = es_Tc - ea_Ta;
    
    return resultVPD_F;
}

double CWSIModel::Delta_F() {
	// The slope of the relationship between saturation vapour pressure (es, kPa) and air temperature (T, �C)
    // s (kPa �C-1), is given by (Tetens, 1930; Murray, 1967):
    double resultDelta_F{ 4098 * (0.6108 * exp(17.27 * data_.Ta / (data_.Ta + 237.3))) / (pow((data_.Ta + 237.3), 2)) };
    
    return resultDelta_F;
}

double CWSIModel::D_F() 
{
    double es_Tc{ 0.0 };         // Saturation vapor pressure (kPa)
    double es_Ta{ 0.0 };         // Saturation vapor pressure (kPa)
    double ea_Tc{ 0.0 };         // Actual vapor pressure (kPa)
    double ea_Ta{ 0.0 };         // Actual vapor pressure (kPa)
   
    es_Ta = 0.6108 * exp(17.27 * data_.Ta / (data_.Ta + 237.3));

    ea_Ta = es_Ta * data_.RH / 100;
   
    double resultD_F{ es_Ta - ea_Ta };
    
    return resultD_F;
} 

double CWSIModel::Qn(double Tc) 
{
	// Ra: Total shortwave irradiance measured by pyranometer
    double Absorb{ 0.0 };  // Absorptivity
    double Sr1{ 0.0 };  // Reflected radiation
    double Sr2{ 0.0 };  // Reflected radiation
    double St1{ 0.0 };   // Transmitted radiation
    double St2{ 0.0 };   // Transmitted radiation
    double St3{ 0.0 }; 
    double Srg{ 0.0 };   // Reflected radiation from ground
    double Stb{ 0.0 };  // Transmitted radiation through but leaf
    double Srt{ 0.0 };  // Reflected radiation from a shaded leaf surface (second time)
    double Lg{ 0.0 }; 
    double la{ 0.0 }; 
    double Lc{ 0.0 };  // emitted thermal radiation
    double Lc_trees{ 0.0 };  // emitted thermal radiation from other trees
    double Lc_trunk{ 0.0 };  // emitted thermal radiation from the trunk
    double Loe{ 0.0 }; 
    double G{ 0.0 }; 
    double An{ 0.0 }; 
    double Qn_Sunlit_Up{ 0.0 }; 
    double Qn_Shaded_Mid{ 0.0 }; 
    double Qn_Shaded_But{ 0.0 }; 
    double Qn_Shaded_Sid{ 0.0 }; 
    double Qn_Ground{ 0.0 }; 
    double emiss_cloudy{ 0.0 }; 
    double SW_Up{ 0.0 }; 
    double SW_Mid{ 0.0 }; 
    double SW_But{ 0.0 }; 
    double SW_Sid{ 0.0 }; 
    double LW_Up{ 0.0 }; 
    double LW_Mid{ 0.0 }; 
    double LW_But{ 0.0 }; 
    double LW_Sid{ 0.0 }; 
    
    emiss_cloudy = emiss_cloudy_F();
           
    Absorb = (1 - ALBEDO_LEAF - TRANS);
    
    la = La_F();
    // Linearize with 1st order Taylor�s Expansion Series:
    // Linearized form of "Lc = ALPHA_LW  * (SIGMA) * ((Tc + 273) ^ 4)"
    Lc = ALPHA_LW * (SIGMA) * (pow((Tc + 273), 4));
     
    // Shaded leaf
    St1 = TRANS * data_.Sgl; // Transmitted radiation
	
	Srg = ALBEDO_LEAF_G * data_.Sgl;
               
    // Leaves at the top of canopy are folded so they receive less direct radiation, however it is assumed
    // that view factor for the upper side of leaf would be 0.5
    SW_Up = Absorb * (0.5 * data_.Sgl); //+ 0.25 * Srg)
    LW_Up = ALPHA_LW * (0.5 * la + 0.5 * Lc);
    Qn_Sunlit_Up = SW_Up + LW_Up; // Qn = SW(in) - SW(out) + LW(in) - LW(out)
        
    SW_Mid = Absorb * (0.5 * St1);
    LW_Mid = ALPHA_LW * (2 * 0.5 * Lc);
    Qn_Shaded_Mid = SW_Mid + LW_Mid;

    double resultQn{ ((0.5 * Qn_Sunlit_Up + 0.5 * Qn_Shaded_Mid)) }; 
	
	return resultQn;       
}

double CWSIModel::La_F() 
{
    double resultLa_F{ emiss_cloudy_F() * (SIGMA) * (pow((data_.Ta + 273), 4)) };
    
    return resultLa_F;
}

double CWSIModel::emiss_cloudy_F() 
{
    double emiss_cleaRaky{ 0.0 }; 
    double emiss_cloudy{ 0.0 }; 
    double c{ 0.0 }; 
    double Ra_Max{ 0.0 }; 
    double Ea{ 0.0 }; 
    
    Ra_Max = RaPot(); // W/m^2
    c = (1 - (data_.Sgl / (Ra_Max * 1)));
    if (c < 0 ) {
    	c = 0;	
	} 
    //c = 1.5 * c;
    if (c > 1) {
		c = 1;
	}
    Ea = ea_F();
    emiss_cleaRaky = 1.72 * pow (Ea / (273.16 + data_.Ta), 0.1428571); // Page 163 & 164
    double resultemiss_cloudy_F{ ((1.00 - 0.84 * c) * emiss_cleaRaky + 0.84 * c) };
    
    return resultemiss_cloudy_F;
}

double CWSIModel::RaPot() 
{
	// Calculates Sgl, extraterrestrial incoming solar radiation, in MJ/m^2/day
    double SmallDelta{ 0.0 };     // solar declination  (Sgldians)
    double dr{ 0.0 };             // Inverse relative distance Earth-Sun
    double omega{ 0.0 };          // sunset hour angle  (Sgldians) 
	double resultRaPot{ 0.0 };   

    double lat_rad{ data_.lat * M_PI / 180 };    // Convert to radians
    SmallDelta = 0.409 * sin(2 * M_PI / 365.0 * data_.DOY - 1.39);
    omega = acos(-tan(lat_rad) * tan(SmallDelta));
    dr = 1 + 0.033 * cos(2 * M_PI * data_.DOY / 365.0);
    resultRaPot = 37.586 * dr * (omega * sin(lat_rad) * sin(SmallDelta) + cos(lat_rad) * cos(SmallDelta) * sin(omega));
    resultRaPot = resultRaPot * 11.574;  // 1 MJ/m^2/d =  1000000 J m-2 d-1 / 86400 s d-1 = 11.574 J m-2 s-1 = 11.574 W m-2
    
    return resultRaPot;
}

double CWSIModel::ea_F() 
{
	double es_Ta{ 0.0 };         // Saturation vapor pressure (kPa)
    double ea_Ta{ 0.0 };         // Actual vapor pressure (kPa)

    es_Ta = 0.6108 * exp(17.27 * data_.Ta / (data_.Ta + 237.3));
    ea_Ta = es_Ta * data_.RH * 0.01;
    double resultea_F{ ea_Ta };
    
    return resultea_F;
}

double CWSIModel::gH_Leaf(double Tc, double dT) 
{
	// dT = canopy and air temp difference (e.g. dTu)
    double gr{ 0.0 }; 
    double Re{ 0.0 }; 
    double Pr{ 0.0 }; 
    double gHa_Forced{ 0.0 }; 
    double Tavg{ 0.0 }; 
    double dTu{ 0.0 }; 
    double gH{ 0.0 }; 
    double g_Free{ 0.0 }; 
    double T{ 0.0 }; 
    double Kv{ 0.0 }; 

    Re = (data_.U * WIND_F * CD / Kv_F());
    Pr = Kv_F() / 2.14e-5;
    
    gH = 0.664 * 41.6 * 2.14e-5 * sqrt(Re) * pow(Pr, 0.33333) / CD;
        
    gHa_Forced = 1.4 * gH;   // mol/m2/s
        
    g_Free = gHa_Free(dT);  // mol/m^2/s
    
    // Factor 1.5 (1 for lower + 0.5 for upper) accounts for the fact that the upper side of leaf is folded cuasing
    // a decreased wind speed and therfore less sensible heat loss
    // Fctor 2 accounts for both sides of the leaf
    double resultgH_Leaf{ 2.0 * gHa_Forced }; // + gHa_Free)
    
    return resultgH_Leaf;
}

double CWSIModel::gV_Leaf(double Tc, double dT) 
{
    double gr{ 0.0 }; 
    double Re{ 0.0 }; 
    double Pr{ 0.0 }; 
    double gVa_Forced{ 0.0 }; 
    double gv{ 0.0 }; 
    double g_Free{ 0.0 }; 
    double T{ 0.0 }; 
    double Kv{ 0.0 }; 
    
    Re = (data_.U * WIND_F * CD / Kv_F());
    Pr = Kv_F() / (2.14 * pow(10, -5));
    
    gv = (0.664 * 41.6 * 2.14e-5 * (pow(Re, 0.5)) * (pow(Pr, 0.333333))) / CD;
       
    gVa_Forced = 1.4 * gv;  // mol/m2/s
    
    g_Free = gHv_Free(dT);  // mol/m^2/s

    double resultgV_Leaf{ 1 * (gVa_Forced) }; // + gHv_Free)
    
    return resultgV_Leaf;
}

double CWSIModel::Kv_F() {	
	double T{ 0.0 }; 
    T = data_.Ta + 273.16;
    double resultKv_F{ -1.1555e-14 * pow(T, 3) + 9.5728e-11 * pow(T, 2) + 3.7604e-8 * T - 3.4484e-6 };
    
    return resultKv_F;
} 

double CWSIModel::gHa_Free(double dT) 
{
    double gr{ 0.0 }; 
    double Re{ 0.0 }; 
    double Pr{ 0.0 }; 
    double gHa_Forced{ 0.0 }; 
    double Tavg{ 0.0 }; 
    double dTu{ 0.0 }; 
    double gH{ 0.0 }; 
    double T{ 0.0 }; 
    double resultgHa_Free{ 0.0 };
       
    if (dT < 0) {
        resultgHa_Free = 0;        
    } else {
    	Tavg = (data_.Ta + dT + data_.Ta) / 2.0;
    	gr = ((9.8 * (pow(CD, 3)) * (dT)) / ((Tavg + 273) * (pow(Kv_F(), 2))));
    	Re = (data_.U * cos(45) * CD / Kv_F());
    	Pr = Kv_F() / (2.14 * pow(10, -5));
    
   	 	resultgHa_Free = (0.54 * 41.6 * (2.14e-5)) * (pow(gr * Pr, 0.25)) / CD; // factor 0.75 (=(1+0.5)/2) is for the heated leaf facing down (average of one side factor 1 and the other side factor 0.5)
	}
	
	return resultgHa_Free;
}

double CWSIModel::gHv_Free(double dT) 
{
    double gr{ 0.0 }; 
    double Re{ 0.0 }; 
    double Pr{ 0.0 }; 
    double gHa_Forced{ 0.0 }; 
    double gHa_Free{ 0.0 }; 
    double Tavg{ 0.0 }; 
    double dTu{ 0.0 }; 
    double gH{ 0.0 }; 
    double T{ 0.0 }; 
    double resultgHv_Free{ 0.0 };
       
    if (dT < 0) {
    	resultgHv_Free = 0;
	} else {
		Tavg = (data_.Ta + dT + data_.Ta) / 2.0;
    	gr = ((9.8 * (pow(CD, 3)) * (dT)) / ((Tavg + 273) * (pow(Kv_F(), 2))));
    	Re = (data_.U * cos(45) * CD / Kv_F());
    	Pr = Kv_F() / (2.14 * pow(10, -5));
    
    	resultgHv_Free = (0.54 * 41.6 * (2.14e-5)) * (pow(gr * Pr, 0.25)) / CD; // factor 0.5 is for the heated leaf facing down
	}
	
	return resultgHv_Free;
}

