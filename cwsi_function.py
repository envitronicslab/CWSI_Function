import math

__author__ = "Y. Osroosh, Ph.D. <yosroosh@gmail.com>"

class ModelsInputData:
    """
    Models input data
    """
    def __init__(self):
        self.Ta = 25.5  # Air Temperature (deg C) -> measured
        self.RH = 60.6  # Relative Humidity (%) -> measured
        self.Sgl = 800.02  # Incoming short wave solar radiation (W/m^2) -> estimated
        self.U = 2.3  # Wind Speed (m/s) -> measured
        self.m = 1  # Calibration coeff.
        self.b1 = 0  # Calibration coeff.
        self.b2 = 3.0  # Calibration Coefficient -> default value = 8.00, 4.00 for Daylight
        self.b0 = 0.1  # Calibration Coefficient -> default value = 0.00
        self.lat = 42.3  # Latitude (degree)
        self.Elev = 2500.0  # Elevation (m)
        self.DOY = 180  # Day of Year (Julian Day)

class CWSIModel:
    """
        Class containing CWSI model definition to calculate theoretical crop water stress index, as well as 
        actual and potential transpiration models
    """
    # Model constants
    GT_MAX = 99  # gT_Max, mmol/m2/s
    GT_MIN = 0.0  # gT_Min, 0.00; mmol/m2/s
    DTL_MIN = -20  # dTl_Min, deg C
    ET_MAX = 99  # ET_Max, mm
    T_ACTUAL_MIN = 0.001  # T_Actual_Min, mm
    CD = 0.036  # Cd, = 0.72 * 0.05  # To Leaf Width = 5 cm
    WIND_F = 1.0  # Wind_f, 0.3 // WS(In-Canopy)/WS(Roza) at the middle of distance between the center of canopy and free air
    ALBEDO_LEAF = 0.09  # Albedo_leaf
    ALBEDO_LEAF_G = 0.1  # Albedo_leaf_g, Albedo of Ground = Thermal Absortptivity of Ground
    TRANS = 0.06  # Trans
    ALPHA_LW = 0.95  # alpha_lw, Longwave absorptivity of apple leaf
    CP_AIR = 29.17  # Cp_Air, Heat capacity of air J/mol/C
    Y = 44000  # y, J/mol
    SIGMA = 5.671e-8  # sigma, Stefan-Boltzmann constant (W/m^2/K^4)
    FE = 1.0  # Fe, 1.00; Fe is the view factor between the entire surface (only one side) of the leaf and the complete sphere of the view

    def __init__(self, data: ModelsInputData):
        self.data = data

    def CWSI(self, Tc):
        dT = Tc - self.data.Ta
        dTu = self.dTu_F()
        dTl = self.dTp_F()

        resultCWSI = (dT - dTl) / (dTu - dTl)

        if dT <= dTl:
            resultCWSI = 0.0001
        if dT >= dTu:
            resultCWSI = 1
        if dTl == dTu:
            resultCWSI = 0.0001

        if resultCWSI < 0.0001:
            resultCWSI = 0.0001

        return resultCWSI

    def Tp_F(self):   
        dT = self.dTp_F()
        Tc = dT + self.data.Ta        

        ga = self.gH_Leaf(0, 0)  # mol/m^2/s

        Qnc = self.Qn(Tc) - self.FE * self.ALPHA_LW * (self.SIGMA) * (math.pow((Tc + 273), 4))
        resultTp_F = (Qnc - ga * self.CP_AIR * (dT)) * (0.018 * 24 * 3600) / self.Y  # mm/day

        # if Tp_F < 0:
        #     Tp_F = self.T_ACTUAL_MIN
        # if Tp_F < self.ETc_Min:
        #     Tp_F = self.T_ACTUAL_MIN
        # if Tp_F > self.ET_MAX:
        #     Tp_F = self.T_Actual_Max

        return resultTp_F

    def T_Actual_F(self, Tc):
        dTu = Tc - self.data.Ta
        ga = self.gH_Leaf(Tc, dTu)  # mol/m^2/s

        Qnc = self.Qn(Tc) - self.FE * self.ALPHA_LW * (self.SIGMA) * (math.pow((Tc + 273), 4))  # Non-Linearized function
        resultT_Actual_F = (Qnc - self.CP_AIR * ga * dTu) * (0.018 * 24 * 3600) / self.Y  # mm/day

        if resultT_Actual_F < 0:
            resultT_Actual_F = self.T_ACTUAL_MIN
        if resultT_Actual_F < self.T_ACTUAL_MIN:
            resultT_Actual_F = self.T_ACTUAL_MIN

        # if resultT_Actual_F > self.ET_MAX:
        #     resultT_Actual_F = self.ET_MAX

        return resultT_Actual_F

    def Rnet(self, Tc):
        resultRnet = self.Qn(Tc) - self.FE * self.ALPHA_LW * (self.SIGMA) * (math.pow((Tc + 273), 4.0))  # Non-Linearized function
        # if resultRnet < 0:
        #     resultRnet = 0
        return resultRnet

    def dTu_F(self):
        ga = self.gH_Leaf(0, 0)  # mol/m^2/s
        Q = self.Q_F()

        resultdTu_F = (Q) / (ga * self.CP_AIR)

        return resultdTu_F

    def dTp_F(self):
        Pa = 101.3 * math.pow(((293 - 0.0065 * self.data.Elev) / 293), 5.26)

        d = self.D_F()
        s = self.Delta_F() / Pa  # 1/C '

        ga = self.gH_Leaf(0, 0)  # mol/m^2/s
        gT = self.gT_F()  # mol/m^2/s

        n = (3 * self.ALPHA_LW - 4 * self.FE) * self.emiss_cloudy_F() * self.SIGMA * (math.pow((self.data.Ta + 273.15), 3))
        Q = self.Q_F()

        resultdTp_F = (Q - gT * self.Y * d / Pa) / (ga * self.CP_AIR - n + self.Y * gT * s)

        if resultdTp_F < self.DTL_MIN:
            resultdTp_F = self.DTL_MIN

        return resultdTp_F

    def gT_F(self):
        d = self.D_F()
        VPD = self.data.m * d + self.data.b1

        dT = ((self.data.m - 1) * d + self.data.b1) / self.Delta_F()

        emiss_cloudy = self.emiss_cloudy_F()
        Pa = 101.3 * math.pow(((293 - 0.0065 * self.data.Elev) / 293), 5.26)

        # YO: check if Tc needs to be 0 or calculated
        ga = 0  # gH_Leaf(Tc, 0)  # mol/m^2/s

        n = (3 * self.ALPHA_LW - 4 * self.FE) * emiss_cloudy * self.SIGMA * (math.pow((self.data.Ta + 273.15), 3))
        Q = self.Q_F()

        resultgT_F = ((Q + (n - ga * self.CP_AIR) * dT) / (self.Y * VPD / Pa))

        # Exit Function
        if resultgT_F <= 0:
            if -0.5 < resultgT_F:
                resultgT_F = self.GT_MIN
            else:
                resultgT_F = self.GT_MAX
        else:
            resultgT_F = self.data.b2 * resultgT_F + self.data.b0

        if resultgT_F > self.GT_MAX:
            resultgT_F = self.GT_MAX

        return resultgT_F

    def Q_F(self):
        # Ra: Total shortwave irradiance measured by pyranometer
        Absorb = (1 - self.ALBEDO_LEAF - self.TRANS)
        la = self.La_F()
        St1 = self.TRANS * self.data.Sgl  # Transmitted radiation

        resultQ_F = 0.25 * (Absorb * self.data.Sgl + Absorb * St1 + 4 * (self.ALPHA_LW - self.FE) * la)

        return resultQ_F

    def VPD_F(self, Tc):
        es_Tc = 0.6108 * math.exp(17.27 * Tc / (Tc + 237.3))
        es_Ta = 0.6108 * math.exp(17.27 * self.data.Ta / (self.data.Ta + 237.3))
        ea_Ta = es_Ta * self.data.RH * 0.01
        resultVPD_F = es_Tc - ea_Ta

        return resultVPD_F

    def Delta_F(self):
        resultDelta_F = 4098 * (0.6108 * math.exp(17.27 * self.data.Ta / (self.data.Ta + 237.3))) / (math.pow((self.data.Ta + 237.3), 2))

        return resultDelta_F

    def D_F(self):
        es_Ta = 0.6108 * math.exp(17.27 * self.data.Ta / (self.data.Ta + 237.3))
        ea_Ta = es_Ta * self.data.RH / 100

        resultD_F = es_Ta - ea_Ta

        return resultD_F

    def Qn(self, Tc):
        # Initialize variables
        Absorb = 0.0
        St1 = 0.0
        la = 0.0
        Lc = 0.0
        G = 0.0
        Qn_Sunlit_Up = 0.0
        Qn_Shaded_Mid = 0.0
        SW_Up = 0.0
        SW_Mid = 0.0
        LW_Up = 0.0
        LW_Mid = 0.0

        # Calculate absorbance
        Absorb = 1 - self.ALBEDO_LEAF - self.TRANS

        # Calculate leaf temperature
        la = self.La_F()  # Assuming you have a function for this

        # Calculate emitted thermal radiation
        Lc = self.ALPHA_LW * self.SIGMA * (math.pow(Tc + 273, 4))

        # Shaded leaf calculations
        St1 = self.TRANS * self.data.Sgl  # Assuming data_.Sgl is define

        # Calculate Qn for sunlit and shaded leaves
        SW_Up = Absorb * (0.5 * self.data.Sgl)  # Adjust as needed
        LW_Up = self.ALPHA_LW * (0.5 * la + 0.5 * Lc)
        Qn_Sunlit_Up = SW_Up + LW_Up

        SW_Mid = Absorb * (0.5 * St1)
        LW_Mid = self.ALPHA_LW * (2 * 0.5 * Lc)
        Qn_Shaded_Mid = SW_Mid + LW_Mid

        resultQn = 0.5 * Qn_Sunlit_Up + 0.5 * Qn_Shaded_Mid

        return resultQn

    def La_F(self):
        resultLa_F = self.emiss_cloudy_F() * (self.SIGMA) * (math.pow((self.data.Ta + 273), 4))

        return resultLa_F

    def emiss_cloudy_F(self):
        emiss_cleaRaky = 0.0
        c = 0.0
        Ra_Max = self.RaPot()  # W/m^2
        Ea = self.ea_F()

        c = 1 - (self.data.Sgl / (Ra_Max * 1))
        if c < 0:
            c = 0
        if c > 1:
            c = 1

        emiss_cleaRaky = 1.72 * math.pow(Ea / (273.16 + self.data.Ta), 0.1428571)  # Page 163 & 164
        resultemiss_cloudy_F = (1.00 - 0.84 * c) * emiss_cleaRaky + 0.84 * c

        return resultemiss_cloudy_F

    def RaPot(self):
        # Calculates Sgl, extraterrestrial incoming solar radiation, in MJ/m^2/day
        SmallDelta = 0.0  # solar declination  (Sgldians)
        dr = 0.0  # Inverse relative distance Earth-Sun
        omega = 0.0  # sunset hour angle  (Sgldians)

        lat_rad = self.data.lat * math.pi / 180  # Convert to radians
        SmallDelta = 0.409 * math.sin(2 * math.pi / 365.0 * self.data.DOY - 1.39)
        omega = math.acos(-math.tan(lat_rad) * math.tan(SmallDelta))
        dr = 1 + 0.033 * math.cos(2 * math.pi * self.data.DOY / 365.0)
        resultRaPot = 37.586 * dr * (omega * math.sin(lat_rad) * math.sin(SmallDelta) + math.cos(lat_rad) * math.cos(SmallDelta) * math.sin(omega))
        resultRaPot = resultRaPot * 11.574  # 1 MJ/m^2/d =  1000000 J m-2 d-1 / 86400 s d-1 = 11.574 J m-2 s-1 = 11.574 W m-2

        return resultRaPot

    def ea_F(self):
        es_Ta = 0.6108 * math.exp(17.27 * self.data.Ta / (self.data.Ta + 237.3))
        ea_Ta = es_Ta * self.data.RH * 0.01
        resultea_F = ea_Ta

        return resultea_F

    def gH_Leaf(self, Tc, dT):
        Re = 0.0
        Pr = 0.0
        gHa_Forced = 0.0
        gH = 0.0

        Re = (self.data.U * self.WIND_F * self.CD / self.Kv_F())
        Pr = self.Kv_F() / 2.14e-5

        gH = 0.664 * 41.6 * 2.14e-5 * math.sqrt(Re) * math.pow(Pr, 0.33333) / self.CD

        gHa_Forced = 1.4 * gH  # mol/m2/s

        # Factor 1.5 (1 for lower + 0.5 for upper) accounts for the fact that the upper side of leaf is folded causing
        # a decreased wind speed and therefore less sensible heat loss
        # Factor 2 accounts for both sides of the leaf
        resultgH_Leaf = 2.0 * gHa_Forced 

        return resultgH_Leaf

    def gV_Leaf(self, Tc, dT):
        Re = 0.0
        Pr = 0.0
        gVa_Forced = 0.0
        gv = 0.0

        Re = (self.data.U * self.WIND_F * self.CD / self.Kv_F())
        Pr = self.Kv_F() / (2.14 * pow(10, -5))

        gv = (0.664 * 41.6 * 2.14e-5 * (math.pow(Re, 0.5)) * (math.pow(Pr, 0.333333))) / self.CD

        gVa_Forced = 1.4 * gv  # mol/m2/s

        resultgV_Leaf = 1 * (gVa_Forced) 

        return resultgV_Leaf

    def Kv_F(self):
        T = 0.0
        T = self.data.Ta + 273.16
        resultKv_F = -1.1555e-14 * math.pow(T, 3) + 9.5728e-11 * math.pow(T, 2) + 3.7604e-8 * T - 3.4484e-6

        return resultKv_F

    def gHa_Free(self, dT):
        gr = 0.0
        Pr = 0.0
        Tavg = 0.0
        resultgHa_Free = 0.0

        if dT < 0:
            resultgHa_Free = 0
        else:
            Tavg = (self.data.Ta + dT + self.data.Ta) / 2.0
            gr = ((9.8 * (math.pow(self.CD, 3)) * (dT)) / ((Tavg + 273) * (math.pow(self.Kv_F(), 2))))
            Pr = self.Kv_F() / (2.14 * pow(10, -5))

            resultgHa_Free = (0.54 * 41.6 * (2.14e-5)) * (math.pow(gr * Pr, 0.25)) / self.CD  # factor 0.75 (=(1+0.5)/2) is for the heated leaf facing down (average of one side factor 1 and the other side factor 0.5)

        return resultgHa_Free

    def gHv_Free(self, dT):
        gr = 0.0
        Pr = 0.0
        Tavg = 0.0
        resultgHv_Free = 0.0

        if dT < 0:
            resultgHv_Free = 0
        else:
            Tavg = (self.data.Ta + dT + self.data.Ta) / 2.0
            gr = ((9.8 * (math.pow(self.CD, 3)) * (dT)) / ((Tavg + 273) * (math.pow(self.Kv_F(), 2))))
            Pr = self.Kv_F() / (2.14 * pow(10, -5))

            resultgHv_Free = (0.54 * 41.6 * (2.14e-5)) * (math.pow(gr * Pr, 0.25)) / self.CD  # factor 0.5 is for the heated leaf facing down

        return resultgHv_Free