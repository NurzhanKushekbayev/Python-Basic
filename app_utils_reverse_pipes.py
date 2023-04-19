import math, sys
import os
from array import *
import numpy as np
import pandas as pd
import sys
# import xlwings as xw

# exl=xw.Book(r'SGtest.xlsx')

# filename = sys.argv[1]
filename = sys.path[0]+"/pipeline_reverse_calc_input_7.xlsx"

#b=xw.Book(r'pipelines_file_reverse_temp.xlsm')
# wb.activate

# f = open('output.txt', 'w')
# print("Distance (m)", "Pin (atm)", "Pout (atm)","Pressure gradient (atm/m)", "Vliq (m/s)", "Vgas (m/s)",sep="   ",file=f)


wb=pd.ExcelFile(filename, engine='openpyxl')
# calc_results=sys.argv[2]
# long_results=sys.argv[3]
calc_results=sys.path[0]+"/calc_results.csv"
long_results=sys.path[0]+"/long_results.csv"

res_arr = np.array([])
prop_arr = np.array([])
#flow_pattern=""


# *************************исходные данные начало*************************************************
# количество итераций расчета всей системы помимо первых двух (первые 2 выполняются всегда)
max_iter=3
# начальное допущение по неизвестным температурам в узлах сети
t_initial_guess=30

# ускорение свободного падения ft/sec2
g=32.174

# SGoil относительная плотность нефти
# SGgas плотность газа
# SGwater плотность воды
SGoil_default=0.860
SGgas_default=0.75
SGwater_default=1.03

# давление и температура при которых проведены измерения свойств (давление сепарации) (атм), (C)
Psep_m=1
Tsep_m=15.56

# # расход нефти м3/сут
# Qosc=5

# расход жидкости м3/сут
Qlsc_m=200


# Газовый фактор м3/м3
Rp_m=50


# обводненность % об.
WC=90
#print('WC',WC)

# точка инверсии
iversion_point=60

# начальные давление и температура (атм), (C)
P_m=5
T_m=50

################################### данные трубопровода #####################################
# внешний диаметр трубопровода (мм)
d_out_m=159

# толщина стенки (мм)
WT_m=6

# внутренний диаметр трубопровода (мм)
d_m=d_out_m-2*WT_m

# шероховатость в мм
roughness_m=0.018

# длина трубопровода
L_m=6000

# угол наклона
incl_angle_deg=0


# количество участков
increment=10



################################## тепловые свойства ####################################
# температурный градиент (геотермальный градиент) (F/ft)
gG=0

# температура окружающей среды (С)
Te_C=15

# глубина укладки трубопровода от поверхности земли до середины трубы (м)
BD_m=1

# *******сделать перевод из W/m*K
# коэф. теплопроводности нефти (Btu/hr*ft*F)
k_o=0.08
# коэф. теплопроводности газа (Btu/hr*ft*F)
k_g=0.02
# коэф. теплопроводности воды (Btu/hr*ft*F)
k_w=0.35
# коэф. теплопроводности стали (Btu/hr*ft*F)
k_s=25
# коэф. теплопроводности грунта может варьироваться в зависимости от типа и влажности (Btu/hr*ft*F)
k_e=1.4

# *************************исходные данные конец*************************************************



def units_conversion(Qlsc_m, Rp_m, Tsep_m, Psep_m, P_m, T_m, d_out_m, WT_m, roughness_m, incl_angle_deg, L_m, Te_C, BD_m):
    global Qlsc, Rp, Tsep, Psep, P, T, d, d_out, WT, roughness, incl_angle_rad, L, Te, BD, WC, d_m
    # **************************переводы единиц******************************************************
    # 
    # перевод расхода нефти в stb/d

    # если дано Q жидкости и обводненность
    if math.isnan(Qlsc_m)==True:
        Qlsc_m=0
    Qlsc= Qlsc_m / 0.158987

    if math.isnan(WC)==True:
        WC=0

    # Qosc=Qosc/0.158987


    # расход жидкости stb/d
    # Qlsc=Qosc/(1-WC/100)


    # Перевод ГФ в scf/STB
    if math.isnan(Rp_m)==True:
        Rp_m=0
    Rp=Rp_m*5.614583333

    # перевод температуры сепарации в град. F
    Tsep=Tsep_m*9/5+32
    #print('Tsep=',Tsep)

    # перевод давления сепарации в psi
    Psep=Psep_m*14.7
    #print ('Psep=',Psep)

    # перевод P и T в psi и F
    P=P_m*14.7
    # print(P_m)
    # print(P)
    T=T_m*9/5+32
    #print ('T=',T)
    #print ('P=',P)
    

    d_m=d_out_m-2*WT_m
    # перевод внутреннего диаметра в дюймы и в ft
    d=d_m/25.4/12


    # перевод внвнешнего диаметра в дюймы и в ft
    d_out=d_out_m/25.4/12

    # перевод толщины стенки в дюймы и в футы
    WT=WT_m/25.4/12

    # перевод в ft
    roughness= roughness_m*0.00328084

    # перевод угла в рад
    incl_angle_rad=incl_angle_deg*math.pi/180

    # перевод длины в ft
    L=L_m*3.28084

    # перевод температуры окружающей среды в (F)
    Te=Te_C*9/5+32

    # перевод глубины укладки трубопровода в ft 
    BD=BD_m*3.28084

    return Qlsc, Rp, Tsep, Psep, P, T, d, d_out, WT, roughness, incl_angle_rad, L, Te, BD
    

def oil_Rs (T, P, SGoil, SGgas, Rp, Tsep, Psep):
    # расчет газосодержания при P и T (scf/STBo)
    global y_api, y_g100, C1_VB, C2_VB, C3_VB, C1_VB2, C2_VB2, C3_VB2, C1_VB3, C2_VB3, C3_VB3, Rs

    # ***************************************************************************************************
    # расчет свойств нефти и газа

    # # поправка плотности нефти на температуру, если плотность дана при +20С приведение к 60 F (15.56 C)
    # dens_T_corr= -0.00131753501400560*SGoil + 0.00181909299719888
    # print(dens_T_corr)
    # # расчет относительной плотности нефти при 60 F (15.56 C)
    # SGoil_15=SGoil+(20-15.56)*dens_T_corr
    # print(SGoil_15)
    # yo=SGoil_15

    ###################################################
    # плотность в град. API
    y_api=141.5/SGoil-131.5
    #print('y_o API=', y_api)

    # корреляция Vasques & Beggs
    if y_api<=30:
        C1_VB=0.0362
        C2_VB=1.0937
        C3_VB=25.7245
        C1_VB2=4.677*10**(-4)
        C2_VB2=1.751*10**(-5)
        C3_VB2=-1.811*10**(-8)
        C1_VB3=27.62
        C2_VB3=0.914328
        C3_VB3=11.172
    else:
        C1_VB=0.0178
        C2_VB=1.1870
        C3_VB=23.9310
        C1_VB2=4.670*10**(-4)
        C2_VB2=1.100*10**(-5)
        C3_VB2=1.337*10**(-9)
        C1_VB3=56.18
        C2_VB3=0.84246
        C3_VB3=10.393

    # плотность газа, который выделился бы при давлении 100 psig
    y_g100=SGgas*(1+5.912*10**(-5)*y_api*Tsep*math.log10(Psep/114.7))
    #print ('y_g100=',y_g100)

    # газосодержание при P и T (scf/STBo)
    Rs=C1_VB*y_g100*P**C2_VB*math.exp(C3_VB*(y_api/(T+459.67)))
    if Rs>Rp:
        Rs=Rp
    #print('Rs=',Rs)
    #print('Rp=',Rp)
    return Rs


def oil_FVF(T, P, y_api, y_g100, Rs, C1_VB2, C2_VB2, C3_VB2, C1_VB3, C2_VB3, C3_VB3):
    global Bob, Pb, c_o, Bo
    # объемный коэф. нефти
    # давление насыщения psi
    
    Pb=(((C1_VB3*Rp)/(y_g100))*(10**(-1*C3_VB3*y_api/(T+459.67))))**C2_VB3
    #print('Pb=',Pb)

    # сжимаемость нефти
    a1_VB=-1433
    a2_VB=5
    a3_VB=17.2
    a4_VB=-1180
    a5_VB=12.61
    a6_VB=10**5
    
    c_o=(a1_VB+a2_VB*Rs+a3_VB*T+a4_VB*y_g100+a5_VB*y_api)/(a6_VB*P)

    # объемный коэффициент при давлении насыщения Bob по корреляции Standing bbl/STBO:
    Bob=0.9759+0.00012*(Rp*(SGgas/SGoil)**0.5+1.25*T)**1.2

    if P<=Pb:
        # для давлений ниже давлений насыщения объемный коэф. нефти Vazques and Beggs bbl/STBO:
        Bo=1+C1_VB2* Rs+(T-60)*(y_api/y_g100)*(C2_VB2+C3_VB2*Rs)
    else:
        # для давлений выше давлений насыщения объемный коэф. нефти Standing bbl/STBO:
        Bo=Bob*math.exp(-c_o*(P-Pb))
        

    #print('c_o=',c_o)
    #print('Bob=',Bob)
    #print('Bo=',Bo)
    return Pb, c_o, Bo, Bob


def gas_free_SG(SGgas, Rs, y_api):
    global y_gd, y_gt, y_gf
    # ***************************************************
    # плотность растворенного газа Katz
    a_katz=-0.000003549307588*y_api-0.000002031557018
    b_katz=0.019987951073597*y_api+0.251483883187716
    y_gd=Rs*a_katz+b_katz
    #print('y_gd=',y_gd)

    y_gt=SGgas
    # относительная плотность свободного газа
    # коммент от Beggs: все гидравлические расчеты должны проводиться с ипользованием y_gf
    if Rp==Rs:
        y_gf=0.56
    else:
        y_gf=(Rp*y_gt-Rs*y_gd)/(Rp-Rs)
        if y_gf<0.56:
            y_gf=0.56
    #print('y_gf=',y_gf)
    return y_gf


def gas_Z_crit_P_T (T, P, y_gf):
    # расчет критических параметров газа и фактора сжимаемости Z
    global Tpc, Ppc, Tpr, Ppr, Z
    # псевдокритическая температура и давление from Standing for Natural Gas Systems deg.R
    Tpc=168+325*y_gf-12.5*y_gf**2
    Ppc=677+15*y_gf-37.5*y_gf**2
    # псевдокритическая температура и давление from Standing for Gas Condensate Systems deg.R
    # T_pc=187+3330*y_gf-71.5*y_gf**2
    # P_pc=706-51.7*y_gf-11.1*y_gf**2

    # pseudoreduced T (deg R) & P
    Tpr=(T+459.67)/Tpc
    Ppr=P/Ppc

    # фактор сжимаемости газа Dranchuk & Abu-Kassem
    A1_DA = 0.3265
    A2_DA = -1.0700
    A3_DA = -0.5339
    A4_DA = 0.01569
    A5_DA = -0.05165
    A6_DA = 0.5475
    A7_DA = -0.7361
    A8_DA = 0.1844
    A9_DA = 0.1056
    A10_DA = 0.6134
    A11_DA = 0.7210


    # итерационный расчет Z
    Z_low = 0.1
    Z_high = 3
    i = 0
    
    while abs(Z_low - Z_high) > 0.000001:

        Z_mid = 0.5 * (Z_high + Z_low)
        # reduced gas density
        ro_r_low = 0.27 * Ppr / (Z_low * Tpr)
        Z_est_low = (1+(A1_DA + A2_DA / Tpr + A3_DA / Tpr ** 3 + A4_DA / Tpr ** 4 + A5_DA / Tpr ** 5) * ro_r_low +
                                (A6_DA + A7_DA / Tpr + A8_DA / Tpr ** 2) * ro_r_low ** 2 -
                                A9_DA * (A7_DA / Tpr + A8_DA / Tpr ** 2) * ro_r_low ** 5 +
                                A10_DA * (1 + A11_DA * ro_r_low ** 2) * ro_r_low ** 2 / Tpr ** 3 * math.exp(-A11_DA * ro_r_low ** 2))
        delta_low=Z_est_low-Z_low
        ro_r_mid = 0.27 * Ppr / (Z_mid * Tpr)
        Z_est_high = (1+(A1_DA + A2_DA / Tpr + A3_DA / Tpr ** 3 + A4_DA / Tpr ** 4 + A5_DA / Tpr ** 5) * ro_r_mid +
                                (A6_DA + A7_DA / Tpr + A8_DA / Tpr ** 2) * ro_r_mid ** 2 -
                                A9_DA * (A7_DA / Tpr + A8_DA / Tpr ** 2) * ro_r_mid ** 5 +
                                A10_DA * (1 + A11_DA * ro_r_mid ** 2) * ro_r_mid ** 2 / Tpr ** 3 * math.exp(-A11_DA * ro_r_mid ** 2))
        delta_high=Z_est_high-Z_mid

        if (delta_high * delta_low < 0):
            Z_high = Z_mid
        else:
            Z_low = Z_mid
        i = i + 1
    Z=Z_mid

    #print('Z=',Z)
    
    return Z


def gas_FVF_density (T, P, Z, y_gf):
    global Bg,ro_g
    # объемный коэффициент газа ft3/scf
    Bg=0.0283* Z *(T+459.67)/P
    #print('Bg=',Bg)
    # плотность свободного газа в lbm/ft3
    ro_g=2.7*y_gf*P/(Z*(T+459.67))
    #print('ro_g=',ro_g)
    return Bg,ro_g


def oil_density(P,SGoil,Rp,Rs,y_gt,Bob,y_gd, Pb):
    global ro_o
    # плотность нефти с растворенным газом lbm/ft3
    # при давлении насыщения
    ro_ob=(62.4*SGoil+0.0136*Rp*y_gt)/Bob

    if P<=Pb:
        # при давлении ниже давления насыщения
        ro_o=(62.4*SGoil+0.0136*Rs*y_gd)/Bo
    else:
        # при давлении выше давления насыщения
        ro_o=ro_ob*math.exp(c_o*(P-Pb))
    #print('ro_ob=',ro_ob)
    #print('ro_o=',ro_o)

    return ro_o


def oil_viscosity(T, y_api, Rs):
    global mu_od, mu_o
    # расчет вязкости дегазированной нефти
    if T<70:
        # если температура нефти меньше 70 F то используется экстраполяция Beal от вязкости при 70 и 80 С по Beggs & densitybinson (cP)
        mu_od70=10**(10**(3.0324-0.02023*y_api)/(70**1.163))-1
        mu_od80=10**(10**(3.0324-0.02023*y_api)/(80**1.163))-1
        Cbeal=math.log(mu_od70/mu_od80)/math.log(80/70)
        Bbeal=mu_od70*(70**Cbeal)
        log_mu_od=math.log(Bbeal)-Cbeal*math.log(T)
        mu_od=math.exp(log_mu_od)
    else:
        # вязкость дегазированной нефти по Beggs & Robinson (cP)
        mu_od=10**(10**(3.0324-0.02023*y_api)/(T**1.163))-1

        
    #print('mu_od=',mu_od)

    # вязкость газонасыщенной нефти по Beggs & Robinson (cP)
    mu_o=(10.715*(Rs+100)**(-0.515))*mu_od**(5.44*(Rs+150)**(-0.338))
    #print('mu_o=',mu_o)

    # # вязкость дегазированной нефти по Glaso
    # mu_od=(3.141*10**10)*T**(-3.444)*(math.log(y_api,10))**(10.313*math.log(T,10)-36.447)
    # print('mu_od=',mu_od)
    
    return mu_od, mu_o


def gas_viscosity(T, y_gf):
    global mu_g
    # вязкость газа по Lee et al. (cP)
    M_g=28.97*y_gf
    K_Lee=((9.4+0.02*M_g)*(T+459.67)**1.5)/(209+19*M_g+(T+459.67))
    X_Lee=3.5+986/(T+459.67)+0.01*M_g
    Y_Lee=2.4-0.2*X_Lee
    mu_g=10**(-4)*K_Lee*math.exp(X_Lee*(ro_g/62.4)**Y_Lee)
    #print ('mu_g=',mu_g)

    return mu_g


def oil_water_surface_tension (T, P, y_api):
    global surf_tens_o, surf_tens_w
    # поверхностное натяжение дегазированное нефти при 60 и 100F (dyne/cm):
    surf_tens_od_60= -0.257330977760622*y_api + 37.3333
    surf_tens_od_100= -0.260528160149153*y_api+ 39.144

    # поверхностное натяжение нефти с учетом растворенного газа (dyne/cm)
    surf_tens_o_60=(surf_tens_od_60*(0.000000000001586*P**4 - 0.000000012329818*P**3 +
                        0.000040595830174*P**2 - 0.082945237808161*P + 100.501979778122000)/100)
    surf_tens_o_100=(surf_tens_od_100*(0.000000000001586*P**4 - 0.000000012329818*P**3 +
                        0.000040595830174*P**2 - 0.082945237808161*P + 100.501979778122000)/100)
    #print('surf_tens_od_60=',surf_tens_od_60)
    #print('surf_tens_od_100=',surf_tens_od_100)
    #print('surf_tens_o_60=',surf_tens_o_60)
    #print('surf_tens_o_100=',surf_tens_o_100)
    if T>=100:
        surf_tens_o=surf_tens_o_100
    else:
        surf_tens_o=surf_tens_o_60
    
    # поверхностное натяжение воды в dyne/cm
    surf_tens_w_280= (-0.00000000000002029688*P**4 + 0.00000000040837191741*P**3 -
                        0.00000194094202273055*P**2 - 0.00388426751703546*P + 51.8826441393669)
    surf_tens_w_74= (0.00000000000001404278*P**4 - 0.00000000034107006638*P**3 +
                        0.00000299215712911354*P**2 - 0.01228475469300550000*P + 73.6922284104043)

    if T==74:
        surf_tens_w=surf_tens_w_74
    else:
        if T==280:
            surf_tens_w=surf_tens_w_280
        else:
            surf_tens_w=(T-74)*(surf_tens_w_280-surf_tens_w_74)/(280-74)+surf_tens_w_74
    
    return surf_tens_o, surf_tens_w


def water_props(T, P, SGwater):
    global Bw, ro_w, mu_w, Rsw

    # расчет свойств воды
    # объемный коэффициент bbl/STBW
    Tx=T-60
    # объемный коэф. bbl/STBW
    Bw=1+1.2*10**(-4)*Tx+1*10**(-6)*Tx**2-3.33*10**(-6)*P
    #print('Bw=',Bw)
    # плотность воды lbm/ft3
    ro_w=62.4*SGwater/Bw
    #print('ro_w=',ro_w)
    # вязкость воды cP
    mu_w=math.exp(1.003-1.479*10**(-2)*T+1.982*10**(-5)*T**2)
    #print('mu_w=',mu_w)
    # газосодержание scf/STBW
    Rsw=(2.12+((3.45*10**(-3))*T-(3.59*10**(-5))*T**2) + (0.0107-(5.26*10**(-5))*T+(1.48*10**(-7))*T**2)*P+
                (-(8.75*10**(-7))+(3.9*10**(-9))*T-(1.02*10**(-11))*T**2)*P**2)
    if Rp==0:
        Rsw=0
    #print('Rsw=',Rsw)

    return Bw, ro_w, mu_w, Rsw


def oil_liq_gas_balance(Qlsc, WC, Rp, Rs, Rsw, Bo, Bg, Bw, ro_o, ro_w, ro_g):
    global Qosc, Qgsc, Qwsc, q_o, q_l, q_w, q_g, f_o, f_w, Xm_o, Xm_g, mass_rate_m, lambda_l
    
    # расчет фазового состояния при P и T

    if math.isnan(WC)==True:
        WC=0
        Qosc=0
    else:
        Qosc=Qlsc*(1-WC/100)
    Qgsc=Rp*Qosc
    #print('Qosc',Qosc)
    #print('Qlsc',Qlsc)
    #print('Qgsc',Qgsc)
    #print('Qosc=',Qosc)
    # расход нефти в раб. условиях bbl/d
    Qo=Qosc* Bo
    #print('Qo=',Qo)
    # расход воды в ст. условиях bbl/d
    Qwsc=Qlsc-Qosc
    #print('Qwsc=',Qwsc)
    # расход воды в рабочих условиях bbl/d
    Qw=Qwsc*Bw
    #print('Qw=',Qw)
    # количество газа в ст. усл. scf/d
    Qgsc=Qosc*Rp
    #print('Qgsc=',Qgsc)
    # свободный газ при рабочих условиях ft3/d
    Qg=(Qgsc - Qosc*Rs - Qwsc* Rsw)* Bg
    # свободный газ при стандартных условиях scf/d
    Qg2=(Qgsc - Qosc*Rs-Qwsc*Rsw)
    # количество газа растворенного в нефти ft3/d
    qgdo=Qosc*Rs*Bg
    # количество газа растворенного в воде ft3/d
    qgdw=Qwsc*Rsw*Bg
    total_diss_gas_sc=(qgdo+qgdw)/Bg
    if Qosc==0 and Qwsc==0:
        diss_glr=0
    else:
        diss_glr=total_diss_gas_sc/(Qosc+Qwsc)

    if Qg<0:
        Qg=0
        Qg2=0
        qgdw=0
        
    

    # расчет расходов по фазам
    Ql=Qo+Qw

    # q_w, q_o, q_l в ft3/sec
    q_w=Qw*5.614/86400
    q_o=Qo*5.614/86400
    q_l=Ql*5.614/86400

    # q_g в ft3/sec
    q_g=Qg/86400

    # доля нефти f_o
    if q_o==0 and q_w==0:
        f_o=0
    else:
        f_o=q_o/(q_o+q_w)

    #print('f_o',f_o)

    # доля воды f_w
    f_w=1-f_o
    #print('f_w',f_w)

    # массовый расход смеси lbm/sec
    mass_rate_m=q_g*ro_g+q_o*ro_o+q_w*ro_w

    # массовая доля нефти в жидкости
    if q_o==0 and q_w==0:
        Xm_o=0
    else:
        Xm_o=q_o*ro_o/(q_o*ro_o+q_w*ro_w)
    # массовая доля газа в смеси
    if q_o==0 and q_w==0 and q_g==0:
        Xm_g=0
    else:
        Xm_g=q_g*ro_g/(q_o*ro_o+q_w*ro_w+q_g*ro_g)

    # доля трубопровода занятая жидкостью
    if q_l==0 and q_g==0:
        lambda_l=1
    else:
        lambda_l=q_l/(q_l+q_g)
    #print('lambda_l',lambda_l)

    return Qosc, Qgsc, Qwsc, q_o, q_l, q_w, q_g, f_o, f_w, Xm_o, Xm_g, mass_rate_m, lambda_l, WC


def thermal_props (T,SGoil, y_api, y_gf):
    global Cp_o,Cp_w, Cp_g

    # тепловые свойства
    # теплоемкость нефти
    # Cpo=(0.388+0.00045*T)/(y_api**0.5)
    #print('Cpo=',Cpo)
    # теплоемкость нефти по Gambill (Btu/lbm*F):

    # # Вариант 0 - по Gambill в чистом виде
    # Cp_o0=(0.388+0.00045*T)/((SGoil)**0.5)

    # # # Вариант 1
    # # *********поправка на температуру по Standing
    # # температурный коэффициент нефти по Standing
    # betta_o=0.000442+0.0000103*y_api
    # # температурный коэффициент воды по Standing
    betta_w=0.0002115+1.32*10**(-6)*((T-32)*5/9)+1.09*10**(-8)*((T-32)*5/9)**2
    # # относительная плотность с поправкой на температуру по Standing
    # SGoil_t1=SGoil*(1+((T-32)*5/9-20)*betta_w)/(1+((T-32)*5/9-20)*betta_o)
    # Cp_o1=(0.388+0.00045*T)/((SGoil_t1)**0.5)

    # Вариант 2
    # *********поправка на температуру воды по Standing
    # а в расчете относительной плотности нефти (SGoil_t2) используется
    # плотность газонасыщенной нефти ro_o

    # поправка плотности чистой воды на T (998.2 - плотность воды при 20С)
    ro_wt=998.2/(1+((T-32)*5/9-20)*betta_w)
    # перевод плотности в lbm/ft3
    ro_wt=ro_wt/16.018
    # относительная плотность газонасыщенной нефти при P и T
    SGoil_t2=ro_o/ro_wt
    Cp_o2=(0.388+0.00045*T)/((SGoil_t2)**0.5)
    Cp_o=Cp_o2
    #print('Cp_o=',Cp_o)

    # теплоемкость воды слабо зависит от Т, принята постоянной (Btu/lbm*F)
    Cp_w=1.001

  
    # # теплоемкость газа по Somerton (1992) !!!выдает неадекватно высокое значение 
    
    # A_Cp_g=-2.031*10**(-7)*(P*0.00689476)**4+5.702*10**(-5)*(P*0.00689476)**3-0.005175*(P*0.00689476)**2+0.1315*(P*0.00689476)-0.1813
    # B_Cp_g=2.84*10**(-7)*(P*0.00689476)**4-7.88*10**(-5)*(P*0.00689476)**3+0.00714*(P*0.00689476)**2-0.188*(P*0.00689476)+0.207
    # C_Cp_g=-1.42*10**(-6)*(P*0.00689476)**4+3.88*10**(-4)*(P*0.00689476)**3-0.03513*(P*0.00689476)**2+0.980*(P*0.00689476)-0.872
    # D_Cp_g=2.95*10**(-7)*(P*0.00689476)**4-7.97*10**(-5)*(P*0.00689476)**3+0.00725*(P*0.00689476)**2-0.222*(P*0.00689476)+0.533
    # E_Cp_g=-2.00*10**(-7)*(P*0.00689476)**4+5.46*10**(-5)*(P*0.00689476)**3-0.00519*(P*0.00689476)**2+0.193*(P*0.00689476)+1.928
    # a=A_Cp_g*((T-32)*5/9)**4
    # b=B_Cp_g*((T-32)*5/9)**3
    # c=C_Cp_g*((T-32)*5/9)**2
    # d=D_Cp_g*((T-32)*5/9)
    # Cp_g=A_Cp_g/(((T-32)*5/9)**4)+B_Cp_g*((T-32)*5/9)**3+C_Cp_g*((T-32)*5/9)**2+D_Cp_g*((T-32)*5/9)+E_Cp_g

    # теплоемкость газа по "Correlation for Natural Gas Heat Capacity" Developed by Dr. Moshfeghian
    a_cp=1.15
    b_cp=1.008
    c_cp=-0.944
    d_cp=0.533
    e_cp=1.110
    f_cp=0.0216
    Cp_g=(a_cp*b_cp**T*T**c_cp+d_cp*e_cp**(P/1000)*(P/1000)**f_cp)*(y_gf/0.6)**0.025
    #print('Cp_g=',Cp_g)


    
    return Cp_o, Cp_w, Cp_g


def mixture_props (WC, ro_o, ro_w, ro_g, f_o, f_w, surf_tens_o, surf_tens_w, Cp_o, Cp_w, Cp_g, Xm_o, Xm_g, lambda_l):
    # свойства смеси
    global ro_l, mu_l, ro_n, mu_n, surf_tens_liq, Cp_l, Cp_m, k_l, k_f

    # плотность жидкости в lbm/ft3 в рабочих условиях
    ro_l=ro_o*f_o+ro_w*f_w
   
    # вязкость жидкости (cP) set visc. of mixture to continious phase visc. (WC=60 inversion point)
    if math.isnan(WC)==True:
        WC=0
    if WC>=iversion_point:
        mu_l=mu_w
    else:
        mu_l=mu_o
        
    # вариант с расчетом вязкости в зависимости от доли нефти и воды в жидкости
    # mu_l=mu_o*f_o+mu_w*f_w

    #print('mu_l=',mu_l)

    # поверхностное натяжение жидкости (dyne/cm)
    surf_tens_liq=surf_tens_o*f_o+surf_tens_w*f_w
    #print('surf_tens_liq=',surf_tens_liq)
    
    # теплоемкость жидкости
    
    Cp_l=Xm_o*Cp_o+(1-Xm_o)*Cp_w
    #print('Cp_l=',Cp_l)

    # теплоемкость смеси жидкости и газа (Btu/lbm*F)
        
    Cp_m=Cp_g*Xm_g+(1-Xm_g)*Cp_l
    #print('Cp_m=',Cp_m)

    # ****************************************************
    # расчет свойств без проскальзывания (non-slip)
    # non-slip вязкость смеси газа и жидкости (cP)
    mu_n=mu_l*lambda_l+mu_g*(1-lambda_l)

    # non-slip плотность смеси газа и жидкости (lbm/ft3)
    ro_n=ro_l*lambda_l+ro_g*(1-lambda_l)
    #print('ro_n',ro_n)
    #print('mu_n',mu_n)


    # коэф. теплопроводности смеси нефти и воды
    k_l=k_o*(Xm_o)+k_w*(1-Xm_o)

    # коэф. теплопроводности смеси нефти и газа
    k_f=k_l*(1-Xm_g)+k_g*(Xm_g)

    return ro_l, mu_l, ro_n, mu_n, surf_tens_liq, Cp_l, Cp_m, k_l, k_f


def non_slip_velosities (d, q_l, q_g, surf_tens_liq, ro_l, g):
    global v_sl, v_sl_ms, v_sg, v_sg_ms, v_m, v_m_ms, Nre, Nfr, Nlv

    # расчет скоростей (superficial velocities)

    # площадь сечения трубопровода в ft2
    Ap=	math.pi/4*(d)**2
    #print('Ap',Ap)

    #print('q_l',q_l)
    # скорость жидкости v_sl (ft/sec)
    v_sl=q_l/Ap
    v_sl_ms=v_sl*0.3048
    #print('v_sl',v_sl)


    # скорость газа v_sg (ft/sec)
    
    v_sg=q_g/Ap
    v_sg_ms=v_sg*0.3048
    #print('q_g',q_g)
    #print('v_sg',v_sg)

    # суммарная скорость смеси (ft/sec)
    v_m=v_sl+v_sg
    v_m_ms=v_m*0.3048
    #print('v_m',v_m)
    # Критерий Рейнольдса для вязкости в cP (1 cp = 1488 lb/ft*s)
    Nre=1488* ro_n * v_m *d/mu_n
    #print('Nre=', Nre)
    
    # расчет числа Фруда
    Nfr=v_m**2/(g * d)
    #print('Nfr',Nfr)

    # критерий скорости по Hagedorn & Brown
    Nlv=1.938*v_sl*(ro_l/(surf_tens_liq))**(1/4)
    
    return v_sl, v_sl_ms, v_sg, v_sg_ms, v_m, v_m_ms, Nre, Nfr, Nlv


def flow_pattern_prediction (lambda_l, Nfr):
    global flow_pattern, fl_ptrn_num, L1, L2, L3, L4
    
    # определение границ структур потока

    L1=316*lambda_l**0.302
    L2=0.000925*lambda_l**(-2.468)
    L3=0.1*lambda_l**(-1.452)
    L4=0.5*lambda_l**(-6.738)

    #print('L1',L1)
    #print('L2',L2)
    #print('L3',L3)
    #print('L4',L4)

    # ***************************************************
    # определение структуры потока по условиям
    if lambda_l<0.01 and Nfr<L1 or lambda_l>=0.01 and Nfr<L2:
        flow_pattern="Segregated"
        fl_ptrn_num=1
    else:
        if lambda_l>=0.01 and L3>=Nfr>=L2:
            flow_pattern="Transition"
            fl_ptrn_num=4
        else:
            if 0.4>lambda_l>=0.01 and L3<Nfr<=L1 or lambda_l>=0.4 and L3<Nfr<=L4:
                flow_pattern="Intermittent"
                fl_ptrn_num=2
            else:
                if 0.4>lambda_l and Nfr>=L1 or lambda_l>=0.4 and Nfr>L4:
                    flow_pattern="Distributed"
                    fl_ptrn_num=3

    #print('Flow Pattern: ', flow_pattern)
    #print('Flow Pattern Number: ', fl_ptrn_num)

    return flow_pattern, fl_ptrn_num, L1, L2, L3, L4


def holdup_prediction (fl_ptrn_num, Nfr, Nlv, lambda_l, incl_angle_deg, incl_angle_rad, ro_l, v_sl, surf_tens_liq,L1, L2, L3, L4):
    global HL_0, HL_incl

    # определение holdup
    if fl_ptrn_num==1 or fl_ptrn_num==2 or fl_ptrn_num==3:
        # a, b и c для segregated
        if fl_ptrn_num==1:
            A_HL=0.980
            B_HL=0.4846
            C_HL=0.0868
        else:
            # a, b и c для intermittent
            if fl_ptrn_num==2:
                A_HL=0.845
                B_HL=0.5351
                C_HL=0.0173
            else:
                # a, b и c для distributed
                if fl_ptrn_num==3:
                    A_HL=1.065
                    B_HL=0.5824
                    C_HL=0.0609
        
        # HL для горизонтальной трубы
        if Nfr==0:
            HL_0=1
        else:
            HL_0=A_HL*lambda_l**B_HL/(Nfr**C_HL)
        
        if HL_0<lambda_l:
            HL_0=lambda_l
        if HL_0>1:
            HL_0=0.924


        # расчет поправок на угол уклона трубопровода
        if fl_ptrn_num==3 and incl_angle_deg>0 or incl_angle_deg==0:
            C_PSI=0
            PSI=1
        else:
            if fl_ptrn_num==1 and incl_angle_deg>0:
                E_C=0.011
                F_C=-3.7680
                G_C=3.5390
                H_C=-1.6140
            else:
                if fl_ptrn_num==2 and incl_angle_deg>0:
                    E_C=2.960
                    F_C=0.3050
                    G_C=-0.4473
                    H_C=0.0978
                else:
                    if fl_ptrn_num==1 or fl_ptrn_num==2 or fl_ptrn_num==3 and incl_angle_deg<0:
                        E_C=4.700
                        F_C=-0.3692
                        G_C=0.1244
                        H_C=-0.5056

            C_PSI=(1-lambda_l)*math.log(E_C*(lambda_l**F_C)*(Nlv**G_C)*(Nfr**H_C))
        
        if C_PSI<0:
            C_PSI=0
        # поправка на угол
        PSI=1+C_PSI*(math.sin(1.8*incl_angle_rad)-0.333*(math.sin(1.8*incl_angle_rad))**3)
        

        #print('HL_0=',HL_0)
        HL_incl=HL_0*PSI
        if HL_incl>1:
            HL_incl=0.924


    else:
        # интерполяция для Transition flow pattern

        if fl_ptrn_num==4:
            # a, b и c для segragated
            A_HL= 0.980
            B_HL= 0.4846
            C_HL= 0.0868
            HL_0_seg=(A_HL * lambda_l ** B_HL)/(Nfr**C_HL)
            if HL_0_seg<lambda_l:
                HL_0_seg=lambda_l
            if HL_0_seg>1:
                HL_0_seg=0.924
            
            # a, b и c для intermittent
            A_HL=0.845
            B_HL=0.5351
            C_HL=0.0173
            HL_0_int=A_HL*lambda_l**B_HL/(Nfr**C_HL)
            if HL_0_int<lambda_l:
                HL_0_int=lambda_l
            if HL_0_int>1:
                HL_0_int=0.924
            
            if incl_angle_deg==0:
                C_PSI_seg=0
                PSI_seg=1
                C_PSI_int=0
                PSI_int=1

            
            else:
                if incl_angle_deg<0:
                    # поправка на угол для segregated & intermittent downhill
                    E_C=4.700
                    F_C=-0.3692
                    G_C=0.1244
                    H_C=-0.5056
                    C_PSI_seg=(1-lambda_l)*math.log(E_C*(lambda_l**F_C)*(Nlv**G_C)*(Nfr**H_C))
                    if C_PSI_seg<0:
                        C_PSI_seg=0
                    PSI_seg=1+C_PSI_seg*(math.sin(1.8*incl_angle_rad)-0.333*(math.sin(1.8*incl_angle_rad))**3)
                    PSI_int=PSI_seg
                    
                
                else:
                    # segregated uphill
                    E_C=0.011
                    F_C=-3.7680
                    G_C=3.5390
                    H_C=-1.6140
                    C_PSI_seg=(1-lambda_l)*math.log(E_C*(lambda_l**F_C)*(Nlv**G_C)*(Nfr**H_C))
                    if C_PSI_seg<0:
                        C_PSI_seg=0
                    # поправка на угол
                    PSI_seg=1+C_PSI_seg*(math.sin(1.8*incl_angle_rad)-0.333*(math.sin(1.8*incl_angle_rad))**3)

                    # intermittent uphill
                    E_C=2.960
                    F_C=0.3050
                    G_C=-0.4473
                    H_C=0.0978
                    C_PSI_int=(1-lambda_l)*math.log(E_C*(lambda_l**F_C)*(Nlv**G_C)*(Nfr**H_C))
                    if C_PSI_int<0:
                        C_PSI_int=0
                    # поправка на угол
                    PSI_int=1+C_PSI_int*(math.sin(1.8*incl_angle_rad)-0.333*(math.sin(1.8*incl_angle_rad))**3)
            
            HL_incl_seg= HL_0_seg * PSI_seg
            if HL_incl_seg>1:
                HL_incl_seg=0.924
            HL_incl_int= HL_0_int * PSI_int
            if HL_incl_int>1:
                HL_incl_int=0.924
            
            A_tr=(L3-Nfr)/(L3-L2)
            #print ('A_tr=',A_tr)
            HL_incl_tr=A_tr*HL_incl_seg+(1-A_tr)*HL_incl_int
            HL_0_tr=A_tr*HL_0_seg+(1-A_tr)*HL_0_int
            HL_0=HL_0_tr
            #print('HL_0_seg=',HL_0_seg)
            #print('HL_0_int=',HL_0_int)
            #print('HL_0=',HL_0_tr)
            HL_incl=HL_incl_tr
            if HL_incl<lambda_l:
                HL_incl=lambda_l
            
    # корректировка Payne:
    if incl_angle_deg>=0:
        HL_incl=0.924*HL_incl
    else:
        if incl_angle_deg<0:
            HL_incl=0.685*HL_incl
    if HL_incl<lambda_l:
        HL_incl=lambda_l

    #print('Holdup',HL_incl)

    return Nlv, HL_0, HL_incl


def slip_velocities (ro_l, ro_g, mu_l, mu_g, HL_incl):
    global ro_s, mu_s, v_l, v_l_ms, v_g, v_g_ms, Nre_s
    # расчет свойств с проскальзыванием (slip)
    # slip вязкость смеси газа и жидкости (cP)
    mu_s=mu_l*HL_incl+mu_g*(1-HL_incl)
    
    # плотность смеси с учетом проскальзывания (lbm/ft3) (slip)
    ro_s=ro_l * HL_incl + ro_g * (1-HL_incl)
    #print('ro_s=',ro_s)

    # скорость жидкости с проскальзыванием ft/sec
    v_l=v_sl/(HL_incl)
    v_l_ms=v_l*0.3048
    #print('v_l_ms=', v_l_ms)


    # скорость газа с проскальзыванием ft/sec
    if (1-HL_incl)==0:
        v_g=0
    else:
        v_g=v_sg/(1-HL_incl)

    # скорость газа в м/с
    v_g_ms=v_g*0.3048

    #print('v_g_ms=', v_g_ms)

    Nre_s=1488* ro_s * (v_l+v_g) *d/mu_s

    return ro_s, mu_s, v_l, v_l_ms, v_g, v_g_ms, Nre_s


def friction_factor (d, roughness, Nre, lambda_l, HL_incl):
    global rel_roughness, f_n, frict_factor, f_f_ratio
    # относительная шероховатость
    rel_roughness=roughness/d
    #print('roughness=', roughness)
    #print('rel_roughness=', rel_roughness)



    if 0<Nre<2000:
        f_n=64/Nre
    
    else:
        if Nre==0:
            f_n=0
        else:
            # расчет нормализующего фактора трения normalizing friction factor (f_n)
            f_low = 0.0001
            f_high = 30
            i = 0

            while abs(f_low - f_high) > 0.00001:
                
                f_mid = 0.5 * (f_high + f_low)
                
                
                f_est_low = (1.74-2*math.log10(2*rel_roughness+18.7/(Nre *(f_low**0.5))))**(-2)
                delta_low=f_est_low-f_low
                
                f_est_high = (1.74-2*math.log10(2*rel_roughness+18.7/(Nre *(f_mid**0.5))))**(-2)
                delta_high=f_est_high-f_mid

                if (delta_high * delta_low < 0):
                    f_high = f_mid
                else:
                    f_low = f_mid
                i = i + 1
            #print('i=',i)
            f_n=f_mid
  
        
    #print('f_n=',f_n)

    # расчет фактора трения

    y=lambda_l/(HL_incl**2)
    #print('y=',y)
    if y==1:
        s=0
    else:
        if 1<y<1.2:
            s=math.log(2.2*y-1.2)
        else:
            s=math.log(y)/(-0.0523+3.182*math.log(y)-0.8725*(math.log(y))**2+0.01853*(math.log(y))**4)
    #print('s=',s)
    frict_factor=f_n*math.exp(s)
    if f_n==0:
        f_f_ratio=0
    else:
        f_f_ratio=frict_factor/f_n
    #print('frict_factor=',frict_factor)

    return rel_roughness, f_n, frict_factor, f_f_ratio


def pressure_gradient (P, d, frict_factor, incl_angle_rad, v_m, v_sg, ro_n, ro_s, g):
    global E_k, pressure_grad, pressure_grad_psift, pressure_grad_atmm
    # расчет градиента давления
  

    # E_k - безразмерная кинетическая энергия
    # *** Brill не использует в расчетах /(1-E_k) отмечено, что E_k>1 некорректно считается по формуле:
    E_k=v_m*v_sg*ro_s/P/g
    if E_k>=1 or E_k>0.01 or E_k<0.0:
        E_k=0
    # E_k=0
    #print('E_k=',E_k)

    
    # градиент давления (psf/ft) *** Brill не использует в расчетах "/(1-E_k)". отмечено, что E_k>1 некорректно считается по формуле
    pressure_grad=((frict_factor * ro_n * v_m**2)/(2 * d * g) + ro_s * g * math.sin(incl_angle_rad)/g)/(1-E_k)
    #print('Pressure Gradient (psf/ft)=',pressure_grad)

    # градиент давления (psi/ft)
    pressure_grad_psift=pressure_grad/144
    #print('Pressure Gradient (psi/ft)=',pressure_grad_psift)

    # перевод градиента давления в атм/м
    pressure_grad_atmm=pressure_grad_psift/14.7*3.28084

    return E_k, pressure_grad, pressure_grad_psift, pressure_grad_atmm


def pressure_loss (P, incr_L, pressure_grad_psift, direction):
    global delta_P, delta_P_atm, Pin, Pin_atm, distance_m

    # расчет перепада давления на расстоянии incr_L метров

    # перепад давления в psi
    delta_P=incr_L*pressure_grad_psift
    delta_P_atm=delta_P/14.7
    #print('delta_P (atm)', delta_P_atm)

    # давление на входе
    Pin=P-direction*delta_P
    Pin_atm=Pin/14.7
    if Pin_atm<1:
        Pin_atm=1
        Pin=14.7
    distance_m=distance/3.28084
    return delta_P, delta_P_atm, Pin, Pin_atm, distance_m


def JT_coefficient (T, q_l, Cp_g, ro_g):
    global JT_mix, dz_dt, Z_T
    # расчет коэффициента Джоуля-Томсона
    # расчет частной производной dz/dT
    delta_T=0.5
    T_dt=T+delta_T
    Z_T=Z

    Rs_T=Rs
    Rs_T_dt = oil_Rs (T_dt, P, SGoil, SGgas, Rp, Tsep, Psep)
    y_gf_T_dt=gas_free_SG (SGgas,Rs_T_dt, y_api)
    Z_dt=gas_Z_crit_P_T (T_dt, P, y_gf_T_dt)

    dZ=Z_dt-Z_T
    dz_dt=dZ/(delta_T)
 
    # расчет коэффициента Джоуля-Томсона для газа JT_gas (F/psi).
    # 
    JT_gas=1/(Cp_g)*(T/(Z_T* ro_g)*dZ/delta_T)
    #print('JT_gas=',JT_gas)

    # расчет коэффициента Джоуля-Томсона для смеси JT_mix (F/psi)
    if q_g==0 and q_l==0:
        JT_mix=0
    else:
        JT_mix=-1/(Cp_m* mass_rate_m)*(q_g *(-T/Z_T*dZ/delta_T)+ q_l)
    #print('JT_mix=',JT_mix)
    #print('Z=',Z)
    #print('Z_dt=',Z_dt)
    # ************************расчет коэффициента Джоуля-Томсона закончен
    # #############################################################################
    return JT_mix


def fl_temp_est (T,Te,BD,d,d_out,distance,incl_angle_rad,mass_rate_m,Cp_m,mu_n,ro_n, k_f, k_e, v_m, Nre, frict_factor, incr_L, g, JT_mix,pressure_grad_psift,gG):
    # расчет ТЕМПЕРАТУРЫ по выкидным линиям обратным ходом
    global Npr, h_f_WmC, h_f, h_soil_WmC, h_soil, OHTC_w_m_K, OHTC, T_f_C, T_f, T_C, Nnu

    # T=T_out_fl_C*9/5+32
    # расчет внутренней конвекции:
    # критерий Прандтля
    Npr=Cp_m*mu_n/k_f
    # Npr=Cp_m*mu_n*0.000671968975*3600/k_f

    
    # критерий нуссельта
    if Nre>10000:
    # Dittus and Boelter proposed the following dimensionless correlation for fully turbulent flow of single-phase fluids:
        if T>=Te:
            n_coef_pr=0.3
        else:
            n_coef_pr=0.4
        # критерий Нуссельта
        Nnu=0.0255*Nre**0.8*Npr**n_coef_pr
        
    else:
        # If the flow is laminar (i.e., NRe < 2100), h_f may be calculated using Hausen’s equation
        if Nre<=2100:
            # проверить корректно ли использовать инкремент (в оригинале Lo - where, Lo is the distance from the pipe inlet to the point of interest)
            if Npr<5:
                Nnu=1.86*(Nre*Npr/(distance/d))**(1/3)
                if Nnu<3.66:
                    Nnu=3.66
            else:
                Nnu=3.66+(0.0668*(d/distance)*Nre*Npr)/(1+0.04*((d/distance)*Nre*Npr)**(2/3))
        else:
            if 2100<Nre<10000:
                # For the transition region (2100 < Rei < 10**4) A correlation proposed by Gnielinski may be used to calculate h_f in this region
                # либо использовать f_n по диаграмме Moody (в оригинале)
                Nnu=((frict_factor/8)*(Nre-1000)*Npr)/(1+((12.7*(frict_factor/8)**0.5)*(Npr**(2/3)-1)))
    
    h_f=Nnu*k_f/d
    # h_f=4.36*k_f/d
    h_f_WmC= h_f*5.678263398

    # расчет U (Btu/hr*ft2*F) - OHTC - overall heat transfer coefficient

    if BD==d_out/2:
        BD=BD*1.01

    # глубина укладки от поверхности до верхней образующей трубы (b_h)
    b_h=BD-d_out/2
    
    if b_h>=d_out:
        # коэффициент теплопередачи грунту для случая когда глубина укладки от поверхности до верхней образующей трубы (b_h) больше, чем d_out
        h_soil=k_e/(math.log(2*BD/d_out+((2*BD/d_out)**2-1)**0.5))
        
    else:
        # коэффициент теплопередачи грунту для случая когда глубина укладки от поверхности до верхней образующей трубы (b_h) меньше, чем d_out
        h_soil=k_e/((d_out/2)*((math.acosh(2*BD/d_out))))

    h_soil_WmC=h_soil*5.678263398
    
    
    # a_hsoil=-0.0000236797*((T-32)*5/9)**3 + 0.0029866553*((T-32)*5/9)**2 - 0.0985424088*((T-32)*5/9) + 3.0189251975
    a_hsoil= 0.0037662338*((T-32)*5/9) + 2.0207792208
    # a_hsoil=1
    
    # расчет (rto*U)**-1
    rto_OHTC_minus_1 = 1/(d/2*h_f)+math.log((d_out/2)/(d/2))/k_s + a_hsoil/(h_soil)
       
    OHTC=1/(rto_OHTC_minus_1*(d_out/2))
    OHTC_w_m_K=OHTC*5.678263398
    #print('OHTC_w_m_K=',OHTC_w_m_K)


    # Расчет длины релаксации A_rel (ft)
    A_rel=Cp_m * (mass_rate_m*60*60) / (OHTC*math.pi*d)

    # Расчет безразмерного параметра Ф (Fi)
    # mechanical equivalent of heat J (foot pound force (4.1550 J·cal−1))
    J=778.24
    
    # dv/dl
    # if incr==1:
    #     dv=0
    # else:
    #     dv=v_m-v_prev
    dv=0
    dL=incr_L
    
    # в Brill&Mukherjee перед ro_n в самом начале формулы стоит домножение на J, у Alves нет этого множителя
    if pressure_grad_psift==0:
        Fi=0
    else:
        Fi=(J*ro_n * JT_mix * Cp_m * pressure_grad_psift - ro_n * 32.174/g * math.sin(incl_angle_rad) - ro_n* v_m /g*dv/dL)/pressure_grad_psift
   
    #print(Fi)
    if A_rel==0:
        T_f=T
    else:
        a_1=(Te-gG* dL *math.sin(incl_angle_rad))
        a_2=(T-Te)*math.exp(-dL/A_rel)
        a_3=gG*math.sin(incl_angle_rad)*A_rel*(1-math.exp(-dL/A_rel))
        a_4=1/(J*ro_n*Cp_m)*pressure_grad_psift*Fi*A_rel*(1-math.exp(-dL/A_rel))

        a_1=(Te-gG* dL *math.sin(incl_angle_rad))
        a_2=(T-Te)*math.exp(dL/A_rel)
        a_3=gG*math.sin(incl_angle_rad)*A_rel*(1-math.exp(dL/A_rel))
        a_4=1/(J*ro_n*Cp_m)*pressure_grad_psift*Fi*A_rel*(1-math.exp(dL/A_rel))
        T_f=a_1+a_2+a_3+a_4
        if math.isnan(T_f):
            T_f=30*9/5+32
        if T_f>140:
            T_f=30*9/5+32
        # TTTTT=(TTTTT-32)*5/9
        # T_f= (Te-gG* dL *math.sin(incl_angle_rad)) + (T-Te)*math.exp(-dL/A_rel)+gG*math.sin(incl_angle_rad)*A_rel*(1-math.exp(-dL/A_rel))+1/(J*ro_n*Cp_m)*pressure_grad_psift*Fi*A_rel*(1-math.exp(-dL/A_rel))
    
    T_C=(T-32)*5/9
    T_f_C=(T_f-32)*5/9
    T_delta=T_C-T_f_C
    
    # print('начальная температура',T_C)
    # print('конечная температура',T_f_C)

    return Npr, h_f_WmC, h_f, h_soil_WmC, h_soil, OHTC_w_m_K, OHTC, T_f_C, T_f, T_C


def temp_estimation (T,Te,BD,d,d_out,distance,incl_angle_rad,mass_rate_m,Cp_m,mu_n,ro_n, k_f, k_e, v_m, Nre, frict_factor, incr_L, g, JT_mix,pressure_grad_psift,gG):
    # расчет ТЕМПЕРАТУРЫ
    global Npr, h_f_WmC, h_f, h_soil_WmC, h_soil, OHTC_w_m_K, OHTC, T_f_C, T_f, T_C, Nnu


    # расчет внутренней конвекции:
    # критерий Прандтля
    Npr=Cp_m*mu_n/k_f
    # Npr=Cp_m*mu_n*0.000671968975*3600/k_f

    
    # критерий нуссельта
    if Nre>10000:
    # Dittus and Boelter proposed the following dimensionless correlation for fully turbulent flow of single-phase fluids:
        if T>=Te:
            n_coef_pr=0.3
        else:
            n_coef_pr=0.4
        # критерий Нуссельта
        Nnu=0.0255*Nre**0.8*Npr**n_coef_pr
        
    else:
        # If the flow is laminar (i.e., NRe < 2100), h_f may be calculated using Hausen’s equation
        if Nre<=2100:
            # проверить корректно ли использовать инкремент (в оригинале Lo - where, Lo is the distance from the pipe inlet to the point of interest)
            if Npr<5:
                Nnu=1.86*(Nre*Npr/(distance/d))**(1/3)
                if Nnu<3.66:
                    Nnu=3.66
            else:
                Nnu=3.66+(0.0668*(d/distance)*Nre*Npr)/(1+0.04*((d/distance)*Nre*Npr)**(2/3))
        else:
            if 2100<Nre<10000:
                # For the transition region (2100 < Rei < 10**4) A correlation proposed by Gnielinski may be used to calculate h_f in this region
                # либо использовать f_n по диаграмме Moody (в оригинале)
                Nnu=((frict_factor/8)*(Nre-1000)*Npr)/(1+((12.7*(frict_factor/8)**0.5)*(Npr**(2/3)-1)))
    
    h_f=Nnu*k_f/d
    # h_f=4.36*k_f/d
    h_f_WmC= h_f*5.678263398

    # расчет U (Btu/hr*ft2*F) - OHTC - overall heat transfer coefficient

    if BD==d_out/2:
        BD=BD*1.01

    # глубина укладки от поверхности до верхней образующей трубы (b_h)
    b_h=BD-d_out/2
    
    if b_h>=d_out:
        # коэффициент теплопередачи грунту для случая когда глубина укладки от поверхности до верхней образующей трубы (b_h) больше, чем d_out
        h_soil=k_e/(math.log(2*BD/d_out+((2*BD/d_out)**2-1)**0.5))
        
    else:
        # коэффициент теплопередачи грунту для случая когда глубина укладки от поверхности до верхней образующей трубы (b_h) меньше, чем d_out
        h_soil=k_e/((d_out/2)*((math.acosh(2*BD/d_out))))

    h_soil_WmC=h_soil*5.678263398
    
    
    # a_hsoil=-0.0000236797*((T-32)*5/9)**3 + 0.0029866553*((T-32)*5/9)**2 - 0.0985424088*((T-32)*5/9) + 3.0189251975
    a_hsoil= 0.0037662338*((T-32)*5/9) + 2.0207792208
    # a_hsoil=1
    
    # расчет (rto*U)**-1
    rto_OHTC_minus_1 = 1/(d/2*h_f)+math.log((d_out/2)/(d/2))/k_s + a_hsoil/(h_soil)
       
    OHTC=1/(rto_OHTC_minus_1*(d_out/2))
    OHTC_w_m_K=OHTC*5.678263398
    #print('OHTC_w_m_K=',OHTC_w_m_K)


    # Расчет длины релаксации A_rel (ft)
    A_rel=Cp_m * (mass_rate_m*60*60) / (OHTC*math.pi*d)

    # Расчет безразмерного параметра Ф (Fi)
    # mechanical equivalent of heat J (foot pound force (4.1550 J·cal−1))
    J=778.24
    
    # dv/dl
    # if incr==1:
    #     dv=0
    # else:
    #     dv=v_m-v_prev
    dv=0
    dL=incr_L
    
    # в Brill&Mukherjee перед ro_n в самом начале формулы стоит домножение на J, у Alves нет этого множителя
    if pressure_grad_psift==0:
        Fi=0
    else:
        Fi=(J*ro_n * JT_mix * Cp_m * pressure_grad_psift - ro_n * 32.174/g * math.sin(incl_angle_rad) - ro_n* v_m /g*dv/dL)/pressure_grad_psift
   
    #print(Fi)
    if A_rel==0:
        T_f=T
    else:
        T_f= (Te-gG* dL *math.sin(incl_angle_rad)) + (T-Te)*math.exp(-dL/A_rel)+gG*math.sin(incl_angle_rad)*A_rel*(1-math.exp(-dL/A_rel))+1/(J*ro_n*Cp_m)*pressure_grad_psift*Fi*A_rel*(1-math.exp(-dL/A_rel))
    if math.isnan(T_f):
        T_f=Te
    # if T_f<140:
    #     T_f=140
    T_C=(T-32)*5/9
    T_f_C=(T_f-32)*5/9
    # print(T_f_C)
    T_delta=T_C-T_f_C
    
    #print('начальная температура',T_C)
    #print('конечная температура',T_f_C)

    return Npr, h_f_WmC, h_f, h_soil_WmC, h_soil, OHTC_w_m_K, OHTC, T_f_C, T_f, T_C


def bottlenecks(ro_n, v_m_ms,Qlsc_m,mu_o,ro_o,WC,surf_tens_w,surf_tens_o,ro_w,mu_w,iversion_point,q_l,q_g):
    global EV, Erosion_flag, comment
    # эрозионная скорость по api rp 14e (m/s) c=80 как у Шелл с тв частицами
    EV=80/(ro_n**0.5)/3.281
    if v_m_ms>=EV:
        Erosion_flag="Эрозия"
    else:
        Erosion_flag="Нет эрозии"
    comment=Erosion_flag
    # определение скорости при которой произойдет расслоение эмульсии -> риск ручейковой коррозии (после точки энверсии при WC>iversion_point риск коррозии имеет место быть при любой скорости)
    # if WC<iversion_point:
        # D_em=((2.32*(Qlsc_m/24/60/60)**0.44)*((mu_o *0.001/(ro_o*16.01846))**0.27)*((ro_o*16.01846)**0.084)*((iversion_point/100-WC/100)**0.21)*((1-WC/100)**0.17))/((9.81**0.24)*((surf_tens_o*0.001)**0.11)*((ro_w*16.01846-ro_o*16.01846)**0.11))*(((3*mu_w*0.001/(ro_w*16.01846)+2*mu_o*0.001/(ro_o*16.01846))/(3*mu_w*0.001/(ro_w*16.01846)+3*mu_o*0.001/(ro_o*16.01846)))**0.35)
        # V_em=0.07*((surf_tens_w-surf_tens_o)*0.001)**0.56*(ro_w*16.01846-ro_o*16.01846)**0.24*(9.81**0.24)/((d_m/1000)**0.21*(ro_o*16.01846)**0.19*(mu_o)**0.61*(1-WC/100)**0.38)*((iversion_point/100-WC/100)/(iversion_point/100))**(-0.48)*(((3*mu_w+3*mu_o)/(3*mu_w+2*mu_o))**0.81)
        # V_em2=(Qlsc_m/24/60/60)/(math.pi*(D_em*0.5)**2/4)
        # V_em=(1.56*(d_m/1000)**0.59*(((surf_tens_w-surf_tens_o)*0.001)**0.44)*(9.81**0.38)*(ro_w*16.01846-ro_o*16.01846)**0.38)/((mu_o)**0.65*(1-WC/100)**0.18)*(((3*mu_w+3*mu_o)/(3*mu_w+2*mu_o))**0.82)*((WC/100**0.38*(1-lambda_l)**1.53*(1-HL_incl)**0.38*(1-(1-HL_incl)*(1+0.1*(1-lambda_l)**(-0.4)))**0.38)/(lambda_l**0.3*((1-lambda_l)-(1-HL_incl))**0.38*(1+((1-lambda_l)-(1-HL_incl))/(lambda_l))**0.15))
        # if V_em>v_m:
        #     comment=comment+'; '+'возможно расслоение эмульсии -> высокий риск ручейковой коррозии'
        # else:
        #     comment=comment+'; '+'скорость потока выше скорости эмульгирования -> низкий риск ручейковой коррозии'

    # else:
    #     V_em=0
    #     comment=comment+'; '+'прямая эмульсия (нефть в воде) -> риск коррозии'
    

    # solids_rate=0
    # solids_velocity=solids_rate/(math.pi*(d_m/1000)**2/4)

    v_rec_ms=0.9
    v_rec=v_rec_ms/0.3048
    if v_m_ms<v_rec_ms:
        # # площадь сечения трубопровода в ft2
        # Ap=	math.pi/4*(d)**2
        # # скорость жидкости v_sl (ft/sec)
        # v_sl=q_l/Ap
        # v_sl_ms=v_sl*0.3048
        # # скорость газа v_sg (ft/sec)
    
        # v_sg=q_g/Ap
        # v_sg_ms=v_sg*0.3048
    
        # ft
        d_rec=((q_l+q_g)/(v_rec*math.pi/4))**0.5
        d_rec=d_rec*304.8
        d_out_set = [89,114,159,219,273,325,426,525,730,820,1020]
        d_rec=min(d_out_set, key=lambda x:abs(x-d_rec))

        comment=comment+'; '+'низкая скорость потока (менее '+str(v_rec_ms)+' м/с)'
        if d_rec!=d_out_m:
            comment=comment+', '+'рекомендуемый диаметр '+str(d_rec)+' мм'
        



# считывание данных с листа pipelines_sort или pipelines
#if wb.sheets["pipelines_sort"].range(('B2')).value != None:
#if wb.parse(wb.sheet_names.index("pipelines_sort")).empty!=True:
#    pipe_list=wb.parse(wb.sheet_names.index("pipelines_sort"))
#    pipe_df=pd.DataFrame(pipe_list)
#    pipes_count=len(pipe_df.index)-1
    

pipe_list=wb.parse(wb.sheet_names==0)
pipe_df=pd.DataFrame(pipe_list)
pipes_count=len(pipe_df.index)-1
    
    #  сортировка dataframe - коллекторы перемещаются в начало
collector_df=pd.DataFrame(columns=pipe_df.columns.values)
i=0
j=0
pipe_df2=pipe_df.copy()
for row in pipe_df.iterrows():
      if math.isnan(pipe_df.iloc[i,pipe_df.columns.get_loc('qliq')])==True:
          collector_df.loc[j] = pipe_df.iloc[i]
          pipe_df2=pipe_df2.drop([i])
          j+=1
      i+=1
i=0
j=0
    
    # collector_df = collector_df.replace({pd.np.nan: None})
    # collector_df =collector_df.where(pd.notnull(collector_df), None)
collector_df4=collector_df.copy()
    
coll_sort_df2=pd.DataFrame(columns=collector_df.columns.values)
    
    # coll_sort_df2.columns = pipe_df2.iloc[0]
for row in collector_df.iterrows():
    # print(collector_df)  
    # if i<=max(collector_df.index.values):    - не нужно, проблема была в том, что после этого цикла строки 1472 и далее были смещены (tab->) и выполнялись в этом цикле.
    if math.isnan(collector_df.iloc[i,collector_df.columns.get_loc('press_end')])==True:
    # if collector_df.iloc[i,collector_df.columns.get_loc('press_end')]==None:  
        coll_sort_df2.loc[j] = collector_df.iloc[i]
        collector_df4=collector_df4.drop([i])
        j+=1
    i+=1
collector_df4=collector_df4.append(coll_sort_df2,ignore_index=True)

collector_df=collector_df4.copy()

# collector_df.columns = pipe_df2.iloc[0]


i=0
j=0
k=0

coll_sort_df=pd.DataFrame(columns=collector_df.columns.values)

collector_df2=collector_df
collector_df3=collector_df2
change=1
#print(collector_df)
while change>0:
    drop_list=array('i')
    change=0
    # a=0
    i=0
    j=0
    k=0
    for row in collector_df.iterrows():
        #print("COOLLECTOR DF: \n",collector_df)
        col_start=collector_df.iloc[i,collector_df.columns.get_loc('start_point')]
        # print(col_start)
        j=0
        for row2 in collector_df2.iterrows():
            if collector_df2.iloc[j,collector_df2.columns.get_loc('end_point')]==col_start and j<i:

                empty_row= pd.DataFrame([], index=[1])
                collector_df3 = pd.concat([collector_df2.iloc[:j], empty_row, collector_df2.iloc[j:]])
                collector_df3.iloc[j]=collector_df.iloc[i]
                collector_df3 = collector_df3.reset_index(drop=True)
                collector_df3=collector_df3.drop([collector_df.index[i]])
                collector_df3 = collector_df3.reset_index(drop=True)
                change=2
                k+=1
                break
            j+=1
        collector_df2=collector_df3
        collector_df=collector_df3
        if change==2:
            break
        i+=1
    
# collector_df3.set_index('№', inplace=True)
# collector_df3.columns = pipe_df2.iloc[0]


#pipe_df2.columns = pipe_df2.iloc[0]
#pipe_df2=pipe_df2.drop([0])


pipe_df=collector_df3.append(pipe_df2,ignore_index=True)

pipe_df3=pipe_df
# print(pipe_df)
#wb.sheets["pipelines_sort"].range((2,2)).value=pipe_df3
#pipe_list=wb.sheets["pipelines_sort"].range((2,2)).expand().value
# pipe_df=pd.DataFrame(pipe_list)

# заполнение пропусков по плотностям (для труб подключенных к источникам), если плотность не задана берем default SG
i=0
for row in pipe_df.iterrows():
    if math.isnan(pipe_df.iloc[i,pipe_df.columns.get_loc('qliq')])!=True:
        if math.isnan(pipe_df.iloc[i,pipe_df.columns.get_loc('SGoil')])==True:
            pipe_df.iloc[i,pipe_df.columns.get_loc('SGoil')]=SGoil_default
        if math.isnan(pipe_df.iloc[i,pipe_df.columns.get_loc('SGgas')])==True:
            pipe_df.iloc[i,pipe_df.columns.get_loc('SGgas')]=SGgas_default
        if math.isnan(pipe_df.iloc[i,pipe_df.columns.get_loc('SGwater')])==True:
            pipe_df.iloc[i,pipe_df.columns.get_loc('SGwater')]=SGwater_default
    i+=1

# определение Qж WC и ГФ для коллекторов
i=0
j=0
# цикл выполняется пока хотя бы одно из значений расхода = None
while pipe_df.iloc[:,pipe_df.columns.get_loc('qliq')].isnull().sum()>0:
    i=0
    for row in pipe_df.iterrows():
        if math.isnan(pipe_df.iloc[i,pipe_df.columns.get_loc('qliq')])==True:
            col_start=pipe_df.iloc[i,pipe_df.columns.get_loc('start_point')]
            j=0
            Qlsc_m=0
            Qosc_m=0
            Qgsc_m=0
            for row2 in pipe_df.iterrows():
                if col_start==pipe_df.iloc[j,pipe_df.columns.get_loc('end_point')] and math.isnan(pipe_df.iloc[j,pipe_df.columns.get_loc('qliq')])!=True:
                    
                    # объемные расходы по фазам
                    Qlsc_m=Qlsc_m+pipe_df.iloc[j,pipe_df.columns.get_loc('qliq')]
                    Qosc_m=Qosc_m+pipe_df.iloc[j,pipe_df.columns.get_loc('qliq')]*(1-pipe_df.iloc[j,pipe_df.columns.get_loc('wc')]/100)
                    Qgsc_m=Qgsc_m+pipe_df.iloc[j,pipe_df.columns.get_loc('qliq')]*(1-pipe_df.iloc[j,pipe_df.columns.get_loc('wc')]/100)*pipe_df.iloc[j,pipe_df.columns.get_loc('gazf')]
                    
                    pipe_df.iloc[i,pipe_df.columns.get_loc('qliq')]=Qlsc_m
                    
                    if Qlsc_m!=0:
                        pipe_df.iloc[i,pipe_df.columns.get_loc('wc')]=(Qlsc_m-Qosc_m)/Qlsc_m*100
                        pipe_df.iloc[i,pipe_df.columns.get_loc('gazf')]=Qgsc_m/Qosc_m
                    
                j+=1
                if Qlsc_m==0:
                    pipe_df.iloc[i,pipe_df.columns.get_loc('qliq')]=np.nan
        i+=1

# расчет плотностей для коллекторов
i=0
j=0
# цикл выполняется пока хотя бы одно из значений плотности = Nan

while pipe_df.iloc[:,pipe_df.columns.get_loc('SGoil')].isnull().sum()>0:
    i=0
    for row in pipe_df.iterrows():
        if math.isnan(pipe_df.iloc[i,pipe_df.columns.get_loc('SGoil')])==True:
            col_start=pipe_df.iloc[i,pipe_df.columns.get_loc('start_point')]
            j=0
            Oil_mass_rate=0
            Gas_mass_rate=0
            Water_mass_rate=0
            SGoil_sum=0
            SGgas_sum=0
            SGwater_sum=0

            for row2 in pipe_df.iterrows():
                if col_start==pipe_df.iloc[j,pipe_df.columns.get_loc('end_point')] and math.isnan(pipe_df.iloc[j,pipe_df.columns.get_loc('SGoil')])!=True:
                    

                    Oil_mass_rate=Oil_mass_rate+pipe_df.iloc[j,pipe_df.columns.get_loc('SGoil')]*pipe_df.iloc[j,pipe_df.columns.get_loc('qliq')]*(1-pipe_df.iloc[j,pipe_df.columns.get_loc('wc')]/100)
                    SGoil_sum=SGoil_sum+pipe_df.iloc[j,pipe_df.columns.get_loc('SGoil')]*(pipe_df.iloc[j,pipe_df.columns.get_loc('SGoil')]*pipe_df.iloc[j,pipe_df.columns.get_loc('qliq')]*(1-pipe_df.iloc[j,pipe_df.columns.get_loc('wc')]/100))
                    
                    Gas_mass_rate=Gas_mass_rate+pipe_df.iloc[j,pipe_df.columns.get_loc('SGgas')]*pipe_df.iloc[j,pipe_df.columns.get_loc('qliq')]*(1-pipe_df.iloc[j,pipe_df.columns.get_loc('wc')]/100)*pipe_df.iloc[j,pipe_df.columns.get_loc('gazf')]*1.2/1000
                    SGgas_sum=SGgas_sum+pipe_df.iloc[j,pipe_df.columns.get_loc('SGgas')]*(pipe_df.iloc[j,pipe_df.columns.get_loc('SGgas')]*pipe_df.iloc[j,pipe_df.columns.get_loc('qliq')]*(1-pipe_df.iloc[j,pipe_df.columns.get_loc('wc')]/100)*pipe_df.iloc[j,pipe_df.columns.get_loc('gazf')]*1.2/1000)
                    
                    Water_mass_rate=Water_mass_rate+pipe_df.iloc[j,pipe_df.columns.get_loc('SGwater')]*pipe_df.iloc[j,pipe_df.columns.get_loc('qliq')]*(pipe_df.iloc[j,pipe_df.columns.get_loc('wc')]/100)
                    SGwater_sum=SGwater_sum+pipe_df.iloc[j,pipe_df.columns.get_loc('SGwater')]*(pipe_df.iloc[j,pipe_df.columns.get_loc('SGwater')]*pipe_df.iloc[j,pipe_df.columns.get_loc('qliq')]*(pipe_df.iloc[j,pipe_df.columns.get_loc('wc')]/100))
                    
                    if Oil_mass_rate!=0:
                        pipe_df.iloc[i,pipe_df.columns.get_loc('SGoil')]=SGoil_sum/Oil_mass_rate
                    else:
                        pipe_df.iloc[i,pipe_df.columns.get_loc('SGoil')]=SGoil_default
                    if Gas_mass_rate !=0:
                        pipe_df.iloc[i,pipe_df.columns.get_loc('SGgas')]=SGgas_sum/Gas_mass_rate
                    else:
                        pipe_df.iloc[i,pipe_df.columns.get_loc('SGgas')]=SGgas_default
                    if Water_mass_rate!=0:
                        pipe_df.iloc[i,pipe_df.columns.get_loc('SGwater')]=SGwater_sum/Water_mass_rate
                    else:
                        pipe_df.iloc[i,pipe_df.columns.get_loc('SGwater')]=SGwater_default
                # print(col_start)
                j+=1       
        i+=1
# exl.sheets["pipe_df"].range((2,2)).value=pipe_df

    

#pipe_df.columns = pipe_df.iloc[0]
#pipe_df=pipe_df.drop([0])
#print (pipe_df)
pipe_df = pipe_df.reset_index(drop=True)
# pipe_df.set_index('№', inplace=True)
#wb.sheets["pipelines_sort"].range((2,2)).value=pipe_df
#pipe_list=wb.sheets["pipelines_sort"].range((2,2)).expand().value
# pipe_df=pd.DataFrame(pipe_list)
initial_pipe_df=pipe_df
#print (pipe_df)
# считывание исходных данных из DataFrame
i=0
k=0
j=0
# print (pipe_df)

while i<=pipes_count:

    d_out_m=pipe_df.iloc[i,pipe_df.columns.get_loc('out_dia')]
    WT_m=pipe_df.iloc[i,pipe_df.columns.get_loc('wall_thick')]
    L_m=pipe_df.iloc[i,pipe_df.columns.get_loc('length')]
    Qlsc_m=pipe_df.iloc[i,pipe_df.columns.get_loc('qliq')]
    WC=pipe_df.iloc[i,pipe_df.columns.get_loc('wc')]
    Rp_m=pipe_df.iloc[i,pipe_df.columns.get_loc('gazf')]
    P_m=pipe_df.iloc[i,pipe_df.columns.get_loc('press_end')]
    T_m=pipe_df.iloc[i,pipe_df.columns.get_loc('temp_start')]
    if math.isnan(pipe_df.iloc[i,pipe_df.columns.get_loc('temp_start')])==True and math.isnan(pipe_df.iloc[i,pipe_df.columns.get_loc('temp_end')])!=True:
        T_m=pipe_df.iloc[i,pipe_df.columns.get_loc('temp_end')]
    if math.isnan(T_m)==True:
        T_m=t_initial_guess
    start=pipe_df.iloc[i,pipe_df.columns.get_loc('start_point')]
    end=pipe_df.iloc[i,pipe_df.columns.get_loc('end_point')]
    ########## для подсчета плотностей начало
    # SGoil=pipe_df.iloc[i,23]
    # SGgas=pipe_df.iloc[i,24]
    # SGwater=pipe_df.iloc[i,25]
    SGoil=pipe_df.iloc[i,pipe_df.columns.get_loc('SGoil')]
    SGgas=pipe_df.iloc[i,pipe_df.columns.get_loc('SGgas')]
    SGwater=pipe_df.iloc[i,pipe_df.columns.get_loc('SGwater')]
    ########## для подсчета плотностей конец
    pipe_name=start+" - "+end
    if math.isnan(pipe_df.iloc[i,pipe_df.columns.get_loc('height_drop')])==True:
        incl_angle_deg=0
    else:
        incl_angle_deg=math.atan(pipe_df.iloc[i,pipe_df.columns.get_loc('height_drop')]/pipe_df.iloc[i,pipe_df.columns.get_loc('length')])*180/math.pi
    

    a=0
    j=0
    if math.isnan(Qlsc_m)==True and math.isnan(WC)==True:
        # материальный баланс для коллектора, расчет температур, определение давления, газового фактора и обводненности для смеси
        Qlsc_m=0
        Qosc_m=0
        Qgsc_m=0
        P_m=0
        total_mass=0
        mass=array('f')
        t=array('f')
        if math.isnan(pipe_df.iloc[i,pipe_df.columns.get_loc('temp_start')])==True:
            T_m=0

        for row in pipe_df.iterrows():
            
            if pipe_df.iloc[j,pipe_df.columns.get_loc('end_point')]==start:
                
                # объемные расходы по фазам
                Qlsc_m=Qlsc_m+pipe_df.iloc[j,pipe_df.columns.get_loc('qliq')]
                Qosc_m=Qosc_m+pipe_df.iloc[j,pipe_df.columns.get_loc('qliq')]*(1-pipe_df.iloc[j,pipe_df.columns.get_loc('wc')]/100)
                Qgsc_m=Qgsc_m+pipe_df.iloc[j,pipe_df.columns.get_loc('qliq')]*(1-pipe_df.iloc[j,pipe_df.columns.get_loc('wc')]/100)*pipe_df.iloc[j,pipe_df.columns.get_loc('gazf')]
                
                # определение массы потоков
                mass.insert(a,(pipe_df.iloc[j,pipe_df.columns.get_loc('qliq')]*(1-pipe_df.iloc[j,pipe_df.columns.get_loc('wc')]/100)*(SGoil+pipe_df.iloc[j,pipe_df.columns.get_loc('gazf')]*SGgas*1.2/1000) + pipe_df.iloc[j,pipe_df.columns.get_loc('qliq')]*pipe_df.iloc[j,pipe_df.columns.get_loc('wc')]/100*SGwater))
                total_mass=total_mass+mass[a]

                if math.isnan(pipe_df.iloc[j,pipe_df.columns.get_loc('temp_end')])==True:
                    t.insert(a,t_initial_guess)
                else:
                    t.insert(a,pipe_df.iloc[j,pipe_df.columns.get_loc('temp_end')])
                
                if math.isnan(pipe_df.iloc[i,pipe_df.columns.get_loc('press_end')])==True or pipe_df.iloc[i,pipe_df.columns.get_loc('press_end')]==0:
                    if P_m!=0:
                        P_m=min(P_m,pipe_df.iloc[j,pipe_df.columns.get_loc('press_end')])
                    else:
                        
                        P_m=0
                
                a+=1
            j+=1
        a=0
        if math.isnan(pipe_df.iloc[i,pipe_df.columns.get_loc('temp_start')])==True:
            for element in mass:
                T_m = T_m + mass[a] / total_mass * t[a]
                a+=1
        if math.isnan(Rp_m)==True:
            Rp_m=Qgsc_m/Qosc_m
        WC=(Qlsc_m-Qosc_m)/Qlsc_m*100
        if math.isnan(pipe_df.iloc[i,pipe_df.columns.get_loc('press_end')])==True:
            pipe_df.iloc[i,pipe_df.columns.get_loc('press_end')]=P_m
        else:
            P_m=pipe_df.iloc[i,pipe_df.columns.get_loc('press_end')]
        

    j=0
    if math.isnan(pipe_df.iloc[i,pipe_df.columns.get_loc('press_end')])==True:
        
        for row in pipe_df.iterrows():
            if pipe_df.iloc[j,pipe_df.columns.get_loc('start_point')]==end:
                P_m=pipe_df.iloc[j,pipe_df.columns.get_loc('press_start')]
            j+=1



    
    
    
    # перевод единиц измерения исходных данных в расчетную систему (field)
    units_conversion(Qlsc_m, Rp_m, Tsep_m, Psep_m, P_m, T_m, d_out_m, WT_m, roughness_m, incl_angle_deg, L_m, Te_C, BD_m)


    ####################################################
    # расчет давлений и температур по длине трубопровода

    # опрелеление шага по длине
    incr_L=L/increment

    distance=incr_L
    incr=1

    direction=-1
    while incr<=increment:


        # T=T_m*9/5+32
        oil_Rs (T, P, SGoil, SGgas, Rp, Tsep, Psep)
        
        oil_FVF (T, P, y_api, y_g100, Rs, C1_VB2, C2_VB2, C3_VB2, C1_VB3, C2_VB3, C3_VB3)
        
        gas_free_SG (SGgas, Rs, y_api)
        
        gas_Z_crit_P_T (T, P, y_gf)

        gas_FVF_density (T, P, Z, y_gf)

        oil_density(P,SGoil,Rp,Rs,y_gt,Bob,y_gd, Pb)
        
        oil_viscosity (T, y_api, Rs)
        
        gas_viscosity (T, y_gf)
        
        oil_water_surface_tension (T, P, y_api)
        
        water_props (T, P, SGwater)

        oil_liq_gas_balance (Qlsc, WC, Rp, Rs, Rsw, Bo, Bg, Bw, ro_o, ro_w, ro_g)
        
        thermal_props (T, SGoil, y_api, y_gf)

        mixture_props (WC, ro_o, ro_w, ro_g, f_o, f_w, surf_tens_o, surf_tens_w, Cp_o, Cp_w, Cp_g, Xm_o, Xm_g, lambda_l)

        non_slip_velosities (d, q_l, q_g, surf_tens_liq, ro_l, g)
        #print(pipe_df)
        #print(start,'-',end)
        flow_pattern_prediction (lambda_l, Nfr)
        
        holdup_prediction (fl_ptrn_num, Nfr, Nlv, lambda_l, incl_angle_deg, incl_angle_rad, ro_l, v_sl, surf_tens_liq,L1, L2, L3, L4)
        
        slip_velocities (ro_l, ro_g, mu_l, mu_g, HL_incl)

        friction_factor (d, roughness, Nre, lambda_l, HL_incl)

        pressure_gradient (P, d, frict_factor, incl_angle_rad, v_m, v_sg, ro_n, ro_s, g)
        
        pressure_loss (P, incr_L, pressure_grad_psift, direction)
        
        JT_coefficient (T, q_l, Cp_g, ro_g)

        if math.isnan(pipe_df.iloc[i,pipe_df.columns.get_loc('temp_start')])==True and math.isnan(pipe_df.iloc[i,pipe_df.columns.get_loc('temp_end')])!=True:
            fl_temp_est (T,Te,BD,d,d_out,distance,incl_angle_rad,mass_rate_m,Cp_m,mu_n,ro_n, k_f, k_e, v_m, Nre, frict_factor, incr_L, g, JT_mix,pressure_grad_psift,gG)
            start_temp_C=T_f_C
        else:
            temp_estimation (T,Te,BD,d,d_out,distance,incl_angle_rad,mass_rate_m,Cp_m,mu_n,ro_n, k_f, k_e, v_m, Nre, frict_factor,incr_L,g,JT_mix,pressure_grad_psift,gG)
            start_temp_C=0

        bottlenecks (ro_n, v_m_ms,Qlsc_m,mu_o,ro_o,WC,surf_tens_w,surf_tens_o,ro_w,mu_w,iversion_point,q_l,q_g)





        if Pin<0:
            print (pipe_name)
            print ("Отрицательное давление на расстоянии (м) =", distance_m)
            print ("Расчет не выполнен")
            break
            
        
        Pout=P
        Pout_atm=Pout/14.7
        

        
        # print(round(distance_m,1), round(Pin_atm,2), round(Pout_atm,2), round(pressure_grad_atmm,3), round(v_l_ms,1), round(v_g_ms,1),sep="  ", file=f)

        # сбор результатов в массив
        if math.isnan(Qlsc_m)==True:
            Qlsc_m=0
        # distance_m=distance*3.28084
        if math.isnan(pipe_df.iloc[i,pipe_df.columns.get_loc('temp_start')])==True and math.isnan(pipe_df.iloc[i,pipe_df.columns.get_loc('temp_end')])!=True:
            res_arr = (np.append(res_arr,[round(distance_m,1), round(d_out*12*25.4,1), round(WT*12*25.4,1), round(Qlsc*0.159,1),round(Rp/5.6,1),round(WC,1), round(Pin_atm,2), round(Pout_atm,2),

                                        round(T_f_C,1),round(T_C,1), flow_pattern, round(v_l_ms,2), round(v_g_ms,2),round(v_m_ms,2), round(pressure_grad_atmm,6),

                                        round(OHTC_w_m_K,1),round(Nre,1), round(HL_incl,2),round(lambda_l,2),round(h_soil_WmC,1),round(h_f_WmC,1), round(Npr,2), 

                                        round(Nnu,2),round(f_f_ratio,2),round(Rs/5.6,1), round(Rsw/5.6,1), round (EV,1), comment]))
        else:
            res_arr = (np.append(res_arr,[round(distance_m,1), round(d_out*12*25.4,1), round(WT*12*25.4,1), round(Qlsc*0.159,1),round(Rp/5.6,1),round(WC,1), round(Pin_atm,2), round(Pout_atm,2),

                                        round(T_C,1),round(T_f_C,1), flow_pattern, round(v_l_ms,2), round(v_g_ms,2),round(v_m_ms,2), round(pressure_grad_atmm,6),

                                        round(OHTC_w_m_K,1),round(Nre,1), round(HL_incl,2),round(lambda_l,2),round(h_soil_WmC,1),round(h_f_WmC,1), round(Npr,2), 

                                        round(Nnu,2),round(f_f_ratio,2),round(Rs/5.6,1), round(Rsw/5.6,1), round (EV,1), comment]))

        prop_arr=np.append(prop_arr,[round(distance_m,1), round(Pin_atm,2), round(Pout_atm,2), round(T_C,1),round(T_f_C,1),Nre,mu_n,mu_s,ro_n,ro_s,v_m,mu_l,mu_g,h_soil,h_f,OHTC, E_k])
        
        P=Pin
        T=T_f
        
        distance=distance+incr_L
        incr=incr+1
        
        
    else:
        
        pipe_df.iloc[i,pipe_df.columns.get_loc('press_start')]=Pin_atm
        if math.isnan(pipe_df.iloc[i,pipe_df.columns.get_loc('temp_start')])==True and math.isnan(pipe_df.iloc[i,pipe_df.columns.get_loc('temp_end')])!=True:
            
            pipe_df.iloc[i,pipe_df.columns.get_loc('temp_start')]=T_f_C
            pipe_df.iloc[i,pipe_df.columns.get_loc('temp_end')]=T_m
        else:
            pipe_df.iloc[i,pipe_df.columns.get_loc('temp_start')]=T_m
            pipe_df.iloc[i,pipe_df.columns.get_loc('temp_end')]=T_f_C
        pipe_df.iloc[i,pipe_df.columns.get_loc('mix_speed_avg')]=v_m_ms
        pipe_df.iloc[i,pipe_df.columns.get_loc('fluid_speed')]=v_l_ms
        pipe_df.iloc[i,pipe_df.columns.get_loc('gaz_speed')]=v_g_ms
        pipe_df.iloc[i,pipe_df.columns.get_loc('flow_type')]=flow_pattern
        pipe_df.iloc[i,pipe_df.columns.get_loc('press_change')]=(pd.to_numeric(pipe_df.iloc[i,pipe_df.columns.get_loc('press_start')])-pd.to_numeric(pipe_df.iloc[i,pipe_df.columns.get_loc('press_end')]))/(pipe_df.iloc[i,pipe_df.columns.get_loc('length')]/1000)

        # print (pipe_df)

        if math.isnan(pipe_df.iloc[i,pipe_df.columns.get_loc('qliq')])==True:
            pipe_df.iloc[i,pipe_df.columns.get_loc('qliq')]=round(Qlsc_m,2)
            pipe_df.iloc[i,pipe_df.columns.get_loc('wc')]=round(WC,2)
            pipe_df.iloc[i,pipe_df.columns.get_loc('gazf')]=round(Rp_m,2)

        
        res_reshape=np.reshape(res_arr,(-1,28))
        prop_reshape=np.reshape(prop_arr,(-1,17))
        Result_df = (pd.DataFrame(data=res_reshape,columns=['Distance (m)','Dout (mm)','WT','Liq.rate(m3/d)','GOR (m3/m3)','WC (%)', 'Pin (atm)', 'Pout (atm)',

                                                            'Tin (C)','Tout (C)','Flow pattern','Vliq (m/s)','Vgas (m/s)','Vm (m/s)','Pressure gradient (atm/m)',

                                                            'OHTC (W/m2*C)','Nre','Holdup','Lambda', 'h_soil (W/m2*C)','h_f (W/m2*C)','Npr','Nnu','f_f_ratio',
                                                            
                                                            'Rs (m3/m3)','Rsw (m3/m3)','EV (m/s)','Comment']))

        if math.isnan(pipe_df.iloc[i,pipe_df.columns.get_loc('temp_start')])==True and math.isnan(pipe_df.iloc[i,pipe_df.columns.get_loc('temp_end')])!=True:
            Result_df_reverse=Result_df
            Result_df_pressure=Result_df_reverse.iloc[:, 6:7]
            Result_df_rest=Result_df_reverse.iloc[:, 8:27]
            Result_df_rest_reverse=Result_df_rest.sort_index(ascending=False, ignore_index=True)

            Result_df_pressure_reverse=Result_df_pressure.sort_index(ascending=False, ignore_index=True)      
            
            Result_df_reverse.iloc[:, 6:7] = Result_df_pressure_reverse
            Result_df_reverse.iloc[:, 8:27] = Result_df_rest_reverse
        else:
            Result_df_reverse=Result_df
        Prop_df = (pd.DataFrame(data=prop_reshape,columns=['Distance (m)', 'Pin (atm)', 'Pout (atm)','Tin (C)', 'Tout (C)','Nre','mu_n','mu_s','ro_n','ro_s','v_m','mu_l','mu_g',
                                                            'h_soil','h_f','OHTC','E_k']))
        
        # print (Result_df_reverse)
        res_arr = np.array([])

        
        # print (Result_df)
        #print (Result_df_reverse)
        #Result_df_reverse['pipe_name']=pipe_name
        
        #Result_df_reverse.to_csv('reverse.csv', mode='a', index=False)

    # print(pipe_name)
    # print(pipe_df[:20])
    # print (Result_df_reverse)

    #print (pipe_df)
    # сразу же пересчитываем трубопровод на второй раз c новыми температурами
    res_arr2 = np.array([])
    prop_arr2 = np.array([])
    incr_L=L/increment
    incr=1
    distance=incr_L
    P=P_m*14.7
    while incr<=increment:
        # print(increment)
        # print(incr)
        
        # T_m=(pd.to_numeric(Result_df_reverse.iloc[increment-incr,8])+pd.to_numeric(Result_df_reverse.iloc[increment-incr,9]))/2

        T_m=pd.to_numeric(Result_df_reverse.iloc[increment-incr,Result_df_reverse.columns.get_loc('Tin (C)')])
        T=T_m*9/5+32

        oil_Rs (T, P, SGoil, SGgas, Rp, Tsep, Psep)
        
        oil_FVF (T, P, y_api, y_g100, Rs, C1_VB2, C2_VB2, C3_VB2, C1_VB3, C2_VB3, C3_VB3)
        
        gas_free_SG (SGgas, Rs, y_api)
        
        gas_Z_crit_P_T (T, P, y_gf)

        gas_FVF_density (T, P, Z, y_gf)

        oil_density(P,SGoil,Rp,Rs,y_gt,Bob,y_gd, Pb)
        
        oil_viscosity (T, y_api, Rs)
        
        gas_viscosity (T, y_gf)
        
        oil_water_surface_tension (T, P, y_api)
        
        water_props (T, P, SGwater)

        oil_liq_gas_balance (Qlsc, WC, Rp, Rs, Rsw, Bo, Bg, Bw, ro_o, ro_w, ro_g)
        
        thermal_props (T, SGoil, y_api, y_gf)

        mixture_props (WC, ro_o, ro_w, ro_g, f_o, f_w, surf_tens_o, surf_tens_w, Cp_o, Cp_w, Cp_g, Xm_o, Xm_g, lambda_l)

        non_slip_velosities (d, q_l, q_g, surf_tens_liq, ro_l, g)

        flow_pattern_prediction (lambda_l, Nfr)
        
        holdup_prediction (fl_ptrn_num, Nfr, Nlv, lambda_l, incl_angle_deg, incl_angle_rad, ro_l, v_sl, surf_tens_liq,L1, L2, L3, L4)
        
        slip_velocities (ro_l, ro_g, mu_l, mu_g, HL_incl)

        friction_factor (d, roughness, Nre, lambda_l, HL_incl)

        pressure_gradient (P, d, frict_factor, incl_angle_rad, v_m, v_sg, ro_n, ro_s, g)
        
        pressure_loss (P, incr_L, pressure_grad_psift, direction)
        
        JT_coefficient (T, q_l, Cp_g, ro_g)
        
        temp_estimation (T,Te,BD,d,d_out,distance,incl_angle_rad,mass_rate_m,Cp_m,mu_n,ro_n, k_f, k_e, v_m, Nre, frict_factor,incr_L,g,JT_mix,pressure_grad_psift,gG)

        bottlenecks (ro_n, v_m_ms,Qlsc_m,mu_o,ro_o,WC,surf_tens_w,surf_tens_o,ro_w,mu_w,iversion_point,q_l,q_g)
        


        if Pin<0:
            print (pipe_name)
            print ("Отрицательное давление на расстоянии (м) =", distance_m)
            print ("Расчет не выполнен")
            break
            
        
        Pout=P
        Pout_atm=Pout/14.7
        

        
        # print(round(distance_m,1), round(Pin_atm,2), round(Pout_atm,2), round(pressure_grad_atmm,3), round(v_l_ms,1), round(v_g_ms,1),sep="  ", file=f)

        # сбор результатов в массив
        if math.isnan(Qlsc_m)==True:
            Qlsc_m=0

        
        res_arr2 = (np.append(res_arr2,[round(distance_m,1), round(d_out*12*25.4,1), round(WT*12*25.4,1), round(Qlsc*0.159,1),round(Rp/5.6,1),round(WC,1), round(Pin_atm,2), round(Pout_atm,2),

                                    round(T_C,1),round(T_f_C,1), flow_pattern, round(v_l_ms,2), round(v_g_ms,2),round(v_m_ms,2), round(pressure_grad_atmm,6),

                                    round(OHTC_w_m_K,1),round(Nre,1), round(HL_incl,2),round(lambda_l,2),round(h_soil_WmC,1),round(h_f_WmC,1), round(Npr,2), 

                                    round(Nnu,2),round(f_f_ratio,2),round(Rs/5.6,1), round(Rsw/5.6,1), round (EV,1), comment]))

        prop_arr2=np.append(prop_arr2,[round(distance_m,1), round(Pin_atm,2), round(Pout_atm,2), round(T_C,1),round(T_f_C,1),Nre,mu_n,mu_s,ro_n,ro_s,v_m,mu_l,mu_g,h_soil,h_f,OHTC, E_k])
        
        P=Pin
        T=T_f
        
        
        distance=distance+incr_L
        incr=incr+1
        
        


    else:
        
        pipe_df.iloc[i,pipe_df.columns.get_loc('press_start')]=Pin_atm
        pipe_df.iloc[i,pipe_df.columns.get_loc('press_end')]=P_m
        # pipe_df.iloc[i,pipe_df.columns.get_loc('temp_start')]=T_m
        if start_temp_C!=0:
            pipe_df.iloc[i,pipe_df.columns.get_loc('temp_start')]=start_temp_C
        else:
           pipe_df.iloc[i,pipe_df.columns.get_loc('temp_start')]=T_m
        
        pipe_df.iloc[i,pipe_df.columns.get_loc('mix_speed_avg')]=v_m_ms
        pipe_df.iloc[i,pipe_df.columns.get_loc('fluid_speed')]=v_l_ms
        pipe_df.iloc[i,pipe_df.columns.get_loc('gaz_speed')]=v_g_ms
        pipe_df.iloc[i,pipe_df.columns.get_loc('flow_type')]=flow_pattern
        pipe_df.iloc[i,pipe_df.columns.get_loc('press_change')]=(pd.to_numeric(pipe_df.iloc[i,pipe_df.columns.get_loc('press_start')])-pd.to_numeric(pipe_df.iloc[i,pipe_df.columns.get_loc('press_end')]))/(pipe_df.iloc[i,pipe_df.columns.get_loc('length')]/1000)


        if math.isnan(pipe_df.iloc[i,pipe_df.columns.get_loc('qliq')])==True:
            pipe_df.iloc[i,pipe_df.columns.get_loc('qliq')]=round(Qlsc_m,2)
            pipe_df.iloc[i,pipe_df.columns.get_loc('wc')]=round(WC,2)
            pipe_df.iloc[i,pipe_df.columns.get_loc('gazf')]=round(Rp_m,2)

        
        res_reshape2=np.reshape(res_arr2,(-1,28))
        prop_reshape2=np.reshape(prop_arr2,(-1,17))
        Result_df2 = (pd.DataFrame(data=res_reshape2,columns=['Distance (m)','Dout (mm)','WT','Liq.rate(m3/d)','GOR (m3/m3)','WC (%)', 'Pin (atm)', 'Pout (atm)',

                                                            'Tin (C)','Tout (C)','Flow pattern','Vliq (m/s)','Vgas (m/s)','Vm (m/s)','Pressure gradient (atm/m)',

                                                            'OHTC (W/m2*C)','Nre','Holdup','Lambda', 'h_soil (W/m2*C)','h_f (W/m2*C)','Npr','Nnu','f_f_ratio',
                                                            
                                                            'Rs (m3/m3)','Rsw (m3/m3)','EV (m/s)','Comment']))

        if math.isnan(pipe_df.iloc[i,pipe_df.columns.get_loc('temp_start')])==True and math.isnan(pipe_df.iloc[i,pipe_df.columns.get_loc('temp_end')])!=True:
            Result_df_reverse2=Result_df2
            Result_df_pressure2=Result_df_reverse2.iloc[:, 6:8]
            Result_df_rest2=Result_df_reverse2.iloc[:, 8:27]
            #print (Result_df2)
            Result_df_pressure_reverse2=Result_df_pressure2.sort_index(ascending=False, ignore_index=True)
            Result_df_rest_reverse2=Result_df_rest2.sort_index(ascending=False, ignore_index=True)
            
            Result_df_reverse2.iloc[:, 6:8] = Result_df_pressure_reverse2
            Result_df_reverse2.iloc[:, 8:27] = Result_df_rest_reverse2
            pipe_df.iloc[i,pipe_df.columns.get_loc('temp_end')]=Result_df_reverse2.iloc[9, 9]
            #print (Result_df_reverse2)
            #print (pipe_df)
        else:
            Result_df_reverse2=Result_df2

        Prop_df2 = (pd.DataFrame(data=prop_reshape2,columns=['Distance (m)', 'Pin (atm)', 'Pout (atm)','Tin (C)', 'Tout (C)','Nre','mu_n','mu_s','ro_n','ro_s','v_m','mu_l','mu_g',
                                                            'h_soil','h_f','OHTC','E_k']))
        
        
        res_arr2 = np.array([])

    
    i=i+1

#print (pipe_df)
pipe_df_reverse=pipe_df
pipe_df_reverse=pipe_df_reverse.iloc[1:, :]
     
pipe_df_reverse=pipe_df_reverse.sort_index(ascending=False, ignore_index=True)

pipe_df=pipe_df.iloc[0]
pipe_df=pd.DataFrame(pipe_df)
pipe_df=pipe_df.T
pipe_df=pipe_df.append(pipe_df_reverse, ignore_index=True)
# wb.sheets["pipe_df"].range((2,2)).value=pipe_df
#print (pipe_df)

###################################################################################
# третий и последующие проходы
iteration=0
while iteration<max_iter:
    res_arr2 = np.array([])
    prop_arr2 = np.array([])
    # проверка и пересчет температуры
    # считывание исходных данных из DataFrame
    frames_ds=[]
    df_columns=['ID','Distance (m)','Dout (mm)','WT','Liq.rate(m3/d)','GOR (m3/m3)','WC (%)', 'Pin (atm)', 'Pout (atm)',

                                                                'Tin (C)','Tout (C)','Flow pattern','Vliq (m/s)','Vgas (m/s)','Vm (m/s)','Pressure gradient (atm/m)',

                                                                'OHTC (W/m2*C)','Nre','Holdup','Lambda', 'h_soil (W/m2*C)','h_f (W/m2*C)','Npr','Nnu','f_f_ratio',
                                                                
                                                                'Rs (m3/m3)','Rsw (m3/m3)','EV (m/s)','Comment']
    i=0
    k=0
    j=0
    

    while i<=pipes_count:
        pipe_ID=pipe_df.iloc[i,pipe_df.columns.get_loc('ID')]
        d_out_m=pipe_df.iloc[i,pipe_df.columns.get_loc('out_dia')]
        WT_m=pipe_df.iloc[i,pipe_df.columns.get_loc('wall_thick')]
        L_m=pipe_df.iloc[i,pipe_df.columns.get_loc('length')]
        Qlsc_m=pipe_df.iloc[i,pipe_df.columns.get_loc('qliq')]
        WC=pipe_df.iloc[i,pipe_df.columns.get_loc('wc')]
        Rp_m=pipe_df.iloc[i,pipe_df.columns.get_loc('gazf')]
        P_m=pipe_df.iloc[i,pipe_df.columns.get_loc('press_end')]
        P_m_1=P_m
        T_m=pipe_df.iloc[i,pipe_df.columns.get_loc('temp_start')]
        T_m_1=T_m
        T_m_final=pipe_df.iloc[i,pipe_df.columns.get_loc('temp_end')]
        
        start=pipe_df.iloc[i,pipe_df.columns.get_loc('start_point')]
        end=pipe_df.iloc[i,pipe_df.columns.get_loc('end_point')]
        pipe_name=start+" - "+end

        if math.isnan(pipe_df.iloc[i,pipe_df.columns.get_loc('height_drop')])==True:
            incl_angle_deg=0
        else:
            incl_angle_deg=math.atan(pipe_df.iloc[i,pipe_df.columns.get_loc('height_drop')]/pipe_df.iloc[i,pipe_df.columns.get_loc('length')])*180/math.pi

        a=0
        j=0
        # if Qlsc_m==None and WC==None:
        #     # материальный баланс для коллектора, расчет температур, определение давления, газового фактора и обводненности для смеси
        #     Qlsc_m=0
        #     Qosc_m=0
        #     Qgsc_m=0
        #     P_m=0
        total_mass=0
        mass=array('f')
        t=array('f')
            # if pipe_df.iloc[i,9]==None:
            #     T_m=0

        # перевод единиц измерения исходных данных в расчетную систему (field)
        units_conversion(Qlsc_m, Rp_m, Tsep_m, Psep_m, P_m, T_m, d_out_m, WT_m, roughness_m, incl_angle_deg, L_m, Te_C, BD_m)


        for row in pipe_df.iterrows():
            
            if pipe_df.iloc[j,pipe_df.columns.get_loc('end_point')]==start:

                T_m_1=0     
                    
                    
                # определение массы потоков
                mass.insert(a,(pipe_df.iloc[j,pipe_df.columns.get_loc('qliq')]*(1-pipe_df.iloc[j,pipe_df.columns.get_loc('wc')]/100)*(SGoil+pipe_df.iloc[j,pipe_df.columns.get_loc('gazf')]*SGgas*1.2/1000) + pipe_df.iloc[j,pipe_df.columns.get_loc('qliq')]*pipe_df.iloc[j,pipe_df.columns.get_loc('wc')]/100*SGwater))
                total_mass=total_mass+mass[a]

                if pipe_df.iloc[j,pipe_df.columns.get_loc('temp_end')]==None:
                    t.insert(a,t_initial_guess)
                else:
                    t.insert(a,pd.to_numeric(pipe_df.iloc[j,pipe_df.columns.get_loc('temp_end')]))
    
                a+=1
            
            if pipe_df.iloc[j,pipe_df.columns.get_loc('start_point')]==end:
                P_m_1=pipe_df.iloc[j,pipe_df.columns.get_loc('press_start')]

            

            j+=1

        a=0
        
        for element in mass:
            T_m_1 = T_m_1 + mass[a] / total_mass * t[a]
            a+=1
        
        # if abs(T_m_1-T_m)/abs(T_m)>0.01 or abs(P_m_1-P_m)/abs(P_m)>0.01 :
        T_m=T_m_1
        P_m=P_m_1
        incr=1
        incr_L=L/increment
        distance=incr_L
        P=P_m*14.7
        T=T_m*9/5+32
        
        while incr<=increment:            

            oil_Rs (T, P, SGoil, SGgas, Rp, Tsep, Psep)
            
            oil_FVF (T, P, y_api, y_g100, Rs, C1_VB2, C2_VB2, C3_VB2, C1_VB3, C2_VB3, C3_VB3)
            
            gas_free_SG (SGgas, Rs, y_api)
            
            gas_Z_crit_P_T (T, P, y_gf)

            gas_FVF_density (T, P, Z, y_gf)

            oil_density(P,SGoil,Rp,Rs,y_gt,Bob,y_gd, Pb)
            
            oil_viscosity (T, y_api, Rs)
            
            gas_viscosity (T, y_gf)
            
            oil_water_surface_tension (T, P, y_api)
            
            water_props (T, P, SGwater)

            oil_liq_gas_balance (Qlsc, WC, Rp, Rs, Rsw, Bo, Bg, Bw, ro_o, ro_w, ro_g)
            
            thermal_props (T, SGoil, y_api, y_gf)

            mixture_props (WC, ro_o, ro_w, ro_g, f_o, f_w, surf_tens_o, surf_tens_w, Cp_o, Cp_w, Cp_g, Xm_o, Xm_g, lambda_l)

            non_slip_velosities (d, q_l, q_g, surf_tens_liq, ro_l, g)

            flow_pattern_prediction (lambda_l, Nfr)
            
            holdup_prediction (fl_ptrn_num, Nfr, Nlv, lambda_l, incl_angle_deg, incl_angle_rad, ro_l, v_sl, surf_tens_liq,L1, L2, L3, L4)
            
            slip_velocities (ro_l, ro_g, mu_l, mu_g, HL_incl)

            friction_factor (d, roughness, Nre, lambda_l, HL_incl)

            pressure_gradient (P, d, frict_factor, incl_angle_rad, v_m, v_sg, ro_n, ro_s, g)
            
            pressure_loss (P, incr_L, pressure_grad_psift, direction)
            
            JT_coefficient (T, q_l, Cp_g, ro_g)
            
            temp_estimation (T,Te,BD,d,d_out,distance,incl_angle_rad,mass_rate_m,Cp_m,mu_n,ro_n, k_f, k_e, v_m, Nre, frict_factor,incr_L,g,JT_mix,pressure_grad_psift,gG)

            bottlenecks (ro_n, v_m_ms,Qlsc_m,mu_o,ro_o,WC,surf_tens_w,surf_tens_o,ro_w,mu_w,iversion_point,q_l,q_g)
            


            if Pin<0:
                print (pipe_name)
                print ("Отрицательное давление на расстоянии (м) =", distance_m)
                print ("Расчет не выполнен")
                break
                
            
            Pout=P
            Pout_atm=Pout/14.7
            

            
            # print(round(distance_m,1), round(Pin_atm,2), round(Pout_atm,2), round(pressure_grad_atmm,3), round(v_l_ms,1), round(v_g_ms,1),sep="  ", file=f)

            # сбор результатов в массив
            if math.isnan(Qlsc_m)==True:
                Qlsc_m=0

            
            res_arr2 = (np.append(res_arr2,[pipe_ID, round(distance_m,1), round(d_out*12*25.4,1), round(WT*12*25.4,1), round(Qlsc*0.159,1),round(Rp/5.6,1),round(WC,1), round(Pin_atm,2), round(Pout_atm,2),

                                        round(T_C,1),round(T_f_C,1), flow_pattern, round(v_l_ms,2), round(v_g_ms,2),round(v_m_ms,2), round(pressure_grad_atmm,6),

                                        round(OHTC_w_m_K,1),round(Nre,1), round(HL_incl,2),round(lambda_l,2),round(h_soil_WmC,1),round(h_f_WmC,1), round(Npr,2), 

                                        round(Nnu,2),round(f_f_ratio,2),round(Rs/5.6,1), round(Rsw/5.6,1), round (EV,1), comment]))

            prop_arr2=np.append(prop_arr2,[round(distance_m,1), round(Pin_atm,2), round(Pout_atm,2), round(T_C,1),round(T_f_C,1),Nre,mu_n,mu_s,ro_n,ro_s,v_m,mu_l,mu_g,h_soil,h_f,OHTC, E_k])
            
            P=Pin
            T=T_f
            
            
            distance=distance+incr_L
            incr=incr+1
            


        else:
            
            pipe_df.iloc[i,pipe_df.columns.get_loc('press_start')]=Pin_atm
            pipe_df.iloc[i,pipe_df.columns.get_loc('press_end')]=P_m
            pipe_df.iloc[i,pipe_df.columns.get_loc('temp_start')]=T_m
            # pipe_df.iloc[i,10]=T_f_C
            pipe_df.iloc[i,pipe_df.columns.get_loc('mix_speed_avg')]=v_m_ms
            pipe_df.iloc[i,pipe_df.columns.get_loc('fluid_speed')]=v_l_ms
            pipe_df.iloc[i,pipe_df.columns.get_loc('gaz_speed')]=v_g_ms
            pipe_df.iloc[i,pipe_df.columns.get_loc('flow_type')]=flow_pattern
            pipe_df.iloc[i,pipe_df.columns.get_loc('press_change')]=(pd.to_numeric(pipe_df.iloc[i,pipe_df.columns.get_loc('press_start')])-pd.to_numeric(pipe_df.iloc[i,pipe_df.columns.get_loc('press_end')]))/(pipe_df.iloc[i,pipe_df.columns.get_loc('length')]/1000)


            if math.isnan(pipe_df.iloc[i,pipe_df.columns.get_loc('qliq')])==True:
                pipe_df.iloc[i,pipe_df.columns.get_loc('qliq')]=round(Qlsc_m,2)
                pipe_df.iloc[i,pipe_df.columns.get_loc('wc')]=round(WC,2)
                pipe_df.iloc[i,pipe_df.columns.get_loc('gazf')]=round(Rp_m,2)

            
            res_reshape2=np.reshape(res_arr2,(-1,29))
            prop_reshape2=np.reshape(prop_arr2,(-1,17))
            Result_df2 = (pd.DataFrame(data=res_reshape2,columns=['ID', 'Distance (m)','Dout (mm)','WT','Liq.rate(m3/d)','GOR (m3/m3)','WC (%)', 'Pin (atm)', 'Pout (atm)',

                                                                'Tin (C)','Tout (C)','Flow pattern','Vliq (m/s)','Vgas (m/s)','Vm (m/s)','Pressure gradient (atm/m)',

                                                                'OHTC (W/m2*C)','Nre','Holdup','Lambda', 'h_soil (W/m2*C)','h_f (W/m2*C)','Npr','Nnu','f_f_ratio',
                                                                
                                                                'Rs (m3/m3)','Rsw (m3/m3)','EV (m/s)','Comment']))

            
            Result_df_reverse2=Result_df2
            Result_df_pressure2=Result_df_reverse2.iloc[:, 7:9]
            Result_df_rest2=Result_df_reverse2.iloc[:, 11:30]
            
            Result_df_pressure_reverse2=Result_df_pressure2.sort_index(ascending=False, ignore_index=True)
            Result_df_rest_reverse2=Result_df_rest2.sort_index(ascending=False, ignore_index=True)
            
            Result_df_reverse2.iloc[:, 7:9] = Result_df_pressure_reverse2
            Result_df_reverse2.iloc[:, 11:30] = Result_df_rest_reverse2
            #pipe_df.iloc[i,10]=Result_df_reverse2.iloc[9, 9]
            pipe_df.iloc[i,pipe_df.columns.get_loc('temp_end')]=Result_df_reverse2.iloc[increment-1,Result_df_reverse2.columns.get_loc('Tout (C)')]
            

            Prop_df2 = (pd.DataFrame(data=prop_reshape2,columns=['Distance (m)', 'Pin (atm)', 'Pout (atm)','Tin (C)', 'Tout (C)','Nre','mu_n','mu_s','ro_n','ro_s','v_m','mu_l','mu_g',
                                                                'h_soil','h_f','OHTC','E_k']))
            
            #wb.sheets["results"].range((2+k,2)).value=Result_df_reverse2
            #wb.sheets["results"].range((2+k,2)).value=pipe_name
            Result_df_reverse2['pipe_name']=Result_df_reverse2.index
            frames_ds.append({'df':Result_df_reverse2, 'pipe_name':pipe_name})
        
            #Result_df_reverse2.to_csv('reverse.csv', mode='a', index=False)
            k+=increment+2
            res_arr2 = np.array([])
            
        # else:
            
        #     res_arr2 = np.array([])
        
        i=i+1

        
        
        


        
    pipe_df_reverse=pipe_df
    pipe_df_reverse=pipe_df_reverse.iloc[1:, :]
        
    pipe_df_reverse=pipe_df_reverse.sort_index(ascending=False, ignore_index=True)

    pipe_df=pipe_df.iloc[0]
    pipe_df=pd.DataFrame(pipe_df)
    pipe_df=pipe_df.T
    pipe_df=pipe_df.append(pipe_df_reverse, ignore_index=True)

    iteration+=1



#pipe_df.columns = pipe_df.iloc[0]
#pipe_df=pipe_df.drop([0])
# pipe_df.set_index('№', inplace=True)

#pipe_df.to_csv("pipe_df.csv")
#print ("prop_df")
#print(Prop_df2)
#Prop_df2.to_csv("prop_df2.csv")
#wb.sheets["summary"].range((2,2)).value=pipe_df
#print(Result_df)
#Result_df.to_csv("result_df.csv")

# long_results=sys.path[0]+"/long_results.csv"
if os.path.isfile(long_results):
    os.remove(long_results)

for entry in frames_ds:
    entry['df'].to_csv(long_results, mode='a', index_label=entry['pipe_name'],
                       header=True, index=True, columns=df_columns)

# calc_results=sys.path[0]+"/calc_results.csv"
pipe_df.set_index('ID', inplace=True)
pipe_df.to_csv(calc_results)
# Result_df2.to_csv(long_results)
# f.close()
