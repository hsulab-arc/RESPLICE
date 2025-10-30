import pandas as pd
import numpy as np
import itertools

df_dict = pd.read_excel("flow/PromoterScreenFinal.xlsx",sheet_name=None,index_col=0)

sheets = list(df_dict.keys())

def ratio_sd(m_g,m_b,var_g,var_b):
    return ((1/3)*((var_g/(m_b**2))+((m_g**2)*var_b**2)/(m_b**4)))**.5

def analyze_df(df):
    green = df.iloc[:,[0,1,2]].values
    blue = df.iloc[:,[3,4,5]].values
    g_mean = np.mean(green,axis=1)
    b_mean = np.mean(blue,axis=1)
    g_var = np.var(green,axis=1,ddof=1)
    b_var = np.var(blue,axis=1,ddof=1)
    ratios = g_mean/b_mean
    ratio_sds = ratio_sd(g_mean,b_mean,g_var,b_var)
    control = ratios[0]
    relative_ratios = ratios/control
    return [ratios,ratio_sds,relative_ratios]
with pd.ExcelWriter("promoterscreen_ratios.xlsx",engine="xlsxwriter") as writer:
    for s in sheets:
        df = df_dict[s]
        index = df.index.tolist()
        analysis = analyze_df(df)
        analyzed_df = pd.DataFrame(np.transpose(analysis),index=index,columns=["Ratio","Std. Dev","Relative Ratio"])
        analyzed_df.to_excel(writer,s,index=True)



