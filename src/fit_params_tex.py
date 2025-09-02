import pandas as pd
from uncertainties import ufloat
from decimal import Decimal


def df_to_tex_df(params_df, cap: str="", lab: str=""):
    if params_df.shape[1] == 2:
        column_format = "r|c"
    elif params_df.shape[1] == 3:
        column_format = "r|c|c"
    else:
        column_format = "r" + "|c" * (params_df.shape[1] - 1)

    tabular = params_df.to_latex(
        index=False,
        escape=False,
        column_format=column_format
    )
    params_tex = (
        r"\begin{table}[H]" "\n"
        r"\centering" "\n"
        + tabular + "\n"
        rf"\caption{{{cap}}}" "\n"
        rf"\label{{{lab}}}" "\n"
        r"\end{table}"
    )
    
    params_tex = (params_tex.replace(r"\toprule", r"")
        .replace(r"\midrule", r"\hline")
        .replace(r"\bottomrule", r"")
        .replace(r"+/-", r" \pm "))
    return params_tex



class FitParametersTex:
    def __init__(self, popt, perr, names, lab: str, cap: str, ndig=3):
        uvals = [ufloat(popt[i], perr[i]) for i in range(len(popt))]

        def sci_notation(val, ndig):
            d = Decimal(str(val))
            exp = d.adjusted()
            coeff = d.scaleb(-exp).normalize()
            coeff_str = f"{float(coeff):.{ndig}g}"
            return coeff_str, exp

        wert_str = []
        for u in uvals:
            coeff_n, exp_n = sci_notation(u.n, ndig)
            coeff_s, exp_s = sci_notation(u.s, ndig)
            wert_str.append(
            rf"$ ({coeff_n} \pm {coeff_s}) \cdot 10^{{{exp_n}}} $"
            )

        self.params_df = pd.DataFrame({
            "Parameter": names,
            "Wert": wert_str
        })
        self.params_tex = df_to_tex_df(self.params_df, cap, lab)
    def params_export_tex(self, path: str, filename: str):
        params_export_tex(self.params_tex, path, filename)
    
    

def params_export_tex(params_tex, path: str, filename: str):
    with open(path+filename, "w", encoding="utf-8") as f:
        f.write(params_tex)
            
    

def combine(param_a, param_b, cap, lab, title_a="Vorl. Strang", title_b="Rückl. Strang"):
    dfs = [param_a.params_df, param_b.params_df]
    combined_df = pd.concat(
        [
            dfs[0].set_index("Parameter")["Wert"].rename(title_a),
            dfs[1].set_index("Parameter")["Wert"].rename(title_b),
        ],
        axis=1
    ).reset_index()

    params_tex = df_to_tex_df(combined_df, cap, lab)
    return params_tex