import pandas as pd
from uncertainties import ufloat
from decimal import Decimal

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

        tabular = self.params_df.to_latex(
            index=False,
            escape=False,
            column_format="@{}r|c@{}"  # optional: Außenabstände entfernen
        )

        self.params_tex = (
            r"\begin{table}[H]" "\n"
            r"\centering" "\n"
            + tabular + "\n"
            rf"\caption{{{cap}}}" "\n"
            rf"\label{{{lab}}}" "\n"
            r"\end{table}"
        )
        
        self.params_tex = (self.params_tex.replace(r"\toprule", r"")
            .replace(r"\midrule", r"\hline")
            .replace(r"\bottomrule", r"")
            .replace(r"+/-", r" \pm "))
    
    def params_export_tex(self, path: str, filename: str):
        with open(path+filename, "w", encoding="utf-8") as f:
            f.write(self.params_tex)