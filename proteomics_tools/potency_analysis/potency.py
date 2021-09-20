"""
This is a script for computing potency values like the ic50
"""
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.optimize import curve_fit


# I Helper functions
def _exp_reg(x, a=100, b=0.1, c=0):
    """Exponential negative regression function

    Parameters
    ----------
    x: array, shape (n_values, 1)
        x values
    a: float
        Maximum (a+c) for x = 0
    b: float, <1
        Rate constant (high -> faster to minimum)
    c: float
        Minimum

    Returns
    -------
    y: array, shape (n_values, 1)
        y values for f(x)

    See Also
    --------
    https://docs.scipy.org/doc/scipy/reference/generated/scipy.optimize.curve_fit.html
    """
    return a * np.exp(-b * x) + c


def _get_fitted_curve(x=None, y=None, f=None, n=1000, upper_bound_ab=1000):
    """Get fitted curve for x and y"""
    xdata = np.linspace(min(x), max(x), n)   # Created x data for curve
    # Lower bound for all parameters=0
    # Upper bound for a=1000, b=1, b=1000
    f_args, pcov = curve_fit(f, x, y, bounds=(0, [upper_bound_ab, 1, upper_bound_ab]))    # Parameter fitting
    ydata = f(xdata, *f_args)   # Predicted y data for curve
    a, b, c = f_args
    str_formula = '$y=%3.7s*e^{-%3.7sx}$'%(round(a, 1), round(b, 3))
    if c < 0:
        str_formula += '$%3.7s$' % round(c, 1)
    else:
        str_formula += '$+%3.7s$' % round(c, 1)
    return xdata, ydata, str_formula


def _get_x50_y50(xdata=None, ydata=None):
    """Get """
    y50 = min(ydata) + (max(ydata) - min(ydata)) / 2
    y_dif = [np.absolute(y_val - y50) for y_val in ydata]
    x50 = xdata[y_dif.index(min(y_dif))]
    return x50, y50


# II Main Functions
class PotencyAnalysis:
    """Class for potency analysis"""

    @staticmethod
    def get_ic50(df=None, x_label=None, y_label=None):
        """Get ic50 for x depending on y"""
        df = df[[x_label, y_label]].dropna()
        x = np.array(df[x_label])
        y = np.array(df[y_label])
        xdata, ydata, str_formula = _get_fitted_curve(f=_exp_reg, x=x, y=y)
        x50, y50 = _get_x50_y50(xdata=xdata, ydata=ydata)
        return x50, y50

    @staticmethod
    def plot_ic50(ax=None, df=None, x_label=None, y_label=None, legend_label=None,
                  fitted=True, xypos_ic50_f=None,
                  log_scale=True, color=None, ic50_unit="µM",
                  size_label=20, size_legend=16, size_info=18):
        """Plot lineplot for concentration with fitted curve used to determine IC50"""
        if xypos_ic50_f is None:
            xypos_ic50_f = [1, 1.4, 0.5, 0.5]
        xpos_ic50, ypos_ic50, xpos_f, ypos_f = xypos_ic50_f
        # Fitted curve
        df = df[[x_label, y_label]].dropna()
        x_array = np.array(df[x_label])
        y_array = np.array(df[y_label])
        xdata, ydata, str_formula = _get_fitted_curve(f=_exp_reg, x=x_array, y=y_array)
        x50, y50 = _get_x50_y50(xdata=xdata, ydata=ydata)
        # Original
        if legend_label is None:
            legend_label = "data"
        ax = sns.lineplot(ax=ax, data=df, x=x_label, y=y_array, color=color, legend=True,
                          linewidth=3, linestyle="-", label=legend_label, err_style="bars")
        if fitted:
            # Fitted curve
            ax = sns.lineplot(ax=ax, x=xdata, y=ydata, linestyle="--", linewidth=3, color="black", label="fitted")
            if log_scale:
                ax.set(xscale="log")
            # Modify plot
            arrow_properties = dict(facecolor="black", width=0.5, headwidth=4, shrink=0.1)
            plt.annotate(f"IC50 = {round(x50, 1)} {ic50_unit}",
                         (x50, y50),
                         (x50*xpos_ic50, y50*ypos_ic50),
                         horizontalalignment="center",
                         arrowprops=arrow_properties, fontsize=size_info)
            #plt.text(x50, y50, "IC50 = {} µM".format(round(x50, 1)), size=TICK_SIZE)
            plt.legend(fontsize=size_legend)
            plt.text(x50 * xpos_f, y50 * ypos_f, str_formula, fontsize=size_info)
        else:
            ax = plt.gca()
            ax.get_legend().remove()
        plt.ylabel(y_label, size=size_label, weight="bold")
        plt.xlabel(x_label, size=size_label, weight="bold")
        sns.despine()
        return ax
