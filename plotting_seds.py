import numpy as np
import matplotlib.pyplot as plt
from read_prep_data import curved_power_law, power_law
from scipy.stats import chi2
plt.rcParams["font.family"] = "serif"
plt.rcParams["axes.grid"] = False

majorticklength=2.5
minorticklength=2
tickwidth=0.5

def plot_sed(freq, flux, fluxerr, fit_params, coord_src, outpath):
    fig, ax = plt.subplots(1,1)

    nu = np.geomspace(
        np.min(freq),
        np.max(freq),
        100
    )

    ax.errorbar(
        freq,
        flux,
        yerr=fluxerr,
        ls='None',
        marker='.'
    )

    legend = False

    if fit_params is not None:
        legend = True
        
        if len(fit_params) == 2:

            ax.plot(
                nu,
                power_law(nu, *fit_params),
                ls='--',
                color='red',
                label=f"Power-law, r$\alpha$  {fit_params[1]:.3f}"
            )
        elif len(fit_params) == 3: 
            ax.plot(
                nu,
                curved_power_law(nu, *fit_params),
                ls=':',
                color='green',
                label=f'Curved power-law, r$q$ {fit_params[2]:.3f}'
            )


    if legend is True:
        ax.legend()

    ax.loglog()
    ax.set(
        xlabel='Frequency (MHz)',
        ylabel='Integrated Flux (Jy)',
        title=coord_src.tex
    )
    ax.tick_params(
        axis="both", which="major", direction="in", length=majorticklength, width=tickwidth, pad=5
    )
    ax.tick_params(
        axis="both", which="minor", direction="in", length=minorticklength, width=tickwidth
    )

    fig.tight_layout()
    fig.savefig(f"{outpath}/{coord_src}.png")
    plt.close(fig)
