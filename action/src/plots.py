from common import *
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

plt.rcParams.update({'font.size':10,'figure.facecolor':'white'})
early=np.load(RES/'benchmark.npz');late=np.load(RES/'late_slow_envelope.npz')
fig,axes=plt.subplots(2,2,figsize=(11,7),layout='constrained')
axes[0,0].plot(early['states'][0],early['order']);axes[0,0].axhline(1,color='gray',ls=':')
axes[0,0].set(xlabel='N desde el final de inflación',ylabel='Tr(D²) / 0.75',title='Orden Q homogéneo; no mide una red')
axes[0,1].semilogy(early['states'][0],np.maximum(abs(early['constraint'])/early['rho'],1e-16))
axes[0,1].set(xlabel='N desde el final de inflación',ylabel='|3h − ρ| / ρ',title='Restricción, sin relajación')
n=late['N_since_matching'];rho=late['densities'];total=rho.sum(axis=0)
for i,label in enumerate(['χ + Q promediados','radiación visible','Goldstones homogéneos']):
    axes[1,0].semilogy(n,rho[i]/total,label=label)
axes[1,0].legend(fontsize=8);axes[1,0].set(xlabel='N desde el empalme',ylabel='fracción de energía',ylim=(1e-13,2),title='Recalentamiento efectivo; sin cuerdas')
dec=json.loads((RES/'decays.json').read_text())
axes[1,1].bar(np.arange(7),dec['widths_GeV']);axes[1,1].axhline(dec['Gamma_chi_GeV'],color='red',ls='--',label='Γχ')
axes[1,1].set(yscale='log',xlabel='Modo Q (resolviendo degeneración)',ylabel='Γ [GeV]',title='Canales de vacío, a nivel árbol');axes[1,1].legend()
fig.suptitle('Candidata efectiva: cierre homogéneo y aproximación tardía',fontsize=14)
fig.savefig(RES/'cadena_efectiva.png',dpi=180);plt.close(fig)
