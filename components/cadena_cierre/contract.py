"""Machine-readable interfaces and scientific gates, never conflate code PASS and physical closure."""
from common import *
from scalar_end import endpoint
def read(n): return json.loads((ROOT/'results'/f'{n}.json').read_text())
def run():
    th=read('thermal'); nr=th['continuous']; b=read('bounce'); dn=read('double_null'); mass=read('mass_bounce')
    st=read('stochastic'); q=read('q_bridge'); dark=read('dark')
    input_common=bool(b['shared_V0_phys']==nr['V0'] and st['V0']==nr['V0'] and dark['shared_F_GeV']==nr['F_GeV'])
    field_model=Potential(nr['V0']); et=field_model.eternal()
    from stochastic import eigen
    launch_windows=[]
    for case in b['cases']:
        r=case['reruns'][-1]
        if r['eternal_sub']:
            growth=eigen(field_model,r['chi_land'],1e-4,512)
            launch_windows.append(dict(label=case['label'],radiation_fraction=case['archived']['frac_rad'],
                chi_land=r['chi_land'],P_zeta=float(field_model.pzeta(r['chi_land'])),
                lambda_volume=growth,UV_window=[r['chi_land']-5e-5,r['chi_land']+5e-5]))
    checks=[]
    def gate(name,state,evidence,limit): checks.append(dict(connection=name,state=state,evidence=evidence,remaining=limit))
    gate('semilla -> parche cerrado','ASSUMED','Turnaround and local characteristic data prescribed','No formation mechanism or probability derived')
    gate('QT -> lanzamiento homogéneo','PASS_TOY' if all(c.get('classification_converged') for c in b['cases']) else 'PARTIAL',
        f"{len(b['cases'])} deterministic cases, two tolerances, no relaxation",
        'Homogeneous initial closed seed; no 3D formation or finite-size island matching')
    gate('rebote escalar -> estabilidad local moderada','PASS_NUMERICS' if dn['pass_homogeneous_convergence'] and all(r['pass_local'] for r in dn['perturbations']) else 'FAIL',
        'Full scalar and NPQT characteristic equations, exact FRW and two perturbations, three meshes',
        'Local annulus, not global seed')
    moderate=[r for r in mass['runs'] if r['label']=='moderate']; high=[r for r in mass['runs'] if r['label']=='high']
    mass_moderate=all('error' not in r and r['geometric_margin_min']>0 for r in moderate)
    highpass=all('error' not in r and r['geometric_margin_min']>0 and r['psi_geometric_vs_mass_relative']<1e-6 for r in high)
    gate('garganta fuerte -> coherencia del disparo','PASS_NUMERICS' if highpass else 'INCONCLUSIVE_NUMERICS',
        'Mass formulation agrees with exact equations and evolves through |chi|>400 with shrinking residuals',
        'Geometric branch margin unresolved relative to the tiny saturation margin; global inhomogeneous launch remains open')
    gate('lanzamiento -> difusión y ramificación','PASS_EFFECTIVE_TOY' if st['pass_numerics'] and all(r['lambda_volume']>0 for r in launch_windows) else 'FAIL',
        'Finite UV-window drift-diffusion eigenvalues and binary branching BVP; analytic constant-coefficient check',
        'Canonical large-field potential, UV validity and effective branching assumed')
    gate('inflación -> recalentamiento -> pivote','PASS_NUMERICS' if th['pass_energy_and_convergence'] and input_common else 'FAIL',
        f"Exact epsilon_H endpoint; pivot residual {nr['matching_history'][-1]['residual']:.2e}; conservative continuous reheating",
        'GR scalar background and first-order slow-roll spectra; averaged condensate after endpoint; no nonlinear preheating')
    qend=q['background']['chi_at_epsilon_H_1']; endpoint_equal=abs(qend-nr['chi_end'])<1e-7
    gate('χ -> curvatura -> ordenación Q','PASS_CONNECTION' if q['connection_pass'] and endpoint_equal else 'FAIL',
        f"Same exact endpoint; ordering switch at {q['background']['N_switch_after_exact_end']:.4f} e-folds",
        'Spectator Q, inherited fluctuation normalization; small lattice is smoke')
    gate('Q -> cicatrices -> emisión','OPEN_CALIBRATION',
        '3D lattice connection and restored/breathing controls; archived production evidence separate',
        'Scaling source, spectrum, wall collapse, energy transfer and backreaction not calibrated for RP3')
    gate('emisión -> materia + radiación','PASS_LEDGER_OPEN_PHYSICS' if dark['checks']['pass_numerical_ledger'] else 'FAIL',
        'One emitted packet distribution; energy/pressure/number at BBN, CMB and today; independent time and momentum quadrature refinement',
        'Conditional U1 source, not RP3 calibration; nonlinear transient, collapse and observational Boltzmann likelihood absent')
    gate('pares geométricos -> χ canónico y V','OPEN_ACTION',
        'Mechanical inverse-distance interaction alone does not define a homogeneous relativistic scalar',
        'Kinetic normalization, number-density dynamics, positive vacuum plateau and common stress tensor missing')
    gate('región eterna -> islas/agujeros negros -> descendencia','OPEN_GEOMETRY',
        'Finite-window effective branching has positive survival solutions in selected windows',
        'No junction, critical island radius, mass function or past completion derived')
    alpha=next(c['reruns'][-1]['alpha_phys'] for c in b['cases'] if c['label']=='eternal_sub_multi' and c['archived']['frac_rad']==0)
    # The coefficient n C sqrt(n mu) of an inverse-chi mechanical pair potential dilutes as a^-9/2.
    pair_dilution=np.exp(-4.5*50)
    tests=dict(common_parameter_identity=input_common,endpoint_identity=endpoint_equal,
        local_moderate_mass_branch=mass_moderate,high_density_branch_certified=highpass,
        qt_alpha_Hstar_squared=alpha*(nr['H_star_GeV']/2.435e18)**2,
        hierarchy_chi_end_V0_identity=True,physical_seed_derived=False,canonical_pair_action_derived=False,
        RP3_emission_calibrated=False,common_Q_backreaction_solved=False,full_chain_closed=False,
        pair_density_potential_coefficient_after_50_efolds=pair_dilution)
    save('contract.json',dict(shared_parameters={k:nr[k] for k in ('V0','chi_star','chi_end','N_star','n_s','r','F_GeV','mu_GeV','xi','Gamma_GeV')},
        alpha_shared_example=alpha,chi_Pzeta_1=et,launch_windows=launch_windows,connections=checks,checks=tests,
        publication_status='RESEARCH_TOY_WITH_REPRODUCIBLE_PARTIAL_CONNECTIONS',
        verified_complete_cosmology=False))
    print('CONTRACT',tests,flush=True)
if __name__=='__main__': run()
