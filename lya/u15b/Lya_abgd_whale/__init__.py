"""Lya_abgd_whale: the Lyman-ABGD likelihood of Hooper et al. 2022 (arXiv:2206.08188, montepython/likelihoods/Lya_abgd,
unmodified) applied to WHALEFALL models through a pre-computed transfer function (docs/preregistros/ESPEC_U15b.md).

Each MCMC step runs CLASS in LCDM.  This likelihood hands the parent a proxy of that CLASS object whose linear P(k) is
P_LCDM(k) * T^2(k; s), with T^2 tabulated at the fiducial cosmology (lya/u15b/prepare_tables.py) and interpolated linearly in
the amplitude s (a nuisance parameter: s = 0 is LCDM, s = 1 the whale row), and whose spectral tilt includes d ln T^2/d ln k.
The LCDM-equivalent model the parent needs is the LCDM run itself (no second CLASS call).

Install: copy this folder into montepython/likelihoods/ (lya/u15b/run_wsl.sh does it)."""
import sys
import numpy as np
from likelihoods.Lya_abgd import Lya_abgd


class _Table:
    def __init__(self, path):
        d = np.load(path)
        self.lk = np.log(d["k_h"]); self.s = d["s"]; self.T2 = d["T2"]
        order = np.argsort(self.s); self.s = self.s[order]; self.T2 = self.T2[order]

    def T2_at(self, s):
        s = float(np.clip(s, self.s[0], self.s[-1]))
        j = int(np.clip(np.searchsorted(self.s, s) - 1, 0, len(self.s) - 2)); w = (s - self.s[j]) / (self.s[j + 1] - self.s[j])
        return (1 - w) * self.T2[j] + w * self.T2[j + 1]


class _ModelCosmo:
    """Proxy of the LCDM classy object: P_lin and the tilt include T^2(k; s); everything else is delegated."""
    def __init__(self, cosmo, lk, T2):
        self._c, self._lk, self._lT2 = cosmo, lk, np.log(np.clip(T2, 1e-300, None))
        self._dlT2 = np.gradient(self._lT2, lk)

    def _T2(self, k_Mpc):
        return float(np.exp(np.interp(np.log(k_Mpc / self._c.h()), self._lk, self._lT2)))

    def pk_lin(self, k_Mpc, z): return self._c.pk_lin(k_Mpc, z) * self._T2(k_Mpc)

    def pk_tilt(self, k_Mpc, z):
        return self._c.pk_tilt(k_Mpc, z) + float(np.interp(np.log(k_Mpc / self._c.h()), self._lk, self._dlT2))

    def __getattr__(self, name): return getattr(self._c, name)


class _Equiv:
    """The LCDM-equivalent object of the parent: the LCDM run itself (set/compute are no-ops)."""
    def __init__(self, cosmo): self._c, self.state = cosmo, False
    def struct_cleanup(self): pass
    def empty(self): pass
    def set(self, *a, **k): pass
    def compute(self, *a, **k): pass
    def h(self): return self._c.h()
    def pk_lin(self, k, z): return self._c.pk_lin(k, z)


class Lya_abgd_whale(Lya_abgd):

    def __init__(self, path, data, command_line):
        Lya_abgd.__init__(self, path, data, command_line)
        if not self.table_file:
            import io_mp
            raise io_mp.ConfigurationError("Lya_abgd_whale: set Lya_abgd_whale.table_file in the .param file")
        self.table = _Table(self.table_file)
        if self.amp_name not in data.mcmc_parameters:
            import io_mp
            raise io_mp.ConfigurationError(f"Lya_abgd_whale: the amplitude parameter '{self.amp_name}' is missing in the .param file")
        self.isLCDM = True

    def loglkl(self, cosmo, data):
        p = data.mcmc_parameters[self.amp_name]
        s = p["current"] * p["scale"]
        proxy = _ModelCosmo(cosmo, self.table.lk, self.table.T2_at(s))
        init = sys.modules.get("initialise")
        if init is None:
            import initialise as init
        orig = init.recover_cosmological_module
        init.recover_cosmological_module = lambda data: _Equiv(cosmo)
        try:
            return Lya_abgd.loglkl(self, proxy, data)
        finally:
            init.recover_cosmological_module = orig
