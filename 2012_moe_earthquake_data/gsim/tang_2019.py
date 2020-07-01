# -*- coding: utf-8 -*-

# The Hazard Library
# Copyright (C) 2012-2020, GEM Foundation
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU Affero General Public License as
# published by the Free Software Foundation, either version 3 of the
# License, or (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU Affero General Public License for more details.
#
# You should have received a copy of the GNU Affero General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.
"""
Module exports :class:`TangEtAl2019`.
"""
from __future__ import division

import numpy as np

from openquake.hazardlib.gsim.base import GMPE, CoeffsTable
from openquake.hazardlib import const
from openquake.hazardlib.imt import PGA, PGV, SA


class TangEtAl2019(GMPE):
    """
    Implements GMPE developed by Tang, Y., N. Lam, H.-H. Tsang, and E. Lumantarna (2019)
    and published as "Use of macroseismic intensity data to validate a regionally adjustable 
    ground motion prediction model" (2019, Geosciences 9, 422, doi: 10.3390/geosciences9100422.
    """
    
    #: Define for tectonic region type stable shallow crust. 
    DEFINED_FOR_TECTONIC_REGION_TYPE = const.TRT.STABLE_CONTINENTAL

    #: Supported intensity measure types are spectral acceleration,
    #: peak ground velocity and peak ground acceleration, see pag 1990.
    DEFINED_FOR_INTENSITY_MEASURE_TYPES = set([
           PGV,
           SA
    ])

    #: Supported intensity measure component is orientation-independent
    #: measure :attr:`~openquake.hazardlib.const.IMC.GMRotI50`.
    DEFINED_FOR_INTENSITY_MEASURE_COMPONENT = const.IMC.AVERAGE_HORIZONTAL
    
    #: Defined for standard deviation types inter-event, intra-event
    #: and total.
    DEFINED_FOR_STANDARD_DEVIATION_TYPES = set([
            const.StdDev.TOTAL,
            const.StdDev.INTER_EVENT,
            const.StdDev.INTRA_EVENT])    

    #: Required site parameters is Vs30.
    #: See paragraph 'Functional Form of the Generic GMPE', pag 1992.
    REQUIRES_SITES_PARAMETERS = set(('vs30', ))
   
    #: Required rupture parameters are magnitude.
    REQUIRES_RUPTURE_PARAMETERS = set(('mag', 'hypo_depth', 'stress_drop', ))

    #: Required distance measure is rrup.
    #: See paragraph 'Functional Form of the Generic GMPE', pag 1991.
    REQUIRES_DISTANCES = set(('rrup', ))

    def get_mean_and_stddevs(self, sites, rup, dists, imt, stddev_types):
        """
        See :meth:`superclass method
        <.base.GroundShakingIntensityModel.get_mean_and_stddevs>`
        for spec of input and result values.
        """
        # extracting dictionary of coefficients specific to required
        # intensity measure type.
        C = self.COEFFS[imt]
        #C_PGV = self.COEFFS_PGV[imt]  # need to fix PGV coeffs
        
        # set period-independent parameters
        stress_drop = 200.
        Q0 = 200.
        kappa = 0.0328
        V8 = 3.5
        Vs30 = sites.vs30 / 1000.
        
        # set distances
        hff = 10 ** (-0.405 + 0.235 * rup.mag)
        R = np.zeros_like(dists.rrup)
        R = np.sqrt(dists.rrup**2 + (hff**2))
        
        # ln mean of either PGV or 5%-damped response spectra
        if imt == PGV():
            mean = np.log(C_PGV['d0']) + self._compute_geometric_factor(R, hff) + \
                   self._compute_generic_source_factor(rup, C_PGV) + \
                   self._compute_anelastic_attenuation_factor(rup, R, C_PGV) + \
                   self._compute_PGV_crustal_modification_factor(rup, Vs30, kappa, V8, C_PGV)
        else:
            mean = np.log(C['d0']) + self._compute_geometric_factor(R, hff) + \
                   self._compute_generic_source_factor(rup, stress_drop, C) + \
                   self._compute_anelastic_attenuation_factor(rup, R, Q0, C) + \
                   self._compute_crustal_modification_factor(rup, Vs30, kappa, stress_drop, V8, C)
            
            stddevs = self._get_stddevs(C, stddev_types, num_sites=len(sites.vs30))
    
            return mean, stddevs
    
    def _get_stddevs(self, C, stddev_types, num_sites):
        """
        No Sigma model given - assume Atkinson & Adams (2013) sigma
        """
        stddevs = []
        for stddev_type in stddev_types:
            stddevs.append(np.ones(num_sites) * 0.63)
        return stddevs

    def _compute_geometric_factor(self, R, hff):
        """
        Compute distance scaling term, equations - not given in paper - from spreadsheet provided by authors
        """
        # Computing geometrical spreading, Fz.
        b1 = -1.0
        b2 = 0.0
        b3 = -0.5
        Rt1 = 70    # Transition distance of 70km.
        Rt2 = 130
                
        G = np.empty_like(R)
        idx = np.where(R<=Rt1)[0]
        G[idx] = R[idx] ** b1
        
        idx = np.where((R>Rt1) & (R<=Rt2))[0]
        G[idx] = (Rt1**b1) * ((R[idx]/Rt1)**b2)
        
        idx = np.where(R>Rt2)[0]
        G[idx] = (Rt1**b1) * ((Rt2/Rt1)**b2) * ((R[idx]/Rt2)**b3)
        
        G30 = 1. / np.sqrt(30.**2 + hff**2)
        
        G_1 = G / G30
        
        #print('R', R, 'G', G, 'G_1', G_1)
        
        return np.log(G_1)
    
    def _compute_generic_source_factor(self, rup, stress_drop, C):
        """
        Compute magnitude-scaling term, equation (3), pag 4/22.
        """
        
        return np.log(10**(C['a1'] * (rup.mag**C['a2']) * (stress_drop**C['a3']) + C['a4']))
        
    def _compute_anelastic_attenuation_factor(self, rup, R, Q0, C):
        """
        Compute distance-scaling terms, equations (4-5), pag 4-5/22.
        """
        
        beta = np.log(10**((C['b1'] * rup.mag + C['b2']) * (Q0**C['b3']) * (np.log10(R)**C['b4']) + C['b5']))
        
        
        
        residual = np.log(10**(C['c1'] * np.log10(R)**4 + C['c2'] * np.log10(R)**3 \
                   + C['c3'] * np.log10(R)**2 + C['c4'] * np.log10(R) + C['c5']))
        
        return beta + residual
        
    def _compute_PGV_crustal_modification_factor(self, rup, Vs30, kappa, V8, C):
        """
        Compute crustal-scaling terms for PGV, equations for PGV (6-9), pag 5-7/22.
        """
        
        # upper crustal amplification, eqn 7
        upperampli = np.log(10**(C['r1'] * rup.mag**C['r2'] * Vs30**C['r3'] + C['r4']*Vs30))
        
        # kappa attenuation, eqn 8
        upperatten = np.log(10**(C['r5'] * (rup.mag**C['r6']) * (kappa**C['r7']) + C['r8'])) 
        
        # adjustment
        ampli_adjustment = np.log((3.8 / V8)**(-0.273 * rup.mag + 3.278))
        
        crustal_modification = upperampli + upperatten + ampli_adjustment
        
        return crustal_modification
    
    
    def _compute_crustal_modification_factor(self, rup, Vs30, kappa, stress_drop, V8, C):
        """
        Compute crustal-scaling terms for response spectral ordinates, 
        as published in Tang et al. (2019) - Other*** 
        """
        
        # upper crustal amplification, eqn 7
        upperampli = np.log(10**(C['r1'] * Vs30**C['r2'] + C['r3'])) 
        
        # kappa attenuation, eqn 8
        upperatten = np.log(10**(C['r4'] * (rup.mag**C['r5']) * (kappa**C['r6']) + C['r7'] * kappa + C['r8']))
        
        # adjustment
        ampli_adjustment = np.log((3.8 / V8)**(C['r9'] + C['r10'] * rup.mag \
                           + C['r11'] * stress_drop \
                           + C['r12'] * rup.mag**2 \
                           + C['r13'] * rup.mag * stress_drop \
                           + C['r14'] * stress_drop**2))
        
        crustal_modification = upperampli + upperatten + ampli_adjustment
        
        return crustal_modification

    # Coefficient table from Tang (2019) PhD thesis

    COEFFS = CoeffsTable(sa_damping=5, table="""\
            imt	d0	a1	a2	a3	a4	b1	b2	b3	b4	b5	c1	c2	c3	c4	c5	r1	r2	r3	r4	r5	r6	r7	r8	r9	r10	r11	r12	r13	r14
            0.01	0.58869	0.658224065	0.5136	0.1852	-4.416	0.03058	-0.5123	-0.2147	3.286	0.3401	0.2288	-1.253	2.192	-1.559	0.4035	-19.15	0.02272	19.68	-7.361	-0.1556	0.1813	0.6193	1.432	1.031	-0.01042	3.27E-05	0.0008653	-4.52E-06	-8.06E-09
            0.02	0.48483	0.637628689	0.5235	0.1874	-4.399	0.02061	-0.3296	-0.2425	3.824	0.2256	0.2989	-1.861	3.946	-3.451	1.045	2.876	-0.1352	-2.414	-10.24	-0.209	0.3884	4.425	0.4725	1.157	-0.05362	1.66E-04	0.004503	-2.37E-05	-2.70E-08
            0.03	0.41744	0.727583324	0.5093	0.1747	-4.575	0.01609	-0.2534	-0.2673	4.183	0.1828	0.2599	-1.661	3.643	-3.302	1.02	1.679	-0.2259	-1.253	-16.02	-0.186	0.5656	10.2	0.2491	1.277	-0.09357	0.0002952	0.007778	-4.09E-05	-7.02E-08
            0.05	0.33436	0.833040414	0.4958	0.1635	-4.804	0.01195	-0.1856	-0.304	4.629	0.1324	0.1421	-0.9229	2.068	-1.941	0.6394	1.221	-0.2949	-0.8125	-83.48	-0.04534	0.9081	77.79	0.09447	1.753	-0.254	0.0008051	0.02111	-1.12E-04	-1.94E-07
            0.1	0.23292	2.708559456	0.3409	0.07316	-7.343	0.009409	-0.1505	-0.3894	5.13	0.09245	0.01014	-0.0614	0.1328	-0.1766	0.1003	0.764	-0.4258	-0.4021	9.581	0.3886	1.65	-15.99	0.01376	4.799	-1.32	0.004104	0.1127	-5.82E-04	-9.78E-07
            0.14	0.18979	4.833944	0.2648	0.04477	-9.833	0.00945	-0.1497	-0.4372	5.296	0.07129	-0.01091	0.08805	-0.233	0.1876	-0.02618	0.489	-0.6369	-0.1646	3.372	0.6984	1.943	-10.41	0.005373	6.825	-1.976	0.005645	0.1653	-0.0007743	-1.59E-06
            0.18	0.16010	6.900794503	0.2218	0.03205	-12.16	0.01146	-0.1801	-0.4902	5.337	0.06476	-0.00699	0.06147	-0.1716	0.1366	-0.01617	0.5184	-0.5764	-0.1989	-169	-0.00536	0.9927	162.7	0.01841	8.172	-2.342	0.005791	0.1893	-0.0007675	-1.68E-06
            0.22	0.13793	10.5686693	0.1742	0.02096	-16.11	0.01195	-0.1943	-0.5262	5.369	0.05658	-0.00713	0.06448	-0.1846	0.1569	-0.03846	0.4157	-0.6863	-0.1173	-291.4	-0.00025	0.9974	285.7	0.01045	8.878	-2.501	0.006269	0.197	-0.0008172	-1.76E-06
            0.26	0.12051	7.122460037	0.2306	0.02865	-12.56	0.01338	-0.21	-0.557	5.417	0.04402	-0.01275	0.1053	-0.2889	0.2671	-0.0637	0.3586	-0.769	-0.07623	-4.092	-0.3354	0.8287	-2.062	0.01251	9.38	-2.603	0.006935	0.2009	-0.000891	-1.97E-06
            0.3	0.10635	6.999140106	0.2399	0.02778	-12.44	0.01498	-0.2274	-0.5788	5.419	0.04241	-0.00679	0.05906	-0.1672	0.1427	-0.01615	0.3446	-0.8133	-0.073	-2.992	-0.4932	0.7827	-2.713	0.005529	9.206	-2.52	0.008759	0.1918	-0.001066	-3.58E-06
            0.4	0.08016	23.00215292	0.11	0.008015	-29.25	0.01951	-0.2755	-0.6289	5.476	0.0447	0.004809	-0.00598	-0.0469	0.06567	-0.06618	0.3081	-0.7762	-0.08746	-8.008	-2.094	0.4837	-3.108	0.01393	7.916	-1.853	0.007233	0.1198	-0.0007119	-4.55E-06
            0.5	0.06219	20.53736333	0.1262	0.007942	-26.94	0.02019	-0.2825	-0.6617	5.556	0.03494	-0.00317	0.05364	-0.1993	0.222	-0.08987	0.2685	-0.8711	-0.0515	-17.21	-2.645	0.4594	-2.603	0.008318	5.509	-0.79	0.004174	0.01307	-9.27E-05	-5.84E-06
            0.6	0.04928	57.9181757	0.05503	0.002595	-64.92	0.02336	-0.3204	-0.6991	5.585	0.03179	0.00436	-0.0024	-0.05127	0.05862	-0.0224	0.2036	-1.07	0.002703	-30.09	-2.984	0.46	-2.216	0.007437	3.823	-0.05803	0.002456	-0.0589	0.0002543	-6.34E-06
            0.7	0.03972	23.99977608	0.118	0.005657	-30.62	0.02829	-0.3764	-0.7411	5.639	0.02613	0.00456	-0.00089	-0.06358	0.08136	-0.00901	0.1646	-1.254	0.03491	-30.09	-2.984	0.46	-2.216	0.007437	2.366	0.562	0.0009248	-0.1189	0.0005519	-6.60E-06
            0.8	0.03250	10.66591622	0.2202	0.0114	-16.9	0.03067	-0.3912	-0.7596	5.696	0.02231	0.0116	-0.04802	0.04639	-0.02111	0.0115	0.1434	-1.366	0.02643	-393.8	-4.77	0.4709	-1.825	0.002573	1.1	1.116	-0.0009826	-0.1706	0.0007256	-4.41E-06
            0.9	0.02697	24.53169812	0.1214	0.004623	-31.4	0.03059	-0.3897	-0.7738	5.741	0.01836	0.01199	-0.05067	0.05426	-0.03273	0.004222	0.1474	-1.303	0.006787	-349.4	-4.629	0.4895	-1.636	-0.00182	0.09228	1.502	-0.001583	-0.2048	0.0008585	-4.71E-06
            1	0.02265	25.96939296	0.1187	0.004031	-32.96	0.02884	-0.3635	-0.7827	5.829	0.009252	0.01524	-0.07654	0.1286	-0.123	0.05115	0.1215	-1.516	0.0204	-115	-3.738	0.4838	-1.415	0.00605	-0.2163	1.6	-1.67E-03	-0.2085	0.0008098	-4.14E-06
            1.5	0.01105	21.50263857	0.1436	0.003456	-28.46	0.04115	-0.4731	-0.8733	6.073	0	0.02104	-0.1085	0.1937	-0.1903	0.1101	0.1703	-0.8065	-0.02665	-483.5	-4.57	0.5104	-1.038	0.001712	-0.7349	1.729	-0.002392	-0.2058	0.000896	-4.35E-06
            2	0.00650	19.11966288	0.1598	0.003136	-26.04	0.0772	-0.7933	-0.9901	6.301	-0.00818	0.0151	-0.07618	0.1392	-0.1644	0.1386	0.1276	-1.14	0.02052	-254.1	-3.952	0.53	-0.7405	0.00877	-0.1029	1.412	-0.002179	-0.1638	0.0007398	-3.22E-06
            3	0.00308	26.56340952	0.1303	0.001744	-34.08	0.1891	-1.713	-1.196	6.859	0.004653	0.00331	0.006942	-0.05006	-0.02251	0.1028	0.1062	-1.356	0.02304	-484.2	-4.315	0.5722	-0.585	0.003078	1.433	0.6956	-0.000875	-0.07835	0.0002907	-1.16E-06
            4	0.00182	23.9546852	0.1432	0.001845	-31.36	0.4871	-4.202	-1.402	7.273	0.01409	-0.00496	0.06972	-0.2007	0.09165	0.1029	0.1039	-1.325	0.01675	-1144	-4.897	0.577	-0.5958	-0.00041	2.011	0.4239	-0.0005047	-0.04617	0.0001649	-6.85E-07
            5	0.00121	9.007683802	0.2892	0.004663	-15.58	2.271	-18.47	-1.597	7.171	0.02252	-0.01405	0.1482	-0.4185	0.308	0.03694	0.109	-1.223	0.01195	-2116	-5.348	0.609	-0.6395	-0.01154	2.618	0.1655	-0.0001854	-0.0184	0.00006758	-3.16E-07
            6	0.00086	4.350575258	0.4345	0.01007	-10.08	2.423	-19.3	-1.504	6.695	0.02346	-0.02284	0.2388	-0.6922	0.5819	-0.0254	0.1084	-1.245	0.01387	-34470	-7.403	0.6309	-0.811	-0.01989	2.785	0.09342	-0.0001036	-0.01046	4.04E-05	-2.05E-07
            8	0.00051	2.060528117	0.5977	0.02473	-6.954	0.8278	-6.66	-1.226	6.161	0.02116	-0.02469	0.2708	-0.8048	0.6959	-0.04405	0.1089	-1.336	0.01741	-60650	-7.77	0.6697	-1.007	-0.02378	2.694	0.1175	-0.0001384	-0.01155	4.01E-05	-1.75E-07
    """)

    # Coefficient for PGV equation
    
    COEFFS_PGV = CoeffsTable(sa_damping=5, table="""\
            imt	d0	a1	a2	a3	a4	b1	b2	b3	b4	b5	c1	c2	c3	c4	c5	r1	r2	r3	r4	r5	r6	r7	r8
            PGV	39	27.79733	0.08414	0.0059	-33.35	0.06287	-0.6326	-0.4963	4.431	0.06135	0.01714	-0.06931	0.08404	-0.09224	0.0876	0.7334	-0.5251	-0.8479	-0.019	-21.35	-1.351	0.5584	-0.03336
    """)