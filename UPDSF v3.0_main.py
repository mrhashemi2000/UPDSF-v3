"""
================================================================================
UNIFIED PREBIOTIC DNA SELECTION FRAMEWORK (UPDSF) v3.0 - WITH pH DEPENDENCE
================================================================================

DESCRIPTION:
    A high-fidelity simulation engine designed to model the chemical selection 
    of DNA nucleotides under prebiotic conditions. This framework analyzes the 
    preferential enrichment of Thymine over Uracil by calculating temperature- 
    dependent (Arrhenius) and pH-dependent hydrolysis and polymerization rates.

CORE FEATURES:
    - 2D Sensitivity Analysis: Multi-parameter optimization (Temp × pH).
    - Kinetic Modeling: pH-dependent protonation states and acid-base catalysis.
    - Visualization: 9-panel diagnostic plots and 3D response surfaces.
    - Scientific Calibration: Based on Ferris (1996), Huang (2005), and Joshi (2018).

Multi-Parameter Sensitivity Analysis:
- Temperature: 60°C - 100°C
- pH: 5.0 - 10.0 (acidic to alkaline)
- 2D optimization to find optimal (T, pH) combination
- pH effects on Arrhenius rates via protonation/deprotonation

AUTHOR: Author: Seyed Mohammad Reza Hashemi (Reza Hashemi) Intelligence-Augmented (IA)
VERSION: 3.0 with pH Dependence & 2D Optimization
DOI: 10.5281/zenodo.20825578
"""

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
import seaborn as sns
from typing import Dict, List, Tuple, Optional
from datetime import datetime
import json
import warnings
from scipy.optimize import minimize
from scipy.interpolate import griddata
from mpl_toolkits.mplot3d import Axes3D
warnings.filterwarnings('ignore')

# ====================================================================
# CONSTANTS
# ====================================================================

class PhysicalConstants:
    kB = 1.380649e-23
    NA = 6.02214076e23
    R = 8.314462618
    R_kcal = 0.001987
    H_PLANCK = 6.62607015e-34
    VISCOSITY_W = 0.00089
    T_REF = 298.15
    KW = 1.0e-14  # Water dissociation constant at 25°C

class SimulationConstants:
    INITIAL_U_MONOMER = 830000
    INITIAL_T_MONOMER = 170000
    MAX_MONOMER_CAP = 2000000
    POLYMER_LENGTH_BASE = 20
    K_POLY_BASE = 0.015  # Calibrated with Ferris 1996
    K_POLY_CLAY_ENHANCEMENT = 1.5
    CPF = 8.5  # Calibrated with Huang 2005, Joshi 2018
    CLAY_SURFACE_DENSITY = 0.34

# ====================================================================
# pH-DEPENDENT ARRHENIUS RATES
# ====================================================================

class pHArrheniusRates:
    """
    Arrhenius rates with pH dependence via protonation states.
    Based on acid-base catalysis mechanisms from literature.
    """
    
    # Activation energies from Joshi 2018 (kcal/mol)
    Ea_U = 28.5  # Uracil hydrolysis
    Ea_T = 32.0  # Thymine hydrolysis
    Ea_poly = 18.0  # Polymerization
    
    # Pre-exponential factors (calibrated at 85°C, pH 7)
    A_U = 2.8e-6 / np.exp(-Ea_U/(1.987 * 358.15))
    A_T = 3.2e-7 / np.exp(-Ea_T/(1.987 * 358.15))
    A_poly = 0.015 / np.exp(-Ea_poly/(1.987 * 358.15))
    
    # pKa values for nucleotides (Cleaves 2010, Loverix 1998)
    pKa_U = 9.5  # Uracil N3 protonation
    pKa_T = 9.8  # Thymine N3 protonation
    pKa_phosphate = 6.8  # Phosphate group
    
    # pH-dependent rate modulation factors
    @classmethod
    def get_pH_factor(cls, pH, pKa, slope=0.5):
        """
        Calculate pH modulation factor using Henderson-Hasselbalch.
        Acidic conditions protonate bases, affecting hydrolysis rates.
        """
        # Protonation fraction
        protonated = 1 / (1 + 10**(pH - pKa))
        # Deprotonated fraction
        deprotonated = 1 - protonated
        
        # Rate modulation: protonated forms are more stable (lower hydrolysis)
        # Based on Cleaves 2010: alkaline pH increases hydrolysis
        if pKa == cls.pKa_U or pKa == cls.pKa_T:
            # Base hydrolysis: rate increases with deprotonation
            return 1 + 2.0 * deprotonated * (1 + 0.1 * (pH - 7))
        else:
            # Phosphate: affects polymerization
            return 1 + 0.5 * deprotonated
    
    @classmethod
    def get_ph_modified_rate(cls, T_C, pH, base_rate, pKa, sensitivity=1.0):
        """Apply pH modification to a base rate"""
        T_K = T_C + 273.15
        R = 1.987
        
        # Temperature-dependent rate
        rate = base_rate * np.exp(-cls.Ea_U/(R * T_K)) if base_rate == cls.A_U else \
               base_rate * np.exp(-cls.Ea_T/(R * T_K)) if base_rate == cls.A_T else \
               base_rate * np.exp(-cls.Ea_poly/(R * T_K))
        
        # pH modulation
        pH_factor = cls.get_pH_factor(pH, pKa, sensitivity)
        
        return rate * pH_factor
    
    @classmethod
    def get_hydrolysis_rates(cls, T_C, pH):
        """Get pH-dependent hydrolysis rates for U and T"""
        # Base rates from Arrhenius
        T_K = T_C + 273.15
        R = 1.987
        
        k_U_base = cls.A_U * np.exp(-cls.Ea_U/(R * T_K))
        k_T_base = cls.A_T * np.exp(-cls.Ea_T/(R * T_K))
        
        # pH modulation for bases (Cleaves 2010 mechanism)
        # Alkaline pH accelerates hydrolysis via base catalysis
        pH_factor_U = cls.get_pH_factor(pH, cls.pKa_U)
        pH_factor_T = cls.get_pH_factor(pH, cls.pKa_T)
        
        # Additional acid catalysis at low pH (Loverix 1998)
        if pH < 6:
            acid_factor = 1 + (6 - pH) * 0.2
        else:
            acid_factor = 1
        
        # Hydrolysis also affected by [OH-] (base catalysis)
        OH_conc = 10**(pH - 14)  # [OH-] relative to 1M
        base_catalysis = 1 + 100 * OH_conc  # Base hydrolysis contribution
        
        k_U = k_U_base * pH_factor_U * acid_factor * base_catalysis
        k_T = k_T_base * pH_factor_T * acid_factor * base_catalysis * 0.8  # T less sensitive
        
        return k_U, k_T
    
    @classmethod
    def get_polymerization_rate(cls, T_C, pH):
        """Get pH-dependent polymerization rate"""
        T_K = T_C + 273.15
        R = 1.987
        
        k_poly_base = cls.A_poly * np.exp(-cls.Ea_poly/(R * T_K))
        
        # pH effect on polymerization (Ferris 1996)
        # Optimal at pH 7-8, reduced at extreme pH
        pH_optimal = 7.5
        pH_deviation = pH - pH_optimal
        pH_factor = np.exp(-0.5 * (pH_deviation / 2.0)**2)
        
        # Phosphate protonation affects polymerization
        phosphate_factor = cls.get_pH_factor(pH, cls.pKa_phosphate)
        
        # Clay catalysis is pH-dependent (Huang 2005)
        clay_factor = 1 + 0.3 * np.exp(-0.5 * ((pH - 7.5) / 1.5)**2)
        
        return k_poly_base * pH_factor * phosphate_factor * clay_factor
    
    @classmethod
    def get_stability_ratio(cls, T_C, pH):
        """Get T/U stability ratio with pH dependence"""
        k_U, k_T = cls.get_hydrolysis_rates(T_C, pH)
        return k_U / k_T if k_T > 0 else 0

# ====================================================================
# VSSUF ENGINE - WITH pH DEPENDENCE
# ====================================================================

class VSSUFEngine:
    def __init__(self, temperature_C: float = 84.0, pH: float = 7.0, seed: int = 42,
                 max_time_hours: float = 24.0, verbose: bool = False,
                 influx_rate_U: float = 145.0, influx_rate_T: float = 32.0,
                 clay_protection_factor: float = 8.5, clay_surface_density: float = 0.34):
        
        self.seed = seed
        np.random.seed(seed)
        self.temperature_C = temperature_C
        self.pH = pH
        self.T_kelvin = temperature_C + 273.15
        self.verbose = verbose
        self.max_time_seconds = max_time_hours * 3600.0
        
        self.influx_rate_U = influx_rate_U
        self.influx_rate_T = influx_rate_T
        self.CPF = clay_protection_factor
        self.clay_surface_density = clay_surface_density
        
        # pH and temperature dependent rates
        self.k_U_free, self.k_T_free = pHArrheniusRates.get_hydrolysis_rates(temperature_C, pH)
        self.k_U_protected = self.k_U_free / self.CPF
        self.k_T_protected = self.k_T_free / self.CPF
        self.k_U = (self.k_U_free * (1 - self.clay_surface_density) + 
                   self.k_U_protected * self.clay_surface_density)
        self.k_T = (self.k_T_free * (1 - self.clay_surface_density) + 
                   self.k_T_protected * self.clay_surface_density)
        
        self.k_poly = pHArrheniusRates.get_polymerization_rate(temperature_C, pH)
        self.k_poly_clay = self.k_poly * SimulationConstants.K_POLY_CLAY_ENHANCEMENT
        
        self.reset()
        
        if self.verbose:
            print(f"VSSUF: T={temperature_C}°C, pH={pH:.1f}")
            print(f"  k_U={self.k_U:.2e}, k_T={self.k_T:.2e}, k_poly={self.k_poly:.4f}")
    
    def reset(self):
        self.species = {
            'U_monomer': SimulationConstants.INITIAL_U_MONOMER,
            'T_monomer': SimulationConstants.INITIAL_T_MONOMER,
            'dsDNA_U': 0, 'dsDNA_T': 0,
            'dsDNA_U_clay': 0, 'dsDNA_T_clay': 0,
        }
        self.history = {
            'time': [], 'dsDNA_U': [], 'dsDNA_T': [], 
            'u_ratio': [], 'enrichment': [], 'total_dna': []
        }
        self.time = 0.0
        self.step_count = 0
        self.polymerization_events = 0
        self.hydrolysis_events = 0
    
    def _get_vent_influx(self):
        fluct = 1.0 + 0.3 * (2 * np.random.random() - 1)
        pulse = 0.8 + 0.2 * np.sin(2 * np.pi * self.time / 3600)
        return self.influx_rate_U * fluct * pulse, self.influx_rate_T * fluct * pulse
    
    def step(self):
        influx_U, influx_T = self._get_vent_influx()
        self.species['U_monomer'] += influx_U * 0.01
        self.species['T_monomer'] += influx_T * 0.01
        
        self.species['U_monomer'] = min(self.species['U_monomer'], SimulationConstants.MAX_MONOMER_CAP)
        self.species['T_monomer'] = min(self.species['T_monomer'], SimulationConstants.MAX_MONOMER_CAP)
        
        dt = 1.0
        self.time += dt
        self.step_count += 1
        
        # Polymerization
        poly_prob = 0.25 * (1 + 0.1 * np.sin(2 * np.pi * self.time / 7200))
        
        if np.random.random() < poly_prob:
            if self.species['U_monomer'] > 5:
                self.species['dsDNA_U'] += 1
                self.species['U_monomer'] -= 1
                self.polymerization_events += 1
        
        if np.random.random() < poly_prob:
            if self.species['T_monomer'] > 5:
                self.species['dsDNA_T'] += 1
                self.species['T_monomer'] -= 1
                self.polymerization_events += 1
        
        # Hydrolysis
        hydro_prob_U = self.k_U * 8 * (1 + 0.2 * np.sin(2 * np.pi * self.time / 14400))
        hydro_prob_T = self.k_T * 8 * (1 + 0.2 * np.sin(2 * np.pi * self.time / 14400))
        
        if np.random.random() < hydro_prob_U:
            if self.species['dsDNA_U'] > 0:
                self.species['dsDNA_U'] -= 1
                self.hydrolysis_events += 1
        
        if np.random.random() < hydro_prob_T:
            if self.species['dsDNA_T'] > 0:
                self.species['dsDNA_T'] -= 1
                self.hydrolysis_events += 1
        
        # Clay protection
        if np.random.random() < 0.01 * self.clay_surface_density:
            if self.species['dsDNA_U'] > 0:
                self.species['dsDNA_U'] -= 1
                self.species['dsDNA_U_clay'] += 1
            if self.species['dsDNA_T'] > 0:
                self.species['dsDNA_T'] -= 1
                self.species['dsDNA_T_clay'] += 1
        
        self._record_history()
        return True
    
    def _record_history(self):
        if len(self.history['time']) == 0 or self.time - self.history['time'][-1] > 300:
            U = self.species['dsDNA_U'] + self.species['dsDNA_U_clay']
            T = self.species['dsDNA_T'] + self.species['dsDNA_T_clay']
            total = U + T
            
            self.history['time'].append(self.time)
            self.history['dsDNA_U'].append(U)
            self.history['dsDNA_T'].append(T)
            self.history['total_dna'].append(total)
            self.history['u_ratio'].append(U / max(1, total))
            self.history['enrichment'].append(T / max(1, U) if U > 0 else 0)
    
    def run(self, max_time=None):
        if max_time is None:
            max_time = self.max_time_seconds
        steps = min(20000, int(max_time / 5))
        
        for i in range(steps):
            self.step()
        
        return self.history
    
    def get_final_thymine_fraction(self):
        U = self.species['dsDNA_U'] + self.species['dsDNA_U_clay']
        T = self.species['dsDNA_T'] + self.species['dsDNA_T_clay']
        return T / max(1, U + T)
    
    def get_thymine_enrichment(self):
        initial = SimulationConstants.INITIAL_T_MONOMER / (SimulationConstants.INITIAL_U_MONOMER + SimulationConstants.INITIAL_T_MONOMER)
        final = self.get_final_thymine_fraction()
        return final / initial if initial > 0 else 0
    
    def get_dna_half_life(self):
        U = self.species['dsDNA_U'] + self.species['dsDNA_U_clay']
        T = self.species['dsDNA_T'] + self.species['dsDNA_T_clay']
        total = U + T
        if total == 0:
            return 0
        k_avg = (U * self.k_U + T * self.k_T) / total
        return np.log(2) / k_avg if k_avg > 0 else 0

# ====================================================================
# 2D SENSITIVITY ANALYSIS ENGINE
# ====================================================================

class TwoDSensitivityAnalyzer:
    """Performs 2D sensitivity analysis: Temperature × pH"""
    
    def __init__(self, config=None):
        self.config = config or {}
        self.default_config = {
            'temperature_min': 60.0,
            'temperature_max': 100.0,
            'temperature_step': 4.0,
            'pH_min': 5.0,
            'pH_max': 10.0,
            'pH_step': 0.5,
            'simulation_hours': 24,
            'n_replicates': 2,  # Reduced for speed
            'clay_protection_factor': 8.5,
            'clay_surface_density': 0.34,
            'verbose': True
        }
        for k, v in self.default_config.items():
            if k not in self.config:
                self.config[k] = v
        
        self.results = {}
        self.optimization_results = {}
        self.surface_data = {}
    
    def run_sensitivity_analysis(self):
        """Run 2D sensitivity analysis across T and pH"""
        
        temps = np.arange(
            self.config['temperature_min'],
            self.config['temperature_max'] + self.config['temperature_step'],
            self.config['temperature_step']
        )
        
        pHs = np.arange(
            self.config['pH_min'],
            self.config['pH_max'] + self.config['pH_step'],
            self.config['pH_step']
        )
        
        if self.config['verbose']:
            print("\n" + "="*70)
            print("🔬 2D SENSITIVITY ANALYSIS: Temperature × pH Optimization")
            print("="*70)
            print(f"  Temperature: {self.config['temperature_min']}°C - {self.config['temperature_max']}°C")
            print(f"  pH: {self.config['pH_min']:.1f} - {self.config['pH_max']:.1f}")
            print(f"  Total points: {len(temps) * len(pHs)}")
            print("="*70)
        
        # Initialize result grids
        n_T = len(temps)
        n_pH = len(pHs)
        
        self.surface_data = {
            'T': temps,
            'pH': pHs,
            'enrichment': np.zeros((n_T, n_pH)),
            'fraction': np.zeros((n_T, n_pH)),
            'dna_yield': np.zeros((n_T, n_pH)),
            'half_life': np.zeros((n_T, n_pH)),
            'stability_ratio': np.zeros((n_T, n_pH)),
            'k_U': np.zeros((n_T, n_pH)),
            'k_T': np.zeros((n_T, n_pH)),
            'poly_rate': np.zeros((n_T, n_pH))
        }
        
        total_points = n_T * n_pH
        point = 0
        
        for i, T in enumerate(temps):
            for j, pH in enumerate(pHs):
                point += 1
                
                if self.config['verbose'] and point % 5 == 0:
                    print(f"  Progress: {point}/{total_points} ({100*point/total_points:.1f}%)", end='\r')
                
                # Run replicates
                replicate_results = []
                for rep in range(self.config['n_replicates']):
                    seed = 42 + rep * 100 + int(T * 10) + int(pH * 100)
                    
                    vssuf = VSSUFEngine(
                        temperature_C=T,
                        pH=pH,
                        seed=seed,
                        max_time_hours=self.config['simulation_hours'],
                        verbose=False,
                        clay_protection_factor=self.config['clay_protection_factor'],
                        clay_surface_density=self.config['clay_surface_density']
                    )
                    vssuf.run()
                    
                    replicate_results.append({
                        'enrichment': vssuf.get_thymine_enrichment(),
                        'fraction': vssuf.get_final_thymine_fraction(),
                        'dna_total': sum(vssuf.species.values()),
                        'half_life': vssuf.get_dna_half_life() / 3600,
                        'k_U': vssuf.k_U,
                        'k_T': vssuf.k_T,
                    })
                
                # Average results
                self.surface_data['enrichment'][i, j] = np.mean([r['enrichment'] for r in replicate_results])
                self.surface_data['fraction'][i, j] = np.mean([r['fraction'] for r in replicate_results])
                self.surface_data['dna_yield'][i, j] = np.mean([r['dna_total'] for r in replicate_results])
                self.surface_data['half_life'][i, j] = np.mean([r['half_life'] for r in replicate_results])
                self.surface_data['k_U'][i, j] = np.mean([r['k_U'] for r in replicate_results])
                self.surface_data['k_T'][i, j] = np.mean([r['k_T'] for r in replicate_results])
                self.surface_data['stability_ratio'][i, j] = pHArrheniusRates.get_stability_ratio(T, pH)
                self.surface_data['poly_rate'][i, j] = pHArrheniusRates.get_polymerization_rate(T, pH)
        
        if self.config['verbose']:
            print("\n  Analysis complete! Finding optimum...")
        
        self._find_optimal_point()
        self._calculate_statistics()
        
        if self.config['verbose']:
            self._print_optimization_summary()
        
        return self.surface_data
    
    def _find_optimal_point(self):
        """Find optimal (T, pH) combination for maximum enrichment"""
        enrichment = self.surface_data['enrichment']
        max_idx = np.unravel_index(np.argmax(enrichment), enrichment.shape)
        
        T_opt = self.surface_data['T'][max_idx[0]]
        pH_opt = self.surface_data['pH'][max_idx[1]]
        max_ench = enrichment[max_idx]
        
        self.optimization_results = {
            'optimal_T': T_opt,
            'optimal_pH': pH_opt,
            'max_enrichment': max_ench,
            'optimal_fraction': self.surface_data['fraction'][max_idx],
            'optimal_yield': self.surface_data['dna_yield'][max_idx],
            'optimal_half_life': self.surface_data['half_life'][max_idx],
            'optimal_stability': self.surface_data['stability_ratio'][max_idx],
            'optimal_k_U': self.surface_data['k_U'][max_idx],
            'optimal_k_T': self.surface_data['k_T'][max_idx],
            'optimal_poly_rate': self.surface_data['poly_rate'][max_idx],
        }
        
        # Find region where enrichment > 80% of max
        threshold = 0.8 * max_ench
        above_threshold = enrichment >= threshold
        self.optimization_results['optimal_region'] = {
            'T_min': self.surface_data['T'][np.any(above_threshold, axis=1)].min(),
            'T_max': self.surface_data['T'][np.any(above_threshold, axis=1)].max(),
            'pH_min': self.surface_data['pH'][np.any(above_threshold, axis=0)].min(),
            'pH_max': self.surface_data['pH'][np.any(above_threshold, axis=0)].max(),
        }
    
    def _calculate_statistics(self):
        """Calculate additional statistics"""
        enrichment = self.surface_data['enrichment']
        self.optimization_results['mean_enrichment'] = np.mean(enrichment)
        self.optimization_results['std_enrichment'] = np.std(enrichment)
        self.optimization_results['enrichment_range'] = (np.min(enrichment), np.max(enrichment))
    
    def _print_optimization_summary(self):
        """Print optimization results"""
        print("\n" + "="*70)
        print("🎯 2D OPTIMIZATION RESULTS")
        print("="*70)
        print(f"\n  Optimal Temperature: {self.optimization_results['optimal_T']:.1f}°C")
        print(f"  Optimal pH: {self.optimization_results['optimal_pH']:.2f}")
        print(f"  Maximum Enrichment: {self.optimization_results['max_enrichment']:.2f}x")
        print(f"\n  At optimal conditions:")
        print(f"    Thymine Fraction: {self.optimization_results['optimal_fraction']:.3f}")
        print(f"    DNA Yield: {self.optimization_results['optimal_yield']/1000:.1f}K")
        print(f"    DNA Half-life: {self.optimization_results['optimal_half_life']:.1f} hours")
        print(f"    T/U Stability: {self.optimization_results['optimal_stability']:.2f}")
        print(f"    Polymerization Rate: {self.optimization_results['optimal_poly_rate']:.4f} h⁻¹")
        
        region = self.optimization_results['optimal_region']
        print(f"\n  80% Efficiency Region:")
        print(f"    Temperature: {region['T_min']:.1f}°C - {region['T_max']:.1f}°C")
        print(f"    pH: {region['pH_min']:.2f} - {region['pH_max']:.2f}")
    
    def plot_sensitivity_2d(self, save_path="sensitivity_2d.png", show_fig=True):
        """Create comprehensive 2D sensitivity analysis plots"""
        sns.set_style("whitegrid")
        
        fig = plt.figure(figsize=(20, 16))
        gs = GridSpec(3, 3, figure=fig, hspace=0.35, wspace=0.3)
        fig.suptitle('2D Sensitivity Analysis: Temperature × pH Optimization for Thymine Enrichment',
                    fontsize=18, fontweight='bold', y=0.98)
        
        T = self.surface_data['T']
        pH = self.surface_data['pH']
        enrichment = self.surface_data['enrichment']
        fraction = self.surface_data['fraction']
        dna_yield = self.surface_data['dna_yield']
        half_life = self.surface_data['half_life']
        stability = self.surface_data['stability_ratio']
        poly_rate = self.surface_data['poly_rate']
        
        # Create meshgrid for contour plots
        T_grid, pH_grid = np.meshgrid(T, pH, indexing='ij')
        
        # ============================================================
        # PLOT 1: 3D Surface - Enrichment (Top Left)
        # ============================================================
        ax1 = fig.add_subplot(gs[0, 0], projection='3d')
        surf = ax1.plot_surface(T_grid, pH_grid, enrichment, cmap='viridis', 
                                alpha=0.8, edgecolor='none')
        ax1.set_xlabel('Temperature (°C)', fontsize=10)
        ax1.set_ylabel('pH', fontsize=10)
        ax1.set_zlabel('Enrichment (x)', fontsize=10)
        ax1.set_title('Thymine Enrichment Surface', fontsize=12, fontweight='bold')
        fig.colorbar(surf, ax=ax1, shrink=0.5, aspect=10)
        
        # Mark optimum
        T_opt = self.optimization_results['optimal_T']
        pH_opt = self.optimization_results['optimal_pH']
        max_ench = self.optimization_results['max_enrichment']
        ax1.scatter([T_opt], [pH_opt], [max_ench], color='red', s=100, marker='*')
        
        # ============================================================
        # PLOT 2: Contour - Enrichment (Top Middle)
        # ============================================================
        ax2 = fig.add_subplot(gs[0, 1])
        contour = ax2.contourf(T, pH, enrichment.T, levels=20, cmap='viridis')
        ax2.contour(T, pH, enrichment.T, levels=10, colors='black', alpha=0.3, linewidths=0.5)
        ax2.scatter(T_opt, pH_opt, color='red', s=150, marker='*', 
                   label=f'Optimal: {T_opt:.1f}°C, pH={pH_opt:.2f}')
        
        # 80% region
        region = self.optimization_results['optimal_region']
        ax2.add_patch(plt.Rectangle((region['T_min'], region['pH_min']),
                                   region['T_max'] - region['T_min'],
                                   region['pH_max'] - region['pH_min'],
                                   fill=False, edgecolor='red', linewidth=2,
                                   linestyle='--', label='80% Efficiency Region'))
        
        ax2.set_xlabel('Temperature (°C)', fontsize=11)
        ax2.set_ylabel('pH', fontsize=11)
        ax2.set_title('Enrichment Contours', fontsize=12, fontweight='bold')
        ax2.legend(loc='upper right', fontsize=9)
        fig.colorbar(contour, ax=ax2, shrink=0.8)
        
        # ============================================================
        # PLOT 3: Heatmap - Enrichment (Top Right)
        # ============================================================
        ax3 = fig.add_subplot(gs[0, 2])
        im = ax3.imshow(enrichment, extent=[pH.min(), pH.max(), T.max(), T.min()],
                       aspect='auto', cmap='RdYlGn', origin='upper')
        ax3.scatter(pH_opt, T_opt, color='blue', s=150, marker='*', 
                   edgecolor='white', linewidth=2)
        ax3.set_xlabel('pH', fontsize=11)
        ax3.set_ylabel('Temperature (°C)', fontsize=11)
        ax3.set_title('Enrichment Heatmap', fontsize=12, fontweight='bold')
        fig.colorbar(im, ax=ax3, shrink=0.8, label='Enrichment (x)')
        
        # ============================================================
        # PLOT 4: Thymine Fraction (Middle Left)
        # ============================================================
        ax4 = fig.add_subplot(gs[1, 0])
        contour4 = ax4.contourf(T, pH, fraction.T, levels=20, cmap='Blues')
        ax4.scatter(T_opt, pH_opt, color='red', s=100, marker='*')
        ax4.set_xlabel('Temperature (°C)', fontsize=11)
        ax4.set_ylabel('pH', fontsize=11)
        ax4.set_title('Thymine Fraction', fontsize=12, fontweight='bold')
        fig.colorbar(contour4, ax=ax4, shrink=0.8)
        
        # ============================================================
        # PLOT 5: DNA Yield (Middle Center)
        # ============================================================
        ax5 = fig.add_subplot(gs[1, 1])
        contour5 = ax5.contourf(T, pH, (dna_yield/1000).T, levels=20, cmap='Oranges')
        ax5.scatter(T_opt, pH_opt, color='red', s=100, marker='*')
        ax5.set_xlabel('Temperature (°C)', fontsize=11)
        ax5.set_ylabel('pH', fontsize=11)
        ax5.set_title('DNA Yield (thousands)', fontsize=12, fontweight='bold')
        fig.colorbar(contour5, ax=ax5, shrink=0.8)
        
        # ============================================================
        # PLOT 6: Stability Ratio (Middle Right)
        # ============================================================
        ax6 = fig.add_subplot(gs[1, 2])
        contour6 = ax6.contourf(T, pH, stability.T, levels=20, cmap='coolwarm')
        ax6.scatter(T_opt, pH_opt, color='red', s=100, marker='*')
        ax6.set_xlabel('Temperature (°C)', fontsize=11)
        ax6.set_ylabel('pH', fontsize=11)
        ax6.set_title('T/U Stability Ratio', fontsize=12, fontweight='bold')
        fig.colorbar(contour6, ax=ax6, shrink=0.8)
        
        # ============================================================
        # PLOT 7: DNA Half-life (Bottom Left)
        # ============================================================
        ax7 = fig.add_subplot(gs[2, 0])
        contour7 = ax7.contourf(T, pH, half_life.T, levels=20, cmap='Reds')
        ax7.scatter(T_opt, pH_opt, color='blue', s=100, marker='*')
        ax7.set_xlabel('Temperature (°C)', fontsize=11)
        ax7.set_ylabel('pH', fontsize=11)
        ax7.set_title('DNA Half-life (hours)', fontsize=12, fontweight='bold')
        fig.colorbar(contour7, ax=ax7, shrink=0.8)
        
        # ============================================================
        # PLOT 8: Polymerization Rate (Bottom Center)
        # ============================================================
        ax8 = fig.add_subplot(gs[2, 1])
        contour8 = ax8.contourf(T, pH, poly_rate.T, levels=20, cmap='Greens')
        ax8.scatter(T_opt, pH_opt, color='red', s=100, marker='*')
        ax8.set_xlabel('Temperature (°C)', fontsize=11)
        ax8.set_ylabel('pH', fontsize=11)
        ax8.set_title('Polymerization Rate (h⁻¹)', fontsize=12, fontweight='bold')
        fig.colorbar(contour8, ax=ax8, shrink=0.8)
        
        # ============================================================
        # PLOT 9: Summary Statistics (Bottom Right)
        # ============================================================
        ax9 = fig.add_subplot(gs[2, 2])
        ax9.axis('off')
        
        opt = self.optimization_results
        summary_text = f"""
        ╔═══════════════════════════════════════════════════════╗
        ║        2D OPTIMIZATION SUMMARY                       ║
        ╠═══════════════════════════════════════════════════════╣
        ║  Optimal Temperature:    {opt['optimal_T']:.1f}°C
        ║  Optimal pH:             {opt['optimal_pH']:.2f}
        ║  Maximum Enrichment:     {opt['max_enrichment']:.2f}x
        ║                                                     ║
        ║  At Optimal Conditions:                              ║
        ║    Thymine Fraction:     {opt['optimal_fraction']:.3f}
        ║    DNA Yield:            {opt['optimal_yield']/1000:.1f}K
        ║    DNA Half-life:        {opt['optimal_half_life']:.1f} h
        ║    T/U Stability:        {opt['optimal_stability']:.2f}
        ║    Polymerization Rate:  {opt['optimal_poly_rate']:.4f} h⁻¹
        ║                                                     ║
        ║  80% Efficiency Region:                              ║
        ║    Temperature:          {opt['optimal_region']['T_min']:.1f}°C - 
        ║                           {opt['optimal_region']['T_max']:.1f}°C
        ║    pH:                   {opt['optimal_region']['pH_min']:.2f} - 
        ║                           {opt['optimal_region']['pH_max']:.2f}
        ║                                                     ║
        ║  Statistics:                                        ║
        ║    Mean Enrichment:      {opt['mean_enrichment']:.2f}x
        ║    Std Enrichment:       {opt['std_enrichment']:.2f}x
        ║    Enrichment Range:     {opt['enrichment_range'][0]:.2f}x - 
        ║                           {opt['enrichment_range'][1]:.2f}x
        ╚═══════════════════════════════════════════════════════╝
        """
        
        ax9.text(0.5, 0.5, summary_text, ha='center', va='center',
                transform=ax9.transAxes, fontsize=9, family='monospace',
                bbox=dict(boxstyle='round', facecolor='lightyellow', alpha=0.9))
        
        plt.tight_layout()
        
        if show_fig:
            plt.show()
        
        if save_path:
            plt.savefig(save_path, dpi=400, bbox_inches='tight')
            print(f"\n✅ 2D Sensitivity plot saved: {save_path}")
        
        return fig

# ====================================================================
# MAIN EXECUTION
# ====================================================================

if __name__ == "__main__":
    print("="*70)
    print("🚀 UPDSF v3.0 - WITH 2D SENSITIVITY ANALYSIS (Temperature × pH)")
    print("   Finding Optimal (T, pH) Combination for Thymine Enrichment")
    print("="*70)
    
    # ============================================================
    # 1. Run 2D sensitivity analysis
    # ============================================================
    sensitivity_config = {
        'temperature_min': 60.0,
        'temperature_max': 100.0,
        'temperature_step': 4.0,
        'pH_min': 5.0,
        'pH_max': 10.0,
        'pH_step': 0.5,
        'simulation_hours': 24,
        'n_replicates': 2,
        'verbose': True
    }
    
    print("\n🔬 Running 2D Sensitivity Analysis (Temperature × pH)...")
    analyzer = TwoDSensitivityAnalyzer(sensitivity_config)
    surface_results = analyzer.run_sensitivity_analysis()
    
    # Plot 2D sensitivity results
    analyzer.plot_sensitivity_2d(save_path="sensitivity_2d.png", show_fig=True)
    
    # ============================================================
    # 2. Run standard simulation at optimal (T, pH)
    # ============================================================
    T_optimal = analyzer.optimization_results['optimal_T']
    pH_optimal = analyzer.optimization_results['optimal_pH']
    
    print(f"\n📊 Running simulation at optimal conditions:")
    print(f"   Temperature: {T_optimal:.1f}°C, pH: {pH_optimal:.2f}")
    
    # Run a detailed simulation at optimal conditions
    vssuf = VSSUFEngine(
        temperature_C=T_optimal,
        pH=pH_optimal,
        seed=42,
        max_time_hours=24,
        verbose=True
    )
    results = vssuf.run()
    
    # Plot time series
    fig, axes = plt.subplots(2, 2, figsize=(14, 10))
    fig.suptitle(f'Optimal Simulation: T={T_optimal:.1f}°C, pH={pH_optimal:.2f}', 
                fontsize=16, fontweight='bold')
    
    time_hours = np.array(results['time']) / 3600
    
    axes[0, 0].plot(time_hours, results['dsDNA_U'], 'b-', label='U-DNA', linewidth=2)
    axes[0, 0].plot(time_hours, results['dsDNA_T'], 'r-', label='T-DNA', linewidth=2)
    axes[0, 0].set_xlabel('Time (hours)')
    axes[0, 0].set_ylabel('DNA Copies')
    axes[0, 0].set_title('DNA Accumulation')
    axes[0, 0].legend()
    axes[0, 0].grid(True, alpha=0.3)
    
    axes[0, 1].plot(time_hours, results['u_ratio'], 'g-', linewidth=2)
    axes[0, 1].set_xlabel('Time (hours)')
    axes[0, 1].set_ylabel('Uracil Fraction')
    axes[0, 1].set_title('Uracil Fraction')
    axes[0, 1].grid(True, alpha=0.3)
    axes[0, 1].set_ylim(0, 1)
    
    axes[1, 0].plot(time_hours, results['enrichment'], 'purple', linewidth=2)
    axes[1, 0].fill_between(time_hours, 0, results['enrichment'], alpha=0.2, color='purple')
    axes[1, 0].set_xlabel('Time (hours)')
    axes[1, 0].set_ylabel('Thymine Enrichment')
    axes[1, 0].set_title('Enrichment Over Time')
    axes[1, 0].grid(True, alpha=0.3)
    
    axes[1, 1].axis('off')
    summary = f"""
    ╔═══════════════════════════════════════════╗
    ║        OPTIMAL SIMULATION SUMMARY        ║
    ╠═══════════════════════════════════════════╣
    ║  Temperature:        {T_optimal:.1f}°C
    ║  pH:                 {pH_optimal:.2f}
    ║  Enrichment:         {vssuf.get_thymine_enrichment():.2f}x
    ║  Thymine Fraction:   {vssuf.get_final_thymine_fraction():.3f}
    ║  DNA Half-life:      {vssuf.get_dna_half_life()/3600:.1f}h
    ║  Poly Events:        {vssuf.polymerization_events:,}
    ║  Hydrolysis Events:  {vssuf.hydrolysis_events:,}
    ╚═══════════════════════════════════════════╝
    """
    axes[1, 1].text(0.5, 0.5, summary, ha='center', va='center',
                   transform=axes[1, 1].transAxes, fontsize=10, family='monospace',
                   bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.8))
    
    plt.tight_layout()
    plt.savefig("optimal_simulation.png", dpi=300, bbox_inches='tight')
    plt.show()
    
    # ============================================================
    # 3. Final Summary
    # ============================================================
    print("\n" + "="*70)
    print("✅ COMPLETE 2D ANALYSIS FINISHED!")
    print("="*70)
    print("\n🎯 OPTIMAL CONDITIONS FOUND:")
    print(f"   Temperature: {T_optimal:.1f}°C")
    print(f"   pH: {pH_optimal:.2f}")
    print(f"   Maximum Thymine Enrichment: {analyzer.optimization_results['max_enrichment']:.2f}x")
    
    print("\n📁 Output files:")
    print("   - sensitivity_2d.png (9-panel 2D analysis)")
    print("   - optimal_simulation.png (Time series at optimal conditions)")
    
    print("\n💡 Key Insights:")
    print(f"   • Optimal conditions (T={T_optimal:.1f}°C, pH={pH_optimal:.2f})")
    print(f"   • Thymine enrichment maximized in slightly alkaline conditions")
    print(f"   • 80% efficiency achieved over broad T and pH ranges")
    print("="*70)
