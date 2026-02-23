import os
from datetime import datetime
import numpy as np

def analyze_circular_dichroism_spectra(
    sample_name,
    sample_type,
    wavelength_data,
    cd_signal_data,
    temperature_data=None,
    thermal_cd_data=None,
    output_dir="./",
):
    """
    Analyzes CD data and returns a dictionary of results.
    """
    # 1. Initialize Results Dictionary
    results = {
        "metadata": {
            "sample_name": sample_name,
            "sample_type": sample_type.lower(),
            "analysis_timestamp": datetime.now().strftime('%Y-%m-%d %H:%M:%S'),
            "output_directory": output_dir
        },
        "spectral_analysis": {
            "wavelength": np.array(wavelength_data),
            "cd_signal": np.array(cd_signal_data),
            "classification": None,
            "features": []
        },
        "thermal_stability": None # Default if no thermal data provided
    }

    # 2. Secondary Structure Logic
    w_data = results["spectral_analysis"]["wavelength"]
    s_data = results["spectral_analysis"]["cd_signal"]

    if results["metadata"]["sample_type"] == "protein":
        # Feature detection
        h_sig = np.sum((w_data >= 190) & (w_data <= 195) & (s_data > 0))
        b_sig = np.sum((w_data >= 215) & (w_data <= 220) & (s_data < 0))
        r_sig = np.sum((w_data >= 195) & (w_data <= 200) & (s_data < 0))

        if h_sig > b_sig and h_sig > r_sig:
            struct = "predominantly alpha-helical"
        elif b_sig > h_sig and b_sig > r_sig:
            struct = "predominantly beta-sheet"
        else:
            struct = "mixed or predominantly random coil"
        
        results["spectral_analysis"]["classification"] = struct
        results["spectral_analysis"]["features"] = ["190-195nm alpha-helix", "215-220nm beta-sheet"]

    elif results["metadata"]["sample_type"] == "nucleic_acid":
        g_sig = np.sum((w_data >= 290) & (w_data <= 300) & (s_data > 0))
        b_sig = np.sum((w_data >= 270) & (w_data <= 280) & (s_data > 0))

        if g_sig > 0:
            struct = "G-quadruplex characteristics"
        elif b_sig > 0:
            struct = "B-form characteristics"
        else:
            struct = "non-standard structure"
            
        results["spectral_analysis"]["classification"] = struct
        results["spectral_analysis"]["features"] = ["290-300nm G-quadruplex peak", "270-280nm B-form peak"]

    # 3. Thermal Stability Logic
    if temperature_data is not None and thermal_cd_data is not None:
        t_data = np.array(temperature_data)
        th_sig = np.array(thermal_cd_data)
        
        # Normalize for unfolded fraction
        min_s, max_s = np.min(th_sig), np.max(th_sig)
        unfolded = (th_sig - min_s) / (max_s - min_s)

        # Calculate Tm (Melting Temp)
        tm_idx = np.argmin(np.abs(unfolded - 0.5))
        tm_val = float(t_data[tm_idx])

        # Cooperativity
        t_range = t_data[-1] - t_data[0]
        width = (t_range / len(t_data) * np.sum((unfolded > 0.2) & (unfolded < 0.8)))

        if width < 0.2 * t_range:
            coop = "highly cooperative"
        elif width < 0.4 * t_range:
            coop = "moderately cooperative"
        else:
            coop = "non-cooperative"

        results["thermal_stability"] = {
            "tm_celsius": tm_val,
            "cooperativity": coop,
            "transition_width": float(width),
            "temperature_array": t_data,
            "unfolded_fraction": unfolded
        }

    return results


def analyze_rna_secondary_structure_features(dot_bracket_structure, sequence=None):
    """
    Calculates structural features of an RNA secondary structure and returns 
    the results as a dictionary.
    """
    # 1. Validation Logic
    if not all(c in "().[]{}" for c in dot_bracket_structure):
        return {"error": "Invalid dot-bracket notation"}

    if sequence and len(sequence) != len(dot_bracket_structure):
        return {"error": "Sequence and structure lengths do not match"}

    # 2. Base Pair Extraction (Stack Algorithm)
    pairs = []
    stack = []
    for i, char in enumerate(dot_bracket_structure):
        if char in "([{":
            stack.append((i, char))
        elif char in ")]}":
            if not stack:
                return {"error": "Unbalanced structure (extra closing bracket)"}
            j, opening_char = stack.pop()
            # Match check
            if ((opening_char == "(" and char != ")") or 
                (opening_char == "[" and char != "]") or 
                (opening_char == "{" and char != "}")):
                return {"error": "Mismatched bracket types"}
            pairs.append((j, i))
    
    if stack:
        return {"error": "Unbalanced structure (extra opening bracket)"}

    pairs.sort()

    # 3. Stem and Loop Identification
    stems = []
    current_stem = []
    for i, (start, end) in enumerate(pairs):
        if (i == 0 or start != pairs[i - 1][0] + 1 or end != pairs[i - 1][1] - 1) and current_stem:
            stems.append(current_stem)
            current_stem = []
        current_stem.append((start, end))
    if current_stem: stems.append(current_stem)

    stem_lengths = [len(s) for s in stems]
    
    loops = []
    for i in range(len(stems)):
        last_pair_end = stems[i][-1][1]
        next_stem_start = stems[i + 1][0][0] if i < len(stems) - 1 else len(dot_bracket_structure)
        loop_size = next_stem_start - last_pair_end - 1
        if loop_size > 0:
            loops.append(loop_size)

    # 4. Energy Calculations
    stem_energies = []
    total_energy = 0.0
    if sequence and stems:
        energy_params = {"AU": -0.9, "UA": -0.9, "GC": -2.1, "CG": -2.1, "GU": -0.5, "UG": -0.5}
        for stem in stems:
            s_energy = sum(energy_params.get(sequence[start] + sequence[end], 0) for start, end in stem)
            stem_energies.append(round(s_energy, 2))
        total_energy = round(sum(stem_energies), 2)

    # 5. Assemble Result Dictionary
    total_len = len(dot_bracket_structure)
    paired_count = len(pairs) * 2

    results = {
        "input": {
            "structure": dot_bracket_structure,
            "sequence": sequence,
            "length": total_len
        },
        "statistics": {
            "total_base_pairs": len(pairs),
            "paired_bases_count": paired_count,
            "unpaired_bases_count": total_len - paired_count,
            "percent_paired": round((paired_count / total_len) * 100, 2)
        },
        "stems": {
            "count": len(stems),
            "lengths": stem_lengths,
            "avg_length": round(sum(stem_lengths)/len(stems), 2) if stems else 0,
            "max_length": max(stem_lengths) if stems else 0,
            "details": [
                {
                    "id": i + 1,
                    "length": len(s),
                    "range": (s[0][0], s[-1][1]),
                    "energy_kcal_mol": stem_energies[i] if sequence else None
                } for i, s in enumerate(stems)
            ]
        },
        "loops": {
            "count": len(loops),
            "sizes": loops,
            "avg_size": round(sum(loops)/len(loops), 2) if loops else 0,
            "max_size": max(loops) if loops else 0
        }
    }

    if sequence:
        results["thermodynamics"] = {
            "total_free_energy": total_energy,
            "zipper_stem_energy": stem_energies[0] if stem_lengths and stem_lengths[0] >= 3 else None
        }

    return results

import numpy as np
from scipy.optimize import curve_fit

def analyze_protease_kinetics(
    time_points,
    fluorescence_data,
    substrate_concentrations,
    enzyme_concentration,
):
    """
    Analyzes protease kinetics and returns a structured dictionary of 
    Michaelis-Menten parameters and initial velocities.
    """
    # 1. Calculate initial velocities (V0)
    # Slope of the first 20% of the reaction time
    initial_velocities = []
    num_points = max(5, int(len(time_points) * 0.2))
    
    for curve in fluorescence_data:
        slope, _ = np.polyfit(time_points[:num_points], curve[:num_points], 1)
        initial_velocities.append(slope)
    
    initial_velocities = np.array(initial_velocities)

    # 2. Michaelis-Menten Model
    def michaelis_menten(s, vmax, km):
        return vmax * s / (km + s)

    # 3. Fit data
    try:
        params, covariance = curve_fit(
            michaelis_menten,
            substrate_concentrations,
            initial_velocities,
            p0=[max(initial_velocities), np.mean(substrate_concentrations)],
            bounds=([0, 0], [np.inf, np.inf]),
        )

        vmax, km = params
        std_dev = np.sqrt(np.diag(covariance))
        vmax_std, km_std = std_dev

        # 4. Derived Kinetic Constants
        kcat = vmax / enzyme_concentration
        kcat_std = vmax_std / enzyme_concentration
        
        cat_efficiency = kcat / km
        # Propagation of error for division
        cat_efficiency_std = cat_efficiency * np.sqrt((kcat_std / kcat) ** 2 + (km_std / km) ** 2)

        # 5. Build Result Dictionary
        return {
            "status": "success",
            "experimental_data": {
                "substrate_concentrations": substrate_concentrations.tolist(),
                "initial_velocities": initial_velocities.tolist()
            },
            "parameters": {
                "vmax": {"value": float(vmax), "std_dev": float(vmax_std), "unit": "a.u./s"},
                "km": {"value": float(km), "std_dev": float(km_std), "unit": "μM"},
                "kcat": {"value": float(kcat), "std_dev": float(kcat_std), "unit": "s^-1"},
                "catalytic_efficiency": {
                    "value": float(cat_efficiency), 
                    "std_dev": float(cat_efficiency_std), 
                    "unit": "μM^-1 s^-1"
                }
            },
            "fit_statistics": {
                "enzyme_concentration": enzyme_concentration,
                "n_points": len(substrate_concentrations)
            }
        }

    except Exception as e:
        return {
            "status": "error",
            "message": str(e)
        }


import numpy as np
from scipy.optimize import curve_fit

def analyze_enzyme_kinetics_assay(
    enzyme_name,
    substrate_concentrations,
    enzyme_concentration,
    modulators=None,
    time_points=None,
):
    """
    Analyzes enzyme kinetics and modulator dose-responses, 
    returning results as a structured dictionary.
    """
    # 1. Setup and Time-Course Simulation
    if time_points is None:
        time_points = np.array([0, 5, 10, 15, 20, 30, 45, 60])
    else:
        time_points = np.array(time_points)

    np.random.seed(42)
    max_act, rate_k = 100, 0.05
    tc_activity = max_act * (1 - np.exp(-rate_k * time_points)) + np.random.normal(0, 3, len(time_points))
    
    # Linear range (approx 30% of max)
    lin_idx = np.where(tc_activity >= 0.3 * max_act)[0][0]
    
    results = {
        "metadata": {
            "enzyme": enzyme_name,
            "enzyme_concentration_nm": enzyme_concentration
        },
        "time_course": {
            "time_points": time_points.tolist(),
            "activity": tc_activity.tolist(),
            "linear_range_minutes": [0, float(time_points[lin_idx])]
        },
        "substrate_kinetics": {},
        "modulators": {}
    }

    # 2. Michaelis-Menten Analysis
    def michaelis_menten(s, vmax, km):
        return vmax * s / (km + s)

    true_vmax, true_km = 120, 25
    s_concentrations = np.array(substrate_concentrations)
    activity_vals = michaelis_menten(s_concentrations, true_vmax, true_km) + np.random.normal(0, 5, len(s_concentrations))

    try:
        p, _ = curve_fit(michaelis_menten, s_concentrations, activity_vals, p0=[100, 20], bounds=([0, 0], [500, 200]))
        results["substrate_kinetics"] = {
            "vmax": float(p[0]),
            "km": float(p[1]),
            "data": {"substrate": s_concentrations.tolist(), "activity": activity_vals.tolist()}
        }
    except Exception:
        results["substrate_kinetics"] = {"error": "MM fit failed"}

    # 3. Modulator Analysis
    if modulators:
        def dose_response(x, ic50, hill):
            return 100 / (1 + (x / ic50) ** hill)

        for name, concs in modulators.items():
            ic50_true = np.random.uniform(1, 50)
            # Simulate sigmoidal activity
            mod_acts = [100 / (1 + (c / ic50_true)**1.0) + np.random.normal(0, 3) if c > 0 else 100 for c in concs]
            
            mod_entry = {"concentrations": concs, "activity_percent": mod_acts, "ic50": None, "hill_coefficient": None}

            if len(concs) >= 4:
                try:
                    nz_c = np.array([c for c in concs if c > 0])
                    nz_a = np.array([a for c, a in zip(concs, mod_acts) if c > 0])
                    p_mod, _ = curve_fit(dose_response, nz_c, nz_a, p0=[10, 1], bounds=([0.1, 0.1], [1000, 10]))
                    mod_entry["ic50"] = float(p_mod[0])
                    mod_entry["hill_coefficient"] = float(p_mod[1])
                except Exception:
                    mod_entry["error"] = "Dose-response fit failed"
            
            results["modulators"][name] = mod_entry

    return results


import numpy as np
import pandas as pd
from scipy.optimize import curve_fit

def analyze_itc_binding_thermodynamics(
    itc_data_path=None,
    itc_data=None,
    temperature=298.15,
    protein_concentration=None,
    ligand_concentration=None,
):
    """
    Analyzes ITC data and returns a structured dictionary of 
    thermodynamic binding parameters.
    """
    # 1. Data Ingestion
    if itc_data_path:
        data_df = pd.read_csv(itc_data_path, sep=None, engine='python')
        data = data_df[["injection", "volume", "heat"]].values
    elif itc_data is not None:
        data = np.array(itc_data)
    else:
        return {"status": "error", "message": "No data provided"}

    # Concentrations (Molar)
    p_conc = protein_concentration if protein_concentration else 1.0
    l_conc = ligand_concentration if ligand_concentration else 10.0
    cell_vol = 1.4  # mL
    
    injections, volumes, heats = data[:, 0], data[:, 1], data[:, 2]
    
    # 2. Pre-processing: Molar Ratios and Dilution
    cum_vol = np.cumsum(volumes)
    dilution = 1 - (cum_vol / cell_vol)
    p_conc_corr = p_conc * dilution

    molar_ratios = np.zeros_like(injections)
    for i in range(len(injections)):
        l_added = volumes[i] * l_conc / 1000
        if i == 0:
            molar_ratios[i] = l_added / (p_conc * cell_vol / 1000)
        else:
            molar_ratios[i] = molar_ratios[i-1] + l_added / (p_conc_corr[i] * (cell_vol - cum_vol[i]) / 1000)

    # 3. Fitting Model
    def one_site_model(x, Kd, dH, n):
        Ka = 1 / Kd
        q = np.zeros_like(x)
        for i in range(len(x)):
            # Simplified binding heat calculation
            term = (n * p_conc_corr[i] * Ka * (x[i] * p_conc_corr[i])) / (1 + Ka * (x[i] * p_conc_corr[i]))
            if i == 0:
                q[i] = term * dH * cell_vol
            else:
                prev_term = (n * p_conc_corr[i-1] * Ka * (x[i-1] * p_conc_corr[i-1])) / (1 + Ka * (x[i-1] * p_conc_corr[i-1]))
                q[i] = term * dH * (cell_vol - cum_vol[i]) - prev_term * dH * (cell_vol - cum_vol[i-1])
        return q

    try:
        popt, pcov = curve_fit(one_site_model, molar_ratios, heats, p0=[1e-6, -5000, 1.0], maxfev=10000)
        Kd, dH, n = popt
        perr = np.sqrt(np.diag(pcov))

        # Thermodynamics
        R = 1.9872
        dG = R * temperature * np.log(Kd)
        dS = (dH - dG) / temperature

        # Goodness of fit
        residuals = heats - one_site_model(molar_ratios, *popt)
        r_sq = 1 - (np.sum(residuals**2) / np.sum((heats - np.mean(heats))**2))

        return {
            "status": "success",
            "parameters": {
                "n_stoichiometry": {"value": float(n), "stderr": float(perr[2])},
                "kd_molar": {"value": float(Kd), "stderr": float(perr[0])},
                "kd_micromolar": {"value": float(Kd * 1e6), "stderr": float(perr[0] * 1e6)},
                "delta_h_cal_mol": {"value": float(dH), "stderr": float(perr[1])},
                "delta_g_cal_mol": {"value": float(dG)},
                "delta_s_cal_mol_k": {"value": float(dS)}
            },
            "fit_quality": {"r_squared": float(r_sq)},
            "experimental_conditions": {"temperature_k": temperature, "cell_volume_ml": cell_vol}
        }

    except Exception as e:
        return {"status": "error", "message": str(e)}


import numpy as np
from Bio import AlignIO, Phylo, SeqIO
from Bio.Align import MultipleSeqAlignment
from Bio.Phylo.TreeConstruction import DistanceCalculator, DistanceTreeConstructor
from io import StringIO

def analyze_protein_conservation(protein_sequences):
    """
    Analyzes protein conservation and returns a dictionary of alignment stats,
    phylogenetic tree data, and conserved positions.
    """
    # 1. Prepare Sequences
    seq_records = []
    if isinstance(protein_sequences, list):
        for i, seq in enumerate(protein_sequences):
            if ">" in seq: # Handle FASTA strings
                record = SeqIO.read(StringIO(seq), "fasta")
            else: # Handle raw strings
                from Bio.Seq import Seq
                from Bio.SeqRecord import SeqRecord
                record = SeqRecord(Seq(seq), id=f"Seq_{i+1}")
            seq_records.append(record)
    
    # 2. Simple Alignment (Padding approach for standalone portability)
    max_len = max(len(s.seq) for s in seq_records)
    for record in seq_records:
        record.seq = record.seq.ljust(max_len, "-")
    
    msa = MultipleSeqAlignment(seq_records)
    
    # 3. Phylogenetic Analysis
    calculator = DistanceCalculator("identity")
    dm = calculator.get_distance(msa)
    constructor = DistanceTreeConstructor()
    tree = constructor.nj(dm)
    
    # Convert tree to Newick string for the dictionary
    tree_io = StringIO()
    Phylo.write(tree, tree_io, "newick")
    newick_string = tree_io.getvalue().strip()

    # 4. Conservation Analysis
    alignment_length = msa.get_alignment_length()
    conservation_data = []
    highly_conserved_indices = []

    for i in range(alignment_length):
        column = msa[:, i]
        most_common = max(column, key=column.count)
        score = column.count(most_common) / len(column)
        
        conservation_data.append({
            "position": i + 1,
            "consensus": most_common,
            "score": round(score, 3)
        })
        
        if score > 0.8:
            highly_conserved_indices.append(i + 1)

    # 5. Assemble Result Dictionary
    return {
        "alignment_metadata": {
            "num_sequences": len(seq_records),
            "alignment_length": alignment_length,
        },
        "phylogeny": {
            "newick_tree": newick_string,
            "method": "Neighbor-Joining (Identity distance)"
        },
        "conservation": {
            "highly_conserved_count": len(highly_conserved_indices),
            "highly_conserved_positions": highly_conserved_indices,
            "full_scores": conservation_data
        }
    }
