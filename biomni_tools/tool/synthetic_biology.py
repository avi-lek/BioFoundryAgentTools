def engineer_bacterial_genome_for_therapeutic_delivery(bacterial_genome_file, genetic_parts):
    """Engineer a bacterial genome by integrating therapeutic genetic parts for therapeutic delivery.

    Parameters
    ----------
    bacterial_genome_file : str
        Path to the file containing the bacterial genome sequence in FASTA format
    genetic_parts : dict
        Dictionary containing genetic parts to be integrated::

        {
            'promoters': list of dict with 'name', 'sequence', and 'position',
            'genes': list of dict with 'name', 'sequence', and 'position',
            'terminators': list of dict with 'name', 'sequence', and 'position',
            'cargo': dict with 'name' and 'sequence' of the therapeutic cargo
        }

    Returns
    -------
    dict
        Dictionary containing engineered genome sequence, features, and summary details.
        No files are written.

    """
    import datetime
    import io

    from Bio import SeqIO
    from Bio.Seq import Seq
    from Bio.SeqFeature import FeatureLocation, SeqFeature
    from Bio.SeqRecord import SeqRecord

    report = {
        "date": datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S"),
        "inputs": {
            "bacterial_genome_file": bacterial_genome_file,
            "genetic_parts": genetic_parts,
        },
        "genome": {},
        "features": [],
        "summary": {},
        "suggested_files": {},
    }

    # Step 1: Load the bacterial genome
    try:
        genome_record = SeqIO.read(bacterial_genome_file, "fasta")
        genome_seq = genome_record.seq
    except Exception as e:
        report["error"] = f"Error loading genome: {str(e)}"
        return report

    # Step 2: Design and integrate genetic parts
    # Create a new genome sequence for modifications
    # In newer BioPython versions, Seq objects are mutable by default
    # or we can convert to a string and modify it, then convert back to Seq
    engineered_seq = str(genome_seq)

    # Track features for visualization
    features = []

    # Keep track of position adjustments as we add sequences
    position_adjustment = 0

    # Add promoters
    if "promoters" in genetic_parts:
        for promoter in genetic_parts["promoters"]:
            position = promoter["position"] + position_adjustment

            # Insert the sequence at the position
            engineered_seq = engineered_seq[:position] + promoter["sequence"] + engineered_seq[position:]
            position_adjustment += len(promoter["sequence"])

            # Adjust positions of subsequent elements
            for part_type in ["genes", "terminators"]:
                if part_type in genetic_parts:
                    for part in genetic_parts[part_type]:
                        if part["position"] > promoter["position"]:
                            part["position"] += len(promoter["sequence"])

            # Add feature for visualization
            feature = SeqFeature(
                FeatureLocation(position, position + len(promoter["sequence"])),
                type="promoter",
                qualifiers={"label": promoter["name"]},
            )
            features.append(feature)
            report["features"].append(
                {
                    "type": "promoter",
                    "name": promoter["name"],
                    "start": int(position),
                    "end": int(position + len(promoter["sequence"])),
                }
            )

    # Add genes
    if "genes" in genetic_parts:
        for gene in genetic_parts["genes"]:
            position = gene["position"] + position_adjustment

            # Insert the sequence at the position
            engineered_seq = engineered_seq[:position] + gene["sequence"] + engineered_seq[position:]
            position_adjustment += len(gene["sequence"])

            # Adjust positions of subsequent elements
            for part_type in ["terminators"]:
                if part_type in genetic_parts:
                    for part in genetic_parts[part_type]:
                        if part["position"] > gene["position"]:
                            part["position"] += len(gene["sequence"])

            # Add feature for visualization
            feature = SeqFeature(
                FeatureLocation(position, position + len(gene["sequence"])),
                type="gene",
                qualifiers={"label": gene["name"]},
            )
            features.append(feature)
            report["features"].append(
                {
                    "type": "gene",
                    "name": gene["name"],
                    "start": int(position),
                    "end": int(position + len(gene["sequence"])),
                }
            )

    # Add terminators
    if "terminators" in genetic_parts:
        for terminator in genetic_parts["terminators"]:
            position = terminator["position"] + position_adjustment

            # Insert the sequence at the position
            engineered_seq = engineered_seq[:position] + terminator["sequence"] + engineered_seq[position:]
            position_adjustment += len(terminator["sequence"])

            # Add feature for visualization
            feature = SeqFeature(
                FeatureLocation(position, position + len(terminator["sequence"])),
                type="terminator",
                qualifiers={"label": terminator["name"]},
            )
            features.append(feature)
            report["features"].append(
                {
                    "type": "terminator",
                    "name": terminator["name"],
                    "start": int(position),
                    "end": int(position + len(terminator["sequence"])),
                }
            )

    # Add therapeutic cargo
    if "cargo" in genetic_parts:
        cargo = genetic_parts["cargo"]
        # Find a suitable position (after the last added element)
        position = max([f.location.end for f in features]) if features else len(engineered_seq) // 2

        # Insert the sequence at the position
        engineered_seq = engineered_seq[:position] + cargo["sequence"] + engineered_seq[position:]

        # Add feature for visualization
        feature = SeqFeature(
            FeatureLocation(position, position + len(cargo["sequence"])),
            type="therapeutic_cargo",
            qualifiers={"label": cargo["name"]},
        )
        features.append(feature)
        report["features"].append(
            {
                "type": "therapeutic_cargo",
                "name": cargo["name"],
                "start": int(position),
                "end": int(position + len(cargo["sequence"])),
            }
        )

    # Step 3: Create the engineered genome record
    engineered_genome = SeqRecord(
        seq=Seq(engineered_seq),
        id=f"{genome_record.id}_engineered",
        name=f"{genome_record.id}_engineered",
        description=f"Engineered {genome_record.id} for therapeutic delivery",
    )

    # Add all features to the engineered genome
    engineered_genome.features = features

    fasta_buf = io.StringIO()
    SeqIO.write(engineered_genome, fasta_buf, "fasta")

    report["genome"] = {
        "id": genome_record.id,
        "length_bp": len(genome_seq),
        "engineered_id": engineered_genome.id,
        "engineered_length_bp": len(engineered_genome.seq),
        "engineered_sequence": engineered_seq,
        "engineered_fasta": fasta_buf.getvalue(),
    }
    report["summary"] = {
        "promoters_added": len(genetic_parts.get("promoters", [])),
        "genes_added": len(genetic_parts.get("genes", [])),
        "terminators_added": len(genetic_parts.get("terminators", [])),
        "cargo_name": genetic_parts.get("cargo", {}).get("name") if "cargo" in genetic_parts else None,
    }
    report["suggested_files"] = {
        "engineered_genome_fasta": f"{genome_record.id}_engineered.fasta",
        "plasmid_map_pdf": f"{genome_record.id}_engineered_map.pdf",
    }

    return report


def analyze_bacterial_growth_rate(time_points, od_measurements, strain_name="Unknown strain", output_dir="./"):
    """Analyze bacterial growth data and extract growth parameters from OD600 measurements.

    Parameters
    ----------
    time_points : list or numpy.ndarray
        Time points at which OD600 measurements were taken (in hours)
    od_measurements : list or numpy.ndarray
        Optical density (OD600) measurements corresponding to each time point
    strain_name : str, optional
        Name of the bacterial strain being analyzed, default is "Unknown strain"
    output_dir : str, optional
        Directory where to save the output files, default is current directory

    Returns
    -------
    dict
        Dictionary containing fitted parameters, model predictions, and summary statistics.
        No files are written; output_dir is accepted for backward compatibility.

    """
    from datetime import datetime

    import numpy as np
    from scipy.optimize import curve_fit

    # Convert inputs to numpy arrays if they aren't already
    time_points = np.array(time_points)
    od_measurements = np.array(od_measurements)

    # Define the Gompertz growth model function
    def gompertz_model(t, lag, mu_max, A):
        """Gompertz growth model.

        Parameters
        ----------
        t: time points
        lag: lag time (hours)
        mu_max: maximum growth rate (per hour)
        A: carrying capacity (maximum OD)

        """
        return A * np.exp(-np.exp(mu_max * np.exp(1) * (lag - t) / A + 1))

    # Initial parameter guesses
    p0 = [
        np.mean(time_points) / 3,  # lag time guess: 1/3 of the mean time
        0.5,  # growth rate guess: 0.5 per hour
        max(od_measurements) * 1.1,  # carrying capacity: slightly above max OD
    ]

    # Fit the model to the data
    try:
        popt, pcov = curve_fit(
            gompertz_model,
            time_points,
            od_measurements,
            p0=p0,
            bounds=([0, 0, 0], [np.inf, np.inf, np.inf]),
        )
        lag_time, mu_max, carrying_capacity = popt

        # Calculate doubling time (ln(2)/μ)
        doubling_time = np.log(2) / mu_max

        # Generate fitted curve for plotting/data return
        t_smooth = np.linspace(min(time_points), max(time_points), 100)
        od_fit = gompertz_model(t_smooth, *popt)
        return {
            "inputs": {
                "strain_name": strain_name,
                "time_points": time_points.tolist(),
                "od_measurements": od_measurements.tolist(),
                "output_dir": output_dir,
            },
            "summary": {
                "analysis_time": datetime.now().strftime("%Y-%m-%d %H:%M:%S"),
                "num_points": len(time_points),
                "duration_hours": float(max(time_points) - min(time_points)),
            },
            "fit": {
                "lag_time_hours": float(lag_time),
                "mu_max_per_hour": float(mu_max),
                "doubling_time_hours": float(doubling_time),
                "carrying_capacity_od": float(carrying_capacity),
                "covariance": pcov.tolist(),
            },
            "model_curve": {
                "time_points": t_smooth.tolist(),
                "od_fit": od_fit.tolist(),
            },
        }

    except RuntimeError as e:
        return {
            "inputs": {
                "strain_name": strain_name,
                "time_points": time_points.tolist(),
                "od_measurements": od_measurements.tolist(),
                "output_dir": output_dir,
            },
            "error": f"Failed to fit the growth model to the provided data: {str(e)}",
            "suggestions": [
                "Check if the data shows a typical growth pattern",
                "Ensure sufficient data points are provided",
                "Try different initial parameter guesses or a different growth model",
            ],
        }


def analyze_barcode_sequencing_data(
    input_file,
    barcode_pattern=None,
    flanking_seq_5prime=None,
    flanking_seq_3prime=None,
    min_count=5,
    output_dir="./results",
):
    """Analyze sequencing data to extract, quantify and determine lineage relationships of barcodes.

    Parameters
    ----------
    input_file : str
        Path to the input sequencing file in FASTQ or FASTA format
    barcode_pattern : str, optional
        Regular expression pattern to identify barcodes. If None, will use flanking sequences
    flanking_seq_5prime : str, optional
        5' flanking sequence of the barcode region
    flanking_seq_3prime : str, optional
        3' flanking sequence of the barcode region
    min_count : int, default=5
        Minimum count threshold for considering a barcode
    output_dir : str, default="./results"
        Directory to save output files

    Returns
    -------
    dict
        Dictionary containing barcode counts, filtered barcodes, and lineage analysis results.
        No files are written; output_dir is accepted for backward compatibility.

    """
    import re
    from collections import Counter

    import numpy as np
    from Bio import SeqIO
    from scipy.cluster.hierarchy import fcluster, linkage
    from scipy.spatial.distance import pdist

    result = {
        "inputs": {
            "input_file": input_file,
            "barcode_pattern": barcode_pattern,
            "flanking_seq_5prime": flanking_seq_5prime,
            "flanking_seq_3prime": flanking_seq_3prime,
            "min_count": min_count,
            "output_dir": output_dir,
        },
        "summary": {},
        "counts": {},
        "filtered_counts": {},
        "lineages": {},
        "suggested_files": {
            "barcode_counts_tsv": "barcode_counts.tsv",
            "barcode_lineages_tsv": "barcode_lineages.tsv",
        },
    }

    # Step 1: Read sequences and extract barcodes
    # Determine file format based on extension
    file_format = "fastq" if input_file.endswith((".fastq", ".fq")) else "fasta"

    barcodes = []
    total_reads = 0

    for record in SeqIO.parse(input_file, file_format):
        total_reads += 1
        seq_str = str(record.seq)

        # Extract barcode using pattern or flanking sequences
        if barcode_pattern:
            match = re.search(barcode_pattern, seq_str)
            if match:
                barcodes.append(match.group(0))
        elif flanking_seq_5prime and flanking_seq_3prime:
            pattern = f"{flanking_seq_5prime}(.*?){flanking_seq_3prime}"
            match = re.search(pattern, seq_str)
            if match:
                barcodes.append(match.group(1))

    if len(barcodes) == 0:
        result["error"] = "No barcodes found. Check your barcode pattern or flanking sequences."
        result["summary"] = {
            "total_reads": total_reads,
            "barcodes_extracted": 0,
        }
        return result

    # Step 2: Quantify barcode abundances
    barcode_counts = Counter(barcodes)

    # Filter low-abundance barcodes
    filtered_barcodes = {bc: count for bc, count in barcode_counts.items() if count >= min_count}

    result["counts"] = dict(barcode_counts)
    result["filtered_counts"] = dict(filtered_barcodes)

    # Only proceed if we have enough barcodes
    if len(filtered_barcodes) < 2:
        result["lineages"] = {
            "status": "insufficient_barcodes",
            "message": "Not enough barcodes for lineage analysis after filtering",
        }
    else:
        # Convert barcodes to numerical representation for distance calculation
        barcode_list = list(filtered_barcodes.keys())

        # Create a simple distance matrix based on Hamming distance
        def hamming_distance(s1, s2):
            # If sequences have different lengths, pad the shorter one
            if len(s1) != len(s2):
                max_len = max(len(s1), len(s2))
                s1 = s1.ljust(max_len)
                s2 = s2.ljust(max_len)
            return sum(c1 != c2 for c1, c2 in zip(s1, s2, strict=False))

        # Create distance matrix
        n = len(barcode_list)
        dist_matrix = np.zeros((n, n))
        for i in range(n):
            for j in range(i + 1, n):
                dist = hamming_distance(barcode_list[i], barcode_list[j])
                dist_matrix[i, j] = dist
                dist_matrix[j, i] = dist

        # Perform hierarchical clustering
        try:
            condensed_dist = pdist(dist_matrix)
            Z = linkage(condensed_dist, method="average")
            max_dist = 3  # Max Hamming distance to consider barcodes in same lineage
            clusters = fcluster(Z, max_dist, criterion="distance")

            # Count clusters
            cluster_counts = Counter(clusters)

            result["lineages"] = {
                "status": "ok",
                "max_hamming_distance": max_dist,
                "lineage_count": len(cluster_counts),
                "largest_lineage_size": max(cluster_counts.values()),
                "assignments": [
                    {"barcode": bc, "count": filtered_barcodes[bc], "lineage": int(clusters[i])}
                    for i, bc in enumerate(barcode_list)
                ],
            }
        except Exception as e:
            result["lineages"] = {
                "status": "error",
                "message": f"Error in lineage analysis: {str(e)}",
            }

    result["summary"] = {
        "total_reads": total_reads,
        "barcodes_extracted": len(barcodes),
        "unique_barcodes": len(barcode_counts),
        "barcodes_passing_threshold": len(filtered_barcodes),
        "min_count": min_count,
    }

    return result


def analyze_bifurcation_diagram(time_series_data, parameter_values, system_name="Dynamical System", output_dir="./"):
    """Performs bifurcation analysis on a dynamical system and generates a bifurcation diagram.

    Parameters
    ----------
    time_series_data : numpy.ndarray
        A 2D array where each row represents a time series for a specific parameter value.
        Shape should be (n_parameter_values, n_time_points).
    parameter_values : numpy.ndarray
        1D array of parameter values corresponding to each time series.
        Shape should be (n_parameter_values,).
    system_name : str, optional
        Name of the dynamical system being analyzed, used for plot titles.
    output_dir : str, optional
        Directory to save the output files.

    Returns
    -------
    dict
        Dictionary containing regime classification, indicators, and attractor points.
        No files are written; output_dir is accepted for backward compatibility.

    """
    import numpy as np
    from scipy.signal import find_peaks

    # Step 1: Verify input data
    if len(parameter_values) != time_series_data.shape[0]:
        return {
            "error": "Number of parameter values does not match number of time series.",
            "inputs": {
                "system_name": system_name,
                "parameter_values": list(parameter_values),
                "output_dir": output_dir,
            },
        }

    # Initialize arrays to store results
    local_maxima_counts = np.zeros(len(parameter_values))
    lyapunov_estimates = np.zeros(len(parameter_values))
    attractor_points = []

    for i, series in enumerate(time_series_data):
        # Use the last 70% of the time series to avoid transients
        steady_state = series[int(0.7 * len(series)) :]

        # Find local maxima (peaks) to identify periodicity
        peaks, _ = find_peaks(steady_state)
        local_maxima_counts[i] = len(peaks)

        # Simple estimate of the largest Lyapunov exponent
        # Positive values suggest chaos, negative values suggest stability
        if len(steady_state) > 10:
            diffs = np.abs(np.diff(steady_state))
            if np.mean(diffs) > 0:
                lyapunov_estimates[i] = np.log(np.mean(diffs))
            else:
                lyapunov_estimates[i] = -1.0  # Stable behavior

        # Collect points for the bifurcation diagram
        # We use local maxima as sampling points for the attractor
        if len(peaks) > 0:
            attractor_points.append((parameter_values[i], steady_state[peaks]))
        else:
            # If no peaks, use the last few points
            attractor_points.append((parameter_values[i], steady_state[-5:]))

    # Step 3: Classify dynamical regimes

    # Initialize arrays for regime classification
    regimes = np.full(len(parameter_values), "unknown", dtype=object)

    # Classify based on Lyapunov exponent and periodicity
    for i in range(len(parameter_values)):
        if lyapunov_estimates[i] > 0.05:
            regimes[i] = "chaotic"
        elif local_maxima_counts[i] == 0:
            regimes[i] = "stable"
        elif local_maxima_counts[i] == 1:
            regimes[i] = "period-1"
        elif local_maxima_counts[i] == 2:
            regimes[i] = "period-2"
        elif local_maxima_counts[i] == 4:
            regimes[i] = "period-4"
        elif local_maxima_counts[i] > 4:
            if lyapunov_estimates[i] > 0:
                regimes[i] = "chaotic"
            else:
                regimes[i] = f"period-{int(local_maxima_counts[i])}"

    # Identify regime changes
    regime_changes = np.where(regimes[:-1] != regimes[1:])[0]
    regime_transitions = [
        {
            "from": regimes[i],
            "to": regimes[i + 1],
            "parameter_value": float(parameter_values[i + 1]),
        }
        for i in regime_changes
    ]

    unique_regimes = np.unique(regimes)
    regime_counts = {regime: int(np.sum(regimes == regime)) for regime in unique_regimes}

    return {
        "inputs": {
            "system_name": system_name,
            "parameter_values": list(parameter_values),
            "output_dir": output_dir,
        },
        "summary": {
            "num_parameter_values": len(parameter_values),
            "time_points_per_series": int(time_series_data.shape[1]),
            "regime_counts": regime_counts,
            "parameter_range": (float(parameter_values[0]), float(parameter_values[-1])),
        },
        "indicators": {
            "local_maxima_counts": local_maxima_counts.tolist(),
            "lyapunov_estimates": lyapunov_estimates.tolist(),
        },
        "attractor_points": [
            {
                "parameter_value": float(param_val),
                "points": points.flatten().tolist(),
            }
            for param_val, points in attractor_points
        ],
        "regimes": regimes.tolist(),
        "regime_transitions": regime_transitions,
        "suggested_files": {
            "bifurcation_diagram_png": f"{system_name.replace(' ', '_')}_bifurcation_diagram.png",
        },
    }


def create_biochemical_network_sbml_model(reaction_network, kinetic_parameters, output_file="biochemical_model.xml"):
    """Generate a mathematical model of a biochemical network in SBML format.

    Parameters
    ----------
    reaction_network : list of dict
        List of dictionaries, each representing a reaction with keys:
        - 'id': Reaction identifier
        - 'name': Reaction name
        - 'reactants': Dict of reactant species IDs and their stoichiometry
        - 'products': Dict of product species IDs and their stoichiometry
        - 'reversible': Boolean indicating if reaction is reversible

    kinetic_parameters : dict
        Dictionary mapping reaction IDs to their kinetic law parameters.
        Each entry should contain:
        - 'law_type': Type of kinetic law (e.g., 'mass_action', 'michaelis_menten')
        - 'parameters': Dict of parameter names and values

    output_file : str, optional
        File path to save the SBML model (default: "biochemical_model.xml")

    Returns
    -------
    dict
        Dictionary containing SBML string, model summary, and validation results.
        No files are written; output_file is accepted for backward compatibility.

    """
    import libsbml

    # Step 1: Create an SBML document
    sbml_ns = libsbml.SBMLNamespaces(3, 2)  # SBML Level 3 Version 2
    document = libsbml.SBMLDocument(sbml_ns)
    model = document.createModel()
    model.setId("biochemical_network_model")
    model.setName("Biochemical Network Model")

    # Step 2: Create default compartment
    compartment = model.createCompartment()
    compartment.setId("default")
    compartment.setConstant(True)
    compartment.setSize(1.0)
    compartment.setSpatialDimensions(3)

    # Step 3: Create species
    species_set = set()

    # Collect all unique species from reactions
    for reaction in reaction_network:
        for species_id in list(reaction["reactants"].keys()) + list(reaction["products"].keys()):
            species_set.add(species_id)

    # Create species in the model
    for species_id in species_set:
        species = model.createSpecies()
        species.setId(species_id)
        species.setCompartment("default")
        species.setInitialConcentration(0.0)  # Default initial concentration
        species.setHasOnlySubstanceUnits(False)
        species.setBoundaryCondition(False)
        species.setConstant(False)

    # Step 4: Create reactions
    for reaction_data in reaction_network:
        reaction = model.createReaction()
        reaction.setId(reaction_data["id"])
        reaction.setName(reaction_data["name"])
        reaction.setReversible(reaction_data["reversible"])
        reaction.setFast(False)

        # Add reactants
        for reactant_id, stoichiometry in reaction_data["reactants"].items():
            species_ref = reaction.createReactant()
            species_ref.setSpecies(reactant_id)
            species_ref.setStoichiometry(stoichiometry)
            species_ref.setConstant(True)

        # Add products
        for product_id, stoichiometry in reaction_data["products"].items():
            species_ref = reaction.createProduct()
            species_ref.setSpecies(product_id)
            species_ref.setStoichiometry(stoichiometry)
            species_ref.setConstant(True)

        # Create kinetic law
        if reaction_data["id"] in kinetic_parameters:
            kinetic_data = kinetic_parameters[reaction_data["id"]]
            kinetic_law = reaction.createKineticLaw()

            # Add parameters to kinetic law
            for param_name, param_value in kinetic_data["parameters"].items():
                parameter = kinetic_law.createParameter()
                parameter.setId(param_name)
                parameter.setValue(param_value)

            # Set the appropriate formula based on the law type
            if kinetic_data["law_type"] == "mass_action":
                # Simple mass action kinetics
                reactants_formula = " * ".join([f"{reactant_id}" for reactant_id in reaction_data["reactants"]])
                formula = f"k * {reactants_formula}" if reactants_formula else "k"
                kinetic_law.setMath(libsbml.parseL3Formula(formula))
            elif kinetic_data["law_type"] == "michaelis_menten":
                # Michaelis-Menten kinetics (assuming single substrate)
                substrate = list(reaction_data["reactants"].keys())[0] if reaction_data["reactants"] else "S"
                formula = f"Vmax * {substrate} / (Km + {substrate})"
                kinetic_law.setMath(libsbml.parseL3Formula(formula))
            # Custom formula provided in the parameters
            elif "formula" in kinetic_data:
                kinetic_law.setMath(libsbml.parseL3Formula(kinetic_data["formula"]))


    # Step 5: Validate the model
    document.setConsistencyChecks(libsbml.LIBSBML_CAT_UNITS_CONSISTENCY, False)
    consistent = document.checkConsistency()
    errors = []
    if consistent != 0:
        for i in range(document.getNumErrors()):
            error = document.getError(i)
            errors.append(error.getMessage())

    sbml_string = libsbml.writeSBMLToString(document)

    return {
        "summary": {
            "species_count": len(species_set),
            "reaction_count": len(reaction_network),
            "consistent": consistent == 0,
        },
        "validation_errors": errors,
        "sbml": sbml_string,
        "output_file": output_file,
    }


def optimize_codons_for_heterologous_expression(target_sequence, host_codon_usage):
    """Analyzes and optimizes a DNA/RNA sequence for improved expression in a heterologous host organism.

    Parameters
    ----------
    target_sequence : str
        The DNA or RNA sequence of the target gene to be optimized.
        Should contain complete codons (length divisible by 3).

    host_codon_usage : dict
        Dictionary mapping codons to their usage frequency in the host organism.
        Format: {'AUG': 0.8, 'GCC': 0.6, ...} or {'ATG': 0.8, 'GCC': 0.6, ...}

    Returns
    -------
    dict
        Dictionary containing optimized sequence, codon statistics, and warnings.

    """
    from Bio.Data import CodonTable

    # Check if sequence is DNA or RNA and standardize to DNA
    is_rna = "U" in target_sequence
    working_seq = target_sequence.replace("U", "T") if is_rna else target_sequence
    warnings = []

    # Verify sequence length is divisible by 3
    if len(working_seq) % 3 != 0:
        warnings.append(
            f"Sequence length ({len(working_seq)}) is not divisible by 3. Optimization may be incomplete."
        )

    # Extract codons from the sequence
    original_codons = [working_seq[i : i + 3] for i in range(0, len(working_seq), 3)]

    # Get standard genetic code
    standard_table = CodonTable.standard_dna_table

    # Create codon to amino acid mapping
    codon_to_aa = {}
    for codon, amino_acid in standard_table.forward_table.items():
        codon_to_aa[codon] = amino_acid

    # Add stop codons
    for stop_codon in standard_table.stop_codons:
        codon_to_aa[stop_codon] = "*"  # * represents stop codon

    # Create amino acid to codon mapping with host frequencies
    aa_to_codons = {}
    for codon, aa in codon_to_aa.items():
        if aa not in aa_to_codons:
            aa_to_codons[aa] = []
        # Use host frequency or default to 0 if not provided
        frequency = host_codon_usage.get(codon, 0)
        aa_to_codons[aa].append((codon, frequency))

    # Sort codons by frequency for each amino acid
    for aa in aa_to_codons:
        aa_to_codons[aa] = sorted(aa_to_codons[aa], key=lambda x: x[1], reverse=True)

    # Optimize codons
    optimized_codons = []
    for codon in original_codons:
        try:
            aa = codon_to_aa[codon]
            # Select the highest frequency codon for this amino acid
            best_codon = aa_to_codons[aa][0][0] if aa_to_codons[aa] else codon
            optimized_codons.append(best_codon)
        except KeyError:
            # If codon not in standard table, keep original
            optimized_codons.append(codon)
            warnings.append(f"Non-standard codon detected: {codon}. Keeping original.")

    # Combine optimized codons into sequence
    optimized_sequence = "".join(optimized_codons)

    # Convert back to RNA if input was RNA
    if is_rna:
        optimized_sequence = optimized_sequence.replace("T", "U")

    # Calculate optimization statistics
    codon_changes = sum(1 for i, j in zip(original_codons, optimized_codons, strict=False) if i != j)
    percent_changed = (codon_changes / len(original_codons)) * 100 if original_codons else 0

    return {
        "inputs": {
            "sequence_type": "RNA" if is_rna else "DNA",
            "sequence_length": len(working_seq),
        },
        "original_sequence": target_sequence,
        "optimized_sequence": optimized_sequence,
        "stats": {
            "num_codons": len(original_codons),
            "codons_modified": codon_changes,
            "percent_modified": percent_changed,
        },
        "warnings": warnings,
        "suggested_files": {
            "original_sequence_txt": "original_sequence.txt",
            "optimized_sequence_txt": "optimized_sequence.txt",
        },
    }


def simulate_gene_circuit_with_growth_feedback(
    circuit_topology,
    kinetic_params,
    growth_params,
    simulation_time=100,
    time_points=1000,
):
    """Simulate gene regulatory circuit dynamics with growth feedback.

    Parameters
    ----------
    circuit_topology : numpy.ndarray
        Adjacency matrix representing the gene circuit topology.
        Positive values indicate activation, negative values indicate repression.
        Shape should be (n_genes, n_genes) where n_genes is the number of genes in the circuit.

    kinetic_params : dict
        Dictionary containing kinetic parameters:
        - 'basal_rates': list of basal expression rates for each gene
        - 'degradation_rates': list of degradation rates for each gene
        - 'hill_coefficients': list of Hill coefficients for regulatory interactions
        - 'threshold_constants': list of threshold constants for regulatory interactions

    growth_params : dict
        Dictionary containing growth-related parameters:
        - 'max_growth_rate': maximum cell growth rate
        - 'growth_inhibition': how gene expression affects growth
        - 'gene_growth_weights': weights for how each gene affects growth

    simulation_time : float, optional
        Total simulation time (default: 100)

    time_points : int, optional
        Number of time points to sample (default: 1000)

    Returns
    -------
    dict
        Dictionary containing simulation parameters, time series results, and summary metrics.
        No files are written.

    """
    import datetime

    import numpy as np
    from scipy.integrate import solve_ivp

    # Extract parameters
    n_genes = circuit_topology.shape[0]
    basal_rates = kinetic_params["basal_rates"]
    degradation_rates = kinetic_params["degradation_rates"]
    hill_coefficients = kinetic_params["hill_coefficients"]
    threshold_constants = kinetic_params["threshold_constants"]

    max_growth_rate = growth_params["max_growth_rate"]
    growth_inhibition = growth_params["growth_inhibition"]
    gene_growth_weights = growth_params["gene_growth_weights"]

    # Define the ODE system
    def gene_circuit_odes(t, state):
        # Extract state variables
        gene_expressions = state[:n_genes]
        cell_mass = state[n_genes]

        # Initialize derivatives
        dxdt = np.zeros(n_genes + 1)

        # Gene expression dynamics
        for i in range(n_genes):
            # Basal expression
            production = basal_rates[i]

            # Regulatory influences
            for j in range(n_genes):
                if circuit_topology[i, j] != 0:
                    # Calculate regulatory effect
                    regulation = gene_expressions[j] ** hill_coefficients[j] / (
                        threshold_constants[j] ** hill_coefficients[j] + gene_expressions[j] ** hill_coefficients[j]
                    )

                    if circuit_topology[i, j] > 0:  # Activation
                        production *= 1 + circuit_topology[i, j] * regulation
                    else:  # Repression
                        production *= 1 + circuit_topology[i, j] * (1 - regulation)

            # Dilution due to growth and degradation
            dilution = degradation_rates[i] + (dxdt[n_genes] / cell_mass if cell_mass > 0 else 0)

            # Final rate equation for gene i
            dxdt[i] = production - dilution * gene_expressions[i]

        # Cell growth dynamics
        growth_burden = sum(gene_growth_weights[i] * gene_expressions[i] for i in range(n_genes))
        dxdt[n_genes] = max_growth_rate * cell_mass / (1 + growth_inhibition * growth_burden)

        return dxdt

    # Initial conditions (starting with low gene expression and unit cell mass)
    initial_state = np.zeros(n_genes + 1)
    initial_state[:n_genes] = 0.1  # Low initial gene expression
    initial_state[n_genes] = 1.0  # Initial cell mass

    # Solve the ODE system
    t_span = (0, simulation_time)
    t_eval = np.linspace(0, simulation_time, time_points)

    solution = solve_ivp(
        gene_circuit_odes,
        t_span,
        initial_state,
        method="LSODA",
        t_eval=t_eval,
        rtol=1e-6,
        atol=1e-9,
    )

    # Extract results
    time_series = solution.t
    gene_expression_time_series = solution.y[:n_genes, :]
    cell_growth_time_series = solution.y[n_genes, :]

    timestamp = datetime.datetime.now().strftime("%Y%m%d_%H%M%S")
    final_gene_expression = gene_expression_time_series[:, -1].tolist()
    final_cell_mass = float(cell_growth_time_series[-1])
    final_growth_rate = float(
        max_growth_rate
        / (1 + growth_inhibition * sum(gene_growth_weights[i] * gene_expression_time_series[i, -1] for i in range(n_genes)))
    )

    return {
        "timestamp": timestamp,
        "inputs": {
            "circuit_topology": circuit_topology.tolist(),
            "kinetic_params": kinetic_params,
            "growth_params": growth_params,
            "simulation_time": simulation_time,
            "time_points": time_points,
        },
        "summary": {
            "num_genes": n_genes,
            "final_gene_expression_levels": final_gene_expression,
            "final_cell_mass": final_cell_mass,
            "final_growth_rate": final_growth_rate,
        },
        "results": {
            "time": time_series.tolist(),
            "gene_expression": gene_expression_time_series.tolist(),
            "cell_growth": cell_growth_time_series.tolist(),
        },
        "suggested_files": {
            "time_series_npz": f"gene_circuit_simulation_{timestamp}.npz",
            "params_json": f"simulation_params_{timestamp}.json",
        },
    }


def identify_fas_functional_domains(sequence, sequence_type="protein", output_file="fas_domains_report.txt"):
    """Identifies functional domains within a Fatty Acid Synthase (FAS) sequence and predicts their roles.

    Parameters
    ----------
    sequence : str
        The nucleotide or protein sequence of a FAS gene
    sequence_type : str
        Type of sequence provided - "protein" or "nucleotide" (default: "protein")
    output_file : str
        Name of the output file to save the detailed domain report (default: "fas_domains_report.txt")

    Returns
    -------
    dict
        Dictionary containing identified domains and FAS-specific domain summary.
        No files are written; output_file is accepted for backward compatibility.

    """
    import json
    import time

    import requests
    from Bio.Seq import Seq

    result = {
        "inputs": {
            "sequence_type": sequence_type,
            "sequence_length": len(sequence),
            "output_file": output_file,
        },
        "translation": {},
        "domains_found": [],
        "fas_domains_found": [],
        "summary": {},
        "suggested_files": {"report_txt": output_file},
    }

    # Convert nucleotide to protein if needed
    if sequence_type == "nucleotide":
        try:
            protein_seq = str(Seq(sequence).translate())
            result["translation"] = {
                "status": "ok",
                "protein_length": len(protein_seq),
            }
        except Exception as e:
            result["translation"] = {"status": "error", "message": str(e)}
            result["error"] = f"Error in translation: {str(e)}"
            return result
    else:
        protein_seq = sequence

    # FAS domain information dictionary
    fas_domains = {
        "ketoacyl-synt": {
            "name": "Ketoacyl Synthase (KS)",
            "function": "Catalyzes the condensation of malonyl-ACP with the growing fatty acid chain",
        },
        "Ketoacyl-synt_C": {
            "name": "Ketoacyl Synthase C-terminal (KS-C)",
            "function": "C-terminal domain of KS, involved in substrate binding",
        },
        "Acyl_transf_1": {
            "name": "Acyltransferase (AT)",
            "function": "Transfers acyl groups from acyl-CoA to ACP",
        },
        "ketoacyl-red": {
            "name": "Ketoacyl Reductase (KR)",
            "function": "Reduces the beta-ketoacyl group to a beta-hydroxyacyl group",
        },
        "Thioesterase": {
            "name": "Thioesterase (TE)",
            "function": "Releases the fatty acid from ACP by hydrolysis",
        },
        "PS-DH": {
            "name": "Dehydratase (DH)",
            "function": "Removes water from beta-hydroxyacyl-ACP to form trans-2-enoyl-ACP",
        },
        "Enoyl_red": {
            "name": "Enoyl Reductase (ER)",
            "function": "Reduces the double bond in enoyl-ACP to a saturated acyl-ACP",
        },
        "ACP": {
            "name": "Acyl Carrier Protein (ACP)",
            "function": "Carries the growing fatty acid chain between enzyme domains",
        },
        "PP-binding": {
            "name": "Phosphopantetheine Binding Domain",
            "function": "Attachment site for the phosphopantetheine prosthetic group in ACP",
        },
    }

    # HMMER web API for domain identification
    url = "https://www.ebi.ac.uk/Tools/hmmer/search/hmmscan"
    headers = {"Content-Type": "application/json", "Accept": "application/json"}
    data = {"hmmdb": "pfam", "seq": protein_seq}

    try:
        response = requests.post(url, data=json.dumps(data), headers=headers)
        if response.status_code != 200:
            result["error"] = f"HMMER API returned status code {response.status_code}"
            return result

        result_url = response.headers.get("Location", "")
        if not result_url:
            result["error"] = "No result URL returned from HMMER API"
            return result

        # Wait for results to be ready
        time.sleep(2)

        # Get results
        result_response = requests.get(f"{result_url}.json")
        if result_response.status_code != 200:
            result["error"] = "Could not retrieve results from HMMER API"
            return result

        results = result_response.json()

        # Process results
        domains_found = []
        if "results" in results and "hits" in results["results"]:
            hits = results["results"]["hits"]
            for hit in hits:
                if "domains" in hit:
                    for domain in hit["domains"]:
                        domain_name = hit.get("name", "").split(".")[0]
                        domain_desc = hit.get("desc", "Unknown")
                        domain_start = domain.get("ali_from", 0)
                        domain_end = domain.get("ali_to", 0)

                        domains_found.append(
                            {
                                "name": domain_name,
                                "description": domain_desc,
                                "start": domain_start,
                                "end": domain_end,
                                "score": domain.get("score", 0),
                            }
                        )

        # Check for FAS domains
        fas_domains_found = []
        for domain in domains_found:
            for fas_key, fas_info in fas_domains.items():
                if fas_key.lower() in domain["name"].lower():
                    fas_domains_found.append(
                        {
                            "name": fas_info["name"],
                            "pfam_id": domain["name"],
                            "start": domain["start"],
                            "end": domain["end"],
                            "function": fas_info["function"],
                        }
                    )

        result["domains_found"] = domains_found
        result["fas_domains_found"] = fas_domains_found
        result["summary"] = {
            "total_domains_found": len(domains_found),
            "fas_domains_found": len(fas_domains_found),
            "fas_domains_present": len(fas_domains_found) > 0,
            "interpretation": (
                "The sequence contains multiple domains characteristic of Fatty Acid Synthase"
                if fas_domains_found
                else "The sequence does not contain typical FAS domains; it may not be FAS or may be partial."
            ),
        }

    except Exception as e:
        result["error"] = f"Error in domain identification: {str(e)}"

    return result
