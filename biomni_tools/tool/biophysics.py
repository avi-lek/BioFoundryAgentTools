def predict_protein_disorder_regions(protein_sequence, threshold=0.5, output_file="disorder_prediction_results.csv"):
    """Predicts intrinsically disordered regions (IDRs) in a protein sequence using IUPred2A.

    Parameters
    ----------
    protein_sequence : str
        The amino acid sequence of the protein to analyze
    threshold : float, optional
        The disorder score threshold above which a residue is considered disordered (default: 0.5)
    output_file : str, optional
        Filename to save the per-residue disorder scores (default: "disorder_prediction_results.csv")

    Returns
    -------
    dict
        Dictionary containing per-residue disorder scores, disordered regions, and summary.
        No files are written; output_file is accepted for backward compatibility.

    """
    import re

    import requests

    # Clean the input sequence
    protein_sequence = "".join(re.findall(r"[A-Za-z]", protein_sequence))

    # Step 1: Submit the sequence to IUPred2A web server
    url = "https://iupred2a.elte.hu/iupred2a"
    payload = {
        "seq": protein_sequence,
        "iupred2": "long",  # Use IUPred2 long disorder prediction
        "anchor2": "no",  # Don't use ANCHOR2 prediction
    }

    try:
        response = requests.post(url, data=payload)
        response.raise_for_status()
    except requests.exceptions.RequestException as e:
        return {"error": f"Error accessing IUPred2A server: {str(e)}"}

    # Step 2: Parse the results to extract disorder scores
    result_lines = response.text.split("\n")
    scores = []

    for line in result_lines:
        if line.startswith("#") or not line.strip():
            continue
        parts = line.split()
        if len(parts) >= 3:
            try:
                position = int(parts[0])
                residue = parts[1]
                score = float(parts[2])
                scores.append((position, residue, score))
            except (ValueError, IndexError):
                continue

    if not scores:
        return {"error": "No valid prediction data was returned from the server."}

    # Step 3: Identify disordered regions
    disordered_regions = []
    current_region = []

    for pos, _, score in scores:
        if score >= threshold:
            if not current_region:
                current_region = [pos]
            elif pos == current_region[-1] + 1:
                current_region.append(pos)
            else:
                if len(current_region) > 1:
                    disordered_regions.append((current_region[0], current_region[-1]))
                current_region = [pos]
        elif current_region and len(current_region) > 1:
            disordered_regions.append((current_region[0], current_region[-1]))
            current_region = []
        elif current_region:
            current_region = []

    # Add the last region if it exists
    if current_region and len(current_region) > 1:
        disordered_regions.append((current_region[0], current_region[-1]))

    total_residues = len(scores)
    disordered_count = sum(1 for _, _, score in scores if score >= threshold)
    disordered_percentage = (disordered_count / total_residues) * 100 if total_residues > 0 else 0
    per_residue = [
        {
            "position": pos,
            "amino_acid": aa,
            "disorder_score": score,
            "is_disordered": score >= threshold,
        }
        for pos, aa, score in scores
    ]
    region_details = [
        {"start": start, "end": end, "length": end - start + 1} for start, end in disordered_regions
    ]

    return {
        "inputs": {
            "threshold": threshold,
            "output_file": output_file,
        },
        "summary": {
            "sequence_length": total_residues,
            "disordered_count": disordered_count,
            "disordered_percentage": disordered_percentage,
            "num_disordered_regions": len(disordered_regions),
        },
        "per_residue": per_residue,
        "disordered_regions": region_details,
        "suggested_files": {
            "per_residue_csv": output_file,
        },
        "notes": "Analysis performed using IUPred2A algorithm (long disorder mode)",
    }


def analyze_cell_morphology_and_cytoskeleton(image_path, output_dir="./results", threshold_method="otsu"):
    """Quantifies cell morphology and cytoskeletal organization from fluorescence microscopy images.

    Parameters
    ----------
    image_path : str
        Path to the fluorescence microscopy image file
    output_dir : str, optional
        Directory to save output files (default: './results')
    threshold_method : str, optional
        Method for cell segmentation ('otsu', 'adaptive', or 'manual') (default: 'otsu')

    Returns
    -------
    dict
        Dictionary containing morphology metrics, cytoskeletal metrics, and segmentation outputs.
        No files are written; output_dir is accepted for backward compatibility.

    """
    from datetime import datetime

    import cv2
    import numpy as np
    import pandas as pd
    from skimage import exposure, feature, filters, io, measure, morphology
    from skimage.color import rgb2gray

    analysis_time = datetime.now().strftime("%Y-%m-%d %H:%M:%S")

    # Load image
    try:
        image = io.imread(image_path)
        # Convert to grayscale if RGB
        if len(image.shape) > 2:
            gray_image = rgb2gray(image)
        else:
            gray_image = image

        # Enhance contrast
        gray_image = exposure.equalize_hist(gray_image)
    except Exception as e:
        return {"error": f"Error loading image: {str(e)}", "image_path": image_path}

    # Segment cells
    if threshold_method == "otsu":
        thresh = filters.threshold_otsu(gray_image)
        binary = gray_image > thresh
    elif threshold_method == "adaptive":
        binary = filters.threshold_local(gray_image, block_size=35, offset=0.05)
        binary = gray_image > binary
    else:  # manual
        thresh = 0.5  # Default value, can be parameterized
        binary = gray_image > thresh

    # Clean up binary image
    binary = morphology.remove_small_objects(binary, min_size=100)
    binary = morphology.remove_small_holes(binary, area_threshold=100)
    binary = morphology.binary_closing(binary, morphology.disk(3))
    # Label cells
    labeled_cells, num_cells = measure.label(binary, return_num=True)
    cell_props = measure.regionprops_table(
        labeled_cells,
        gray_image,
        properties=(
            "area",
            "perimeter",
            "major_axis_length",
            "minor_axis_length",
            "eccentricity",
            "orientation",
            "solidity",
        ),
    )

    # Convert to DataFrame for easier manipulation
    cell_df = pd.DataFrame(cell_props)

    # Calculate additional metrics
    if len(cell_df) > 0:
        cell_df["aspect_ratio"] = cell_df["major_axis_length"] / cell_df["minor_axis_length"]
        cell_df["circularity"] = (4 * np.pi * cell_df["area"]) / (cell_df["perimeter"] ** 2)

    # Edge detection to highlight cytoskeletal fibers
    edges = feature.canny(gray_image, sigma=2)

    # Use Hough transform to detect lines (cytoskeletal fibers)
    if np.any(edges):
        lines = cv2.HoughLinesP(
            edges.astype(np.uint8),
            1,
            np.pi / 180,
            threshold=10,
            minLineLength=10,
            maxLineGap=5,
        )

        if lines is not None:
            # Calculate line orientations
            orientations = []
            for line in lines:
                x1, y1, x2, y2 = line[0]
                if x2 - x1 != 0:  # Avoid division by zero
                    angle = np.arctan2(y2 - y1, x2 - x1) * 180 / np.pi
                    orientations.append(angle)

            if orientations:
                # Convert to numpy array for calculations
                orientations = np.array(orientations)

                # Calculate alignment metrics
                mean_orientation = np.mean(orientations)
                # Normalize angles to -90 to 90 degrees
                norm_angles = np.mod(orientations + 90, 180) - 90
                std_orientation = np.std(norm_angles)

                # Order parameter (measure of alignment, 1 = perfectly aligned, 0 = random)
                # Convert angles to radians for calculation
                rad_angles = np.radians(norm_angles)
                order_parameter = np.sqrt(np.mean(np.cos(2 * rad_angles)) ** 2 + np.mean(np.sin(2 * rad_angles)) ** 2)

                # Add fiber data to dataframe
                fiber_df = pd.DataFrame({"fiber_orientation": orientations})
            else:
                fiber_df = pd.DataFrame()
        else:
            fiber_df = pd.DataFrame()
    else:
        fiber_df = pd.DataFrame()

    morphology_summary = {
        "num_cells": int(num_cells),
    }
    if len(cell_df) > 0:
        morphology_summary.update(
            {
                "average_area": float(cell_df["area"].mean()),
                "average_aspect_ratio": float(cell_df["aspect_ratio"].mean()),
                "average_circularity": float(cell_df["circularity"].mean()),
                "average_eccentricity": float(cell_df["eccentricity"].mean()),
                "area_min": float(cell_df["area"].min()),
                "area_max": float(cell_df["area"].max()),
            }
        )

    cytoskeleton_summary = {}
    if "order_parameter" in locals():
        cytoskeleton_summary = {
            "num_fibers": int(len(orientations)) if "orientations" in locals() else 0,
            "mean_orientation_deg": float(mean_orientation),
            "std_orientation_deg": float(std_orientation),
            "order_parameter": float(order_parameter),
        }
    elif len(fiber_df) == 0:
        cytoskeleton_summary = {"num_fibers": 0}

    return {
        "inputs": {
            "image_path": image_path,
            "output_dir": output_dir,
            "threshold_method": threshold_method,
        },
        "summary": {
            "analysis_time": analysis_time,
            "morphology": morphology_summary,
            "cytoskeleton": cytoskeleton_summary,
        },
        "threshold": {
            "method": threshold_method,
            "value": float(thresh) if "thresh" in locals() else None,
        },
        "cell_morphology": cell_df.to_dict("records") if len(cell_df) > 0 else [],
        "fiber_orientations": fiber_df.to_dict("records") if len(fiber_df) > 0 else [],
        "segmentation": {
            "binary_mask": binary,
            "labeled_cells": labeled_cells,
        },
        "suggested_files": {
            "cell_morphology_csv": "cell_morphology_data.csv",
            "fiber_orientation_csv": "fiber_orientation_data.csv",
            "cell_segmentation_png": "cell_segmentation.png",
        },
    }


def analyze_tissue_deformation_flow(image_sequence, output_dir="results", pixel_scale=1.0):
    """Quantify tissue deformation and flow dynamics from microscopy image sequence.

    Parameters
    ----------
    image_sequence : list or numpy.ndarray
        Sequence of microscopy images (either a list of file paths or a 3D numpy array [time, height, width])
    output_dir : str, optional
        Directory to save results (default: "results")
    pixel_scale : float, optional
        Physical scale of pixels (e.g., μm/pixel) for proper scaling of metrics (default: 1.0)

    Returns
    -------
    dict
        Dictionary containing flow fields, deformation metrics, and summary statistics.
        No files are written; output_dir is accepted for backward compatibility.

    """
    from datetime import datetime

    import cv2
    import numpy as np

    # Load images if paths are provided
    if isinstance(image_sequence[0], str):
        loaded_images = []
        for img_path in image_sequence:
            img = cv2.imread(img_path, cv2.IMREAD_GRAYSCALE)
            if img is None:
                return {"error": f"Could not load image {img_path}", "image_path": img_path}
            loaded_images.append(img)
        frames = np.array(loaded_images)
    else:
        frames = image_sequence
        # Convert to grayscale if needed
        if len(frames.shape) > 3:  # Has color channels
            frames = np.array(
                [cv2.cvtColor(frame, cv2.COLOR_RGB2GRAY) if frame.shape[-1] == 3 else frame for frame in frames]
            )

    # Parameters for optical flow
    lk_params = {
        "winSize": (15, 15),
        "maxLevel": 2,
        "criteria": (cv2.TERM_CRITERIA_EPS | cv2.TERM_CRITERIA_COUNT, 10, 0.03),
    }

    # Metrics storage
    num_frames = len(frames)
    flow_fields = []
    divergence_maps = []
    curl_maps = []
    strain_maps = []

    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")

    # Create feature points grid (evenly spaced points)
    y, x = np.mgrid[0 : frames[0].shape[0] : 20, 0 : frames[0].shape[1] : 20]
    feature_points = np.stack((x.flatten(), y.flatten()), axis=1).astype(np.float32)

    for i in range(num_frames - 1):
        # Calculate optical flow using Lucas-Kanade
        prev_frame = frames[i]
        next_frame = frames[i + 1]

        # Compute flow for the grid points
        next_points, status, _ = cv2.calcOpticalFlowPyrLK(prev_frame, next_frame, feature_points, None, **lk_params)

        # Filter valid points
        valid_idx = status.flatten() == 1
        valid_prev_points = feature_points[valid_idx]
        valid_next_points = next_points[valid_idx]

        # Calculate displacement vectors
        displacement = valid_next_points - valid_prev_points
        points = valid_prev_points

        # Interpolate flow field to full image resolution
        flow_field = np.zeros((frames[0].shape[0], frames[0].shape[1], 2), dtype=np.float32)

        # Simple nearest-neighbor interpolation for demonstration
        # In a production system, you might use more sophisticated interpolation
        for j, (x, y) in enumerate(points.astype(int)):
            if 0 <= y < flow_field.shape[0] and 0 <= x < flow_field.shape[1]:
                flow_field[y, x] = displacement[j]

        # Calculate derivatives for deformation analysis
        u = flow_field[:, :, 0]  # x-component of flow
        v = flow_field[:, :, 1]  # y-component of flow

        # Calculate divergence (expansion/contraction)
        # Using central differences for derivatives
        u_x = cv2.Sobel(u, cv2.CV_64F, 1, 0, ksize=3) / (8.0 * pixel_scale)
        v_y = cv2.Sobel(v, cv2.CV_64F, 0, 1, ksize=3) / (8.0 * pixel_scale)
        divergence = u_x + v_y

        # Calculate curl (rotation)
        u_y = cv2.Sobel(u, cv2.CV_64F, 0, 1, ksize=3) / (8.0 * pixel_scale)
        v_x = cv2.Sobel(v, cv2.CV_64F, 1, 0, ksize=3) / (8.0 * pixel_scale)
        curl = v_x - u_y

        # Calculate strain tensor components
        strain_xx = u_x
        strain_yy = v_y
        strain_xy = 0.5 * (u_y + v_x)

        # Magnitude of strain tensor (Frobenius norm)
        strain_magnitude = np.sqrt(strain_xx**2 + strain_yy**2 + 2 * strain_xy**2)

        # Store results
        flow_fields.append(flow_field)
        divergence_maps.append(divergence)
        curl_maps.append(curl)
        strain_maps.append(strain_magnitude)

        # Optional visualization data (vector sampling)
        step = 20
        sample_vectors = []
        for y in range(0, flow_field.shape[0], step):
            for x in range(0, flow_field.shape[1], step):
                dx, dy = flow_field[y, x]
                if abs(dx) > 0.5 or abs(dy) > 0.5:
                    sample_vectors.append({"x": int(x), "y": int(y), "dx": float(dx), "dy": float(dy)})

    # Calculate summary statistics
    mean_divergence = np.mean([np.mean(div) for div in divergence_maps])
    max_divergence = np.max([np.max(div) for div in divergence_maps])
    mean_curl = np.mean([np.mean(np.abs(c)) for c in curl_maps])
    mean_strain = np.mean([np.mean(s) for s in strain_maps])

    return {
        "timestamp": timestamp,
        "inputs": {
            "output_dir": output_dir,
            "pixel_scale": pixel_scale,
        },
        "summary": {
            "num_frames": num_frames,
            "image_shape": (int(frames[0].shape[0]), int(frames[0].shape[1])),
            "mean_divergence": float(mean_divergence),
            "max_divergence": float(max_divergence),
            "mean_abs_curl": float(mean_curl),
            "mean_strain": float(mean_strain),
        },
        "results": {
            "flow_fields": flow_fields,
            "divergence_maps": divergence_maps,
            "curl_maps": curl_maps,
            "strain_maps": strain_maps,
            "sample_vectors": sample_vectors,
        },
        "suggested_files": {
            "flow_viz_png_glob": f"{output_dir}/flow_viz_*.png",
            "divergence_npy_glob": f"{output_dir}/divergence_*.npy",
            "curl_npy_glob": f"{output_dir}/curl_*.npy",
            "strain_npy_glob": f"{output_dir}/strain_*.npy",
            "summary_npy": f"{output_dir}/deformation_summary.npy",
        },
    }
