from paper.extractedData.load_csvs import *
from utils.imports import *
def plot_erosion_rings():
    data = load_unperturbed_monolayer_transcripts()
    xedges, yedges, binary_mask_cleaned, extent = compute_unperturbed_monolayer_spatial_mask(save_to_pickle=False)
    plot_binary_mask_cleaned(binary_mask_cleaned, extent)
    check_mask_fidelity(data, binary_mask_cleaned, extent)
    ring_masks = calculate_ring_masks(binary_mask_cleaned, xedges, yedges, num_iterations=NUM_ITERATIONS, plot_rings=False, save_rings=True)
    result_dict = compute_transcript_density_in_rings_all_genes(
        binary_mask=binary_mask_cleaned,
        erosion_step=5,
        num_iterations=30,
        data=data[data['name'] == 'Nupr1'],
        xedges=xedges,
        yedges=yedges,
        xy_spacing=XY_SPACING
    )
    plot_erosion_steps(ring_masks, xedges, yedges, binary_mask_cleaned, erosion_step=EROSION_STEP, num_iterations=NUM_ITERATIONS, plot_rings=True,
                       image_x_range=None, image_y_range=None)
    plot_erosion_steps(ring_masks, xedges, yedges, binary_mask_cleaned, image_x_range=[1500, 2000],
                       image_y_range=[1500, 2000])

def plot_erosion_rings_from_saved_components():
    pass

def plot_zoomed_in_erosion_rings():
    xedges, yedges, binary_mask_cleaned, extent = compute_unperturbed_monolayer_spatial_mask(save_to_pickle=False)
    ring_masks = calculate_ring_masks(binary_mask_cleaned, xedges, yedges, num_iterations=NUM_ITERATIONS,
                                      plot_rings=False, save_rings=True)
    plot_erosion_steps(ring_masks, xedges, yedges, binary_mask_cleaned, image_x_range=[1500, 2000],
                       image_y_range=[1500, 2000])
def calculate_ring_masks(binary_mask, xedges, yedges, num_iterations, plot_rings=True, save_rings=True):
    data = load_unperturbed_monolayer_transcripts()
    densities = []
    areas = []
    counts = []
    ring_masks = []

    # Get the positions of transcripts for the specified gene
    gene_data = data[data['name'] == 'Nupr1']
    x_transcripts = gene_data['x'].values
    y_transcripts = gene_data['y'].values

    # Map the transcript positions to the mask grid indices
    x_indices = np.searchsorted(xedges, x_transcripts, side='right') - 1
    y_indices = np.searchsorted(yedges, y_transcripts, side='right') - 1

    # Remove indices that are out of bounds
    valid_indices = (x_indices >= 0) & (x_indices < binary_mask.shape[1]) & \
                    (y_indices >= 0) & (y_indices < binary_mask.shape[0])

    x_indices = x_indices[valid_indices]
    y_indices = y_indices[valid_indices]

    # Initialize previous mask as the initial mask
    previous_mask = binary_mask.copy()

    for i in range(num_iterations):
        print(f"Iteration {i + 1} of {num_iterations}")
        # Erode the mask
        selem = disk(EROSION_STEP)
        eroded_mask = erosion(previous_mask, selem)

        # Compute the ring region
        ring_region = previous_mask & (~eroded_mask)

        # If the ring region is empty, break the loop
        if not ring_region.any():
            break

        # Save the ring mask for plotting
        ring_masks.append(ring_region.copy())

        # For each transcript, check if it is in the ring region
        transcript_in_ring = ring_region[y_indices, x_indices]

        # Count the number of transcripts in the ring region
        count_in_ring = np.sum(transcript_in_ring)

        # Compute the area of the ring region (assuming each pixel represents 10x10 units)
        area = ring_region.sum() * (XY_SPACING * XY_SPACING)

        # Compute the density
        density = count_in_ring / area if area > 0 else 0

        # Store the results
        counts.append(count_in_ring)
        areas.append(area)
        densities.append(density)

        # Update previous mask
        previous_mask = eroded_mask.copy()
    if save_rings:
        with open(
                r'C:\Users\micha\thesis\code\data\intestinal_organoid\non_sprinkled_july23_pasadena\monolayer_ring_masks.pkl',
                'wb') as f:
            pickle.dump(ring_masks, f)
    # Plot the rings with colors representing densities
    if plot_rings and ring_masks:
        # Create an array to hold the density values for each pixel
        density_image = np.zeros_like(binary_mask, dtype=float)

        # Normalize densities for coloring
        max_density = max(densities) if densities else 1
        norm_densities = [d / max_density for d in densities]

        # Assign density values to the pixels in each ring
        for ring_mask, density_norm in zip(ring_masks, norm_densities):
            density_image[ring_mask] = density_norm

        # Plot the density image
        plt.figure(figsize=(8, 6))
        extent = [xedges[0], xedges[-1], yedges[0], yedges[-1]]
        plt.imshow(density_image, extent=extent, origin='lower', cmap='viridis', aspect='auto')
        plt.colorbar(label='Normalized Density')
        plt.title(f'Density of Nupr1 Transcripts in Rings')
        plt.xlabel('x')
        plt.ylabel('y')
        plt.gca().set_aspect('equal', adjustable='box')
        plt.show()

    return ring_masks #densities, areas, counts


def compute_unperturbed_monolayer_spatial_mask(save_to_pickle=True):
    #Step 1: load transcripts
    data = load_unperturbed_monolayer_transcripts()


    # Step 2: Get x and y positions
    x = data['x'].values
    y = data['y'].values

    # Step 3: Create x and y bins
    edge_vec_x = np.arange(0, x.max() + XY_SPACING, XY_SPACING)
    edge_vec_y = np.arange(0, y.max() + XY_SPACING, XY_SPACING)

    # Step 4: Compute 2D histogram (note y and x are swapped to match image coordinates)
    hist, yedges, xedges = np.histogram2d(y, x, bins=[edge_vec_y, edge_vec_x])

    # Step 5: Apply Gaussian filter
    hist_filtered = gaussian_filter(hist, sigma=3)

    # Step 6: Normalize histogram
    hist_filtered_norm = hist_filtered / hist_filtered.max()

    # Step 7: Compute threshold
    threshold_density = hist_filtered_norm.min() + 0.001 * (hist_filtered_norm.max() - hist_filtered_norm.min())

    # Step 8: Create binary mask based on threshold
    binary_mask = hist_filtered_norm > threshold_density

    # Step 9: Set margins to False
    margin = 5  # Adjust as needed
    binary_mask[:margin, :] = False
    binary_mask[:, :margin] = False
    binary_mask[-margin:, :] = False
    binary_mask[:, -margin:] = False

    # Step 10: Apply morphological operations
    # Close small holes
    selem = disk(20)  # Structuring element of size 20 to close small holes
    # binary_mask_closed = binary_closing(binary_mask, selem)
    binary_mask_opened = binary_opening(binary_mask, selem)

    # Remove small objects
    binary_mask_cleaned = morphology.remove_small_objects(binary_mask_opened, min_size=200)

    extent = [xedges[0], xedges[-1], yedges[0], yedges[-1]]
    if save_to_pickle:
        save_to_pickle_monolayer_masking_components(xedges, yedges, binary_mask_cleaned, extent)

    return xedges, yedges, binary_mask_cleaned, extent
def plot_binary_mask_cleaned(binary_mask_cleaned,extent):
    plt.imshow(binary_mask_cleaned, extent=extent, origin='lower', cmap='gray', aspect='auto')
    plt.title('binary_mask_cleaned')
    plt.xlabel('x')
    plt.ylabel('y')
    plt.show()

def check_mask_fidelity(data, binary_mask_cleaned, extent):
    """
    # Here, we plot the spots of a highly expressed gene (Nupr1) on top of the binary mask to check for fidelity of the mask.
    """
    gene_name = 'Nupr1'
    gene_data = data[data['name'] == gene_name]

    # Plot the scatter plot of Nupr1 transcripts on top of the binary_mask_cleaned
    plt.figure(figsize=(8, 6))
    plt.imshow(binary_mask_cleaned, extent=extent, origin='lower', cmap='gray', aspect='auto')
    plt.scatter(gene_data['x'], gene_data['y'], s=1, c='red', alpha=0.5)
    plt.title(f'Transcripts and Mask Overlay for {gene_name}')
    plt.xlabel('x')
    plt.ylabel('y')
    plt.gca().set_aspect('equal', adjustable='box')
    plt.show()


def plot_erosion_steps(ring_masks, xedges, yedges, binary_mask, erosion_step=5, num_iterations=30, plot_rings=True,
                       image_x_range=None, image_y_range=None):
    density_image = np.zeros_like(binary_mask, dtype=float)
    color_per_ring = np.arange(len(ring_masks))
    for ring_mask, ring_color in zip(ring_masks, color_per_ring):
        density_image[ring_mask] += ring_color
    # Plot the density image
    plt.figure(figsize=(8, 6), facecolor='white')  # Set figure background to white
    extent = [xedges[0], xedges[-1], yedges[0], yedges[-1]]

    if image_x_range is not None:
        density_image = density_image[image_x_range[0]:image_x_range[1], image_y_range[0]:image_y_range[1]]
        extent = [xedges[image_x_range[0]], xedges[image_x_range[1]], yedges[image_y_range[0]],
                  yedges[image_y_range[1]]]
        binary_mask = binary_mask[image_x_range[0]:image_x_range[1], image_y_range[0]:image_y_range[1]]
    # Set the colormap and vmin so 0 values appear as white
    cmap = plt.get_cmap('hsv', len(ring_masks))
    cmap.set_under('white')
    background_image = np.where(binary_mask, 0.5, 1.0)
    plt.imshow(background_image, extent=extent, origin='lower', cmap='gray', aspect='auto', vmin=0, vmax=1)
    img = plt.imshow(density_image, extent=extent, origin='lower', cmap=cmap, aspect='auto', vmin=0.01, alpha=0.6)
    plt.title(f'Transcripts Density Rings')
    plt.xlabel('x')
    plt.ylabel('y')
    plt.gca().set_aspect('equal', adjustable='box')
    plt.tight_layout(rect=[0, 0, 0.85, 1])
    plt.show()


def compute_transcript_density_in_rings(gene_name, binary_mask, erosion_step, num_iterations, data, xedges, yedges,
                                        plot_rings=True):
    """
    Computes the density of transcripts in successive rings of the mask.
    Optionally plots the rings with colors representing densities.

    Parameters:
    - gene_name: Name of the gene (e.g., 'Nupr1')
    - binary_mask: The initial binary mask (numpy array)
    - erosion_step: The amount of erosion per iteration (in pixels)
    - num_iterations: Number of iterations
    - data: The original data DataFrame with 'x', 'y', and 'name' columns
    - xedges: The bin edges used for x (from the mask histogram)
    - yedges: The bin edges used for y (from the mask histogram)
    - plot_rings: Boolean flag to plot the rings with densities (default True)

    Returns:
    - densities: List of densities (number of transcripts per unit area) in each ring
    - areas: List of areas of each ring
    - counts: List of counts of transcripts in each ring
    """
    # Initialize lists to store results
    densities = []
    areas = []
    counts = []
    ring_masks = []

    # Get the positions of transcripts for the specified gene
    gene_data = data[data['name'] == gene_name]
    x_transcripts = gene_data['x'].values
    y_transcripts = gene_data['y'].values

    # Map the transcript positions to the mask grid indices
    x_indices = np.searchsorted(xedges, x_transcripts, side='right') - 1
    y_indices = np.searchsorted(yedges, y_transcripts, side='right') - 1

    # Remove indices that are out of bounds
    valid_indices = (x_indices >= 0) & (x_indices < binary_mask.shape[1]) & \
                    (y_indices >= 0) & (y_indices < binary_mask.shape[0])

    x_indices = x_indices[valid_indices]
    y_indices = y_indices[valid_indices]

    # Initialize previous mask as the initial mask
    previous_mask = binary_mask.copy()

    for i in range(num_iterations):
        print(f"Iteration {i + 1} of {num_iterations}")
        # Erode the mask
        selem = disk(erosion_step)
        eroded_mask = erosion(previous_mask, selem)

        # Compute the ring region
        ring_region = previous_mask & (~eroded_mask)

        # If the ring region is empty, break the loop
        if not ring_region.any():
            break

        # Save the ring mask for plotting
        ring_masks.append(ring_region.copy())

        # For each transcript, check if it is in the ring region
        transcript_in_ring = ring_region[y_indices, x_indices]

        # Count the number of transcripts in the ring region
        count_in_ring = np.sum(transcript_in_ring)

        # Compute the area of the ring region (assuming each pixel represents 10x10 units)
        area = ring_region.sum() * (XY_SPACING * XY_SPACING)

        # Compute the density
        density = count_in_ring / area if area > 0 else 0

        # Store the results
        counts.append(count_in_ring)
        areas.append(area)
        densities.append(density)

        # Update previous mask
        previous_mask = eroded_mask.copy()

    with open(
            r'C:\Users\micha\thesis\code\data\intestinal_organoid\non_sprinkled_july23_pasadena\monolayer_ring_masks.pkl',
            'wb') as f:
        pickle.dump(ring_masks, f)
    # Plot the rings with colors representing densities
    if plot_rings and ring_masks:
        # Create an array to hold the density values for each pixel
        density_image = np.zeros_like(binary_mask, dtype=float)

        # Normalize densities for coloring
        max_density = max(densities) if densities else 1
        norm_densities = [d / max_density for d in densities]

        # Assign density values to the pixels in each ring
        for ring_mask, density_norm in zip(ring_masks, norm_densities):
            density_image[ring_mask] = density_norm

        # Plot the density image
        plt.figure(figsize=(8, 6))
        extent = [xedges[0], xedges[-1], yedges[0], yedges[-1]]
        plt.imshow(density_image, extent=extent, origin='lower', cmap='viridis', aspect='auto')
        plt.colorbar(label='Normalized Density')
        plt.title(f'Density of {gene_name} Transcripts in Rings')
        plt.xlabel('x')
        plt.ylabel('y')
        plt.gca().set_aspect('equal', adjustable='box')
        plt.show()

    return densities, areas, counts


def compute_transcript_density_in_rings_all_genes(binary_mask, erosion_step, num_iterations, data, xedges, yedges,
                                                  xy_spacing=XY_SPACING):
    """
    Computes the density of transcripts in successive rings of the mask for all genes.

    Parameters:
    - binary_mask: The initial binary mask (numpy array)
    - erosion_step: The amount of erosion per iteration (in pixels)
    - num_iterations: Number of iterations
    - data: The original data DataFrame with 'x', 'y', and 'name' columns
    - xedges: The bin edges used for x (from the mask histogram)
    - yedges: The bin edges used for y (from the mask histogram)
    - xy_spacing: The spacing between x and y coordinates (used to compute area)

    Returns:
    - result_dict: A dictionary where keys are gene names, and values are dictionaries with keys 'densities', 'areas', 'counts'
    """
    result_dict = {}

    # Get unique gene names
    gene_names = data['name'].unique()

    # Map transcript positions to indices for all genes
    gene_indices = {}

    for i, gene_name in enumerate(gene_names):
        print(f"Preprocessing gene: {gene_name}, gene number {i} out of {len(gene_names)}")
        # Initialize result storage for each gene
        result_dict[gene_name] = {'densities': [], 'areas': [], 'counts': []}

        # Get the positions of transcripts for the gene
        gene_data = data[data['name'] == gene_name]
        x_transcripts = gene_data['x'].values
        y_transcripts = gene_data['y'].values

        # Map the transcript positions to the mask grid indices
        x_indices = np.searchsorted(xedges, x_transcripts, side='right') - 1
        y_indices = np.searchsorted(yedges, y_transcripts, side='right') - 1

        # Remove indices that are out of bounds
        valid_indices = (
                (x_indices >= 0) & (x_indices < binary_mask.shape[1]) &
                (y_indices >= 0) & (y_indices < binary_mask.shape[0])
        )
        x_indices = x_indices[valid_indices]
        y_indices = y_indices[valid_indices]

        # Store indices for each gene
        gene_indices[gene_name] = (x_indices, y_indices)

    # Initialize previous mask as the initial mask
    previous_mask = binary_mask.copy()

    # Create structuring element once (fixed size)
    selem = disk(erosion_step)

    for i in range(num_iterations):
        print(f"Iteration {i + 1} of {num_iterations}")

        # Erode the mask
        eroded_mask = erosion(previous_mask, selem)

        # Compute the ring region
        ring_region = previous_mask & (~eroded_mask)
        # if num_iterations % 10 == 0:
        #     ring_mask_uint8 = (ring_region * 255).astype(np.uint8)
        #     find_ring_width_in_pixels(ring_mask_uint8)
        #     plt.imshow(ring_region)
        #     plt.show()

        # If the ring region is empty, break the loop
        if not ring_region.any():
            print("Ring region is empty. Stopping iterations.")
            break

        # Compute the area of the ring region
        area = ring_region.sum() * (xy_spacing ** 2)

        # Loop over all genes
        for gene_name in gene_names:
            x_indices, y_indices = gene_indices[gene_name]

            # For each transcript, check if it is in the ring region
            transcript_in_ring = ring_region[y_indices, x_indices]

            # Count the number of transcripts in the ring region
            count_in_ring = np.sum(transcript_in_ring)

            # Compute the density
            density = count_in_ring / area if area > 0 else 0

            # Store the results
            result_dict[gene_name]['counts'].append(count_in_ring)
            result_dict[gene_name]['areas'].append(area)
            result_dict[gene_name]['densities'].append(density)

        # Update previous mask
        previous_mask = eroded_mask.copy()

    return result_dict