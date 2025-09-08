import numpy as np

from paper.extractedData.load_csvs import *
from utils.imports import *
def calculate_eroded_transcription_densities(plot_erosion=True, save_components=True):
    """
    Thie function calculates the erosion components for measuring transcript densities in rings from the monolayer
    edge, inwards. These transcript profiles
    :param plot_erosion: boolean variable signifying if to plot the erosion rings
    :param save_components: boolean variable signifying if to save the erosion components
    :return:  result_dict: a dictionary holding the transcription profiles of the genes across the erosion rings
    used for downstream analysis and comparison to in vivo villus zonation expression profiles
    #TODO add the result dict structure more in detail
    """
    data = load_unperturbed_monolayer_transcripts()
    xedges, yedges, binary_mask_cleaned, extent = compute_unperturbed_monolayer_spatial_mask(save_to_pickle=save_components)
    #view the binary mask of the full monolayer that with be eroded, uncomment to plot
    #plot_binary_mask_cleaned(binary_mask_cleaned, extent)
    #sanity check to view the fidelity of the mask, uncomment to plot
    #check_mask_fidelity(data, binary_mask_cleaned, extent)
    ring_masks, avg_ring_width, densities, areas, counts = calculate_ring_masks(binary_mask_cleaned, xedges, yedges, num_iterations=NUM_ITERATIONS, plot_rings=False, save_rings=save_components)
    result_dict = compute_transcript_density_in_rings_all_genes(
        binary_mask=binary_mask_cleaned,
        erosion_step=5,
        num_iterations=30,
        data=data[data['name'] != 'GFP'],
        xedges=xedges,
        yedges=yedges,
        xy_spacing=XY_SPACING, save=save_components)

    if plot_erosion:
        plot_erosion_rings_from_calculated_components(xedges, yedges, binary_mask_cleaned, ring_masks)
    return result_dict

def plot_erosion_rings_from_calculated_components(xedges, yedges, binary_mask_cleaned, ring_masks):
    """
    This function plots erosion rings , can be called directly following their calculation
    """
    #plots the erosion rings over the full monolayer
    plot_erosion_steps(ring_masks, xedges, yedges, binary_mask_cleaned, erosion_step=EROSION_STEP, num_iterations=NUM_ITERATIONS, plot_rings=True,
                       image_x_range=None, image_y_range=None)
    #plot the erosion on a zoomed in area
    plot_erosion_steps(ring_masks, xedges, yedges, binary_mask_cleaned, image_x_range=EROSION_RINGS_ZOOM_IN,
                       image_y_range=EROSION_RINGS_ZOOM_IN)

def plot_erosion_rings_from_saved_components():
    """
    This function plots the erosion rings in which transcript density of each of the paneled genes
    were measured. It can be used to plot the rings if the erosion components have already been calculated and saved
    """

    #load the components necessary for plotting the erosion rings
    xedges, yedges, binary_mask_cleaned, ring_masks, extent = load_erosion_components()

    #plot full monolayer erosion rings
    plot_erosion_steps(ring_masks, xedges, yedges, binary_mask_cleaned, erosion_step=EROSION_STEP,
                       num_iterations=NUM_ITERATIONS, plot_rings=True,
                       image_x_range=None, image_y_range=None)
    #plot erosion rings on zoomed in region
    plot_erosion_steps(ring_masks, xedges, yedges, binary_mask_cleaned, image_x_range=EROSION_RINGS_ZOOM_IN,
                       image_y_range=EROSION_RINGS_ZOOM_IN)


def calculated_and_plot_zoomed_in_erosion_rings():
    """
    This function computes the erosion components for measuring transcript density in rings from
    monolayer edge inwards and plots the erosion rings on a zoomed in region
    """
    #compute the erosion components
    xedges, yedges, binary_mask_cleaned, extent = compute_unperturbed_monolayer_spatial_mask(save_to_pickle=False)
    #compute ring masks
    ring_masks, densities, areas, counts = calculate_ring_masks(binary_mask_cleaned, xedges, yedges, num_iterations=NUM_ITERATIONS,
                                      plot_rings=False, save_rings=True)
    #plot the erosion steps on a zoomed in region
    plot_erosion_steps(ring_masks, xedges, yedges, binary_mask_cleaned, image_x_range=EROSION_RINGS_ZOOM_IN,
                       image_y_range=EROSION_RINGS_ZOOM_IN)

def calculate_ring_masks(binary_mask, xedges, yedges, num_iterations, plot_rings=True, save_rings=True,output_path = AUTONOMOUS_ZONATION_PLOTS_FOLDER_PATH):
    """
    This function calculates the ring masks in which in ring, the transcript densities for each gene are measured.
    :param plot_rings: boolean: whether to plot the ring masks
    :param save_rings: boolean: whether to save the ring masks
    :return: returns the ring masks and the average ring width per iteration
    """
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

    avg_ring_widths = []
    for i in range(num_iterations):
        print(f"Iteration {i + 1} of {num_iterations}")
        # Erode the mask
        selem = disk(EROSION_STEP)
        eroded_mask = erosion(previous_mask, selem)
        # Compute the ring region
        ring_region = previous_mask & (~eroded_mask)
        avg = calculate_iteration_width(previous_mask, eroded_mask)
        avg_ring_widths.append(avg)
        df = pd.DataFrame({'avg_ring_width': avg_ring_widths})
        df.to_csv(os.path.join(output_path,"avg_ring_widths.csv"), index=False)
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
    arr_avg = np.array(avg_ring_widths)
    print(np.mean(arr_avg[np.isfinite(avg_ring_widths)]))

    if save_rings:
        save_one_monolayer_masking_component_to_pickle(ring_masks, 'monolayer_ring_masks')
        save_one_monolayer_masking_component_to_pickle(arr_avg[np.isfinite(avg_ring_widths)],'iteration_widths')

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
        plt.close()

    print(avg_ring_widths)
    arr_avg = np.array(avg_ring_widths)
    print(np.mean(arr_avg[np.isfinite(avg_ring_widths)]))
    return ring_masks, np.array(avg_ring_widths) ,densities, areas, counts

def calculate_iteration_width(mask1, mask2):
    """
    This function calculates the average width in micrometers of one ring
    """
    ring_mask = mask1 & ~mask2
    labeled_ring = label(ring_mask)
    labeled_outer = label(mask1)

    # Find the properties of the labeled regions
    region_ring = regionprops(labeled_ring)[0]
    region_mask = regionprops(labeled_outer)[0]

    # Area and perimeter of the ring
    ring_area = region_ring.area
    ring_perimeter = region_mask.perimeter

    # Correct the area and perimeter back to original dimensions
    corrected_ring_area = ring_area*(XY_SPACING ** 2) #binned and spaced every XY_SPACING beforehand in reality,
    corrected_ring_perimeter = ring_perimeter*(XY_SPACING)
        # guard against 0 perimeter or 0 area
    if corrected_ring_perimeter <= 0 or corrected_ring_area <= 0:
        return 0
    # Calculate the average width of the ring based on the corrected perimeter
    average_ring_width = corrected_ring_area / corrected_ring_perimeter

    # Convert to micrometers
    average_ring_width_um = average_ring_width * PIXEL2NM / 1000
    return average_ring_width_um


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
    plt.close()

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
    plt.close()


def plot_erosion_steps(ring_masks, xedges, yedges, binary_mask, erosion_step=5, num_iterations=30, plot_rings=True,
                       image_x_range=None, image_y_range=None, output_path = AUTONOMOUS_ZONATION_PLOTS_FOLDER_PATH):
    density_image = np.zeros_like(binary_mask, dtype=float)
    color_per_ring = np.arange(len(ring_masks))
    for ring_mask, ring_color in zip(ring_masks, color_per_ring):
        density_image[ring_mask] += ring_color
    extent = [xedges[0], xedges[-1], yedges[0], yedges[-1]]
    fig, ax = plt.subplots(figsize=(8, 8))
    if image_x_range is not None:
        density_image = density_image[image_x_range[0]:image_x_range[1], image_y_range[0]:image_y_range[1]]

        extent = [xedges[image_x_range[0]], xedges[image_x_range[1]], yedges[image_y_range[0]],
                  yedges[image_y_range[1]]]
        binary_mask = binary_mask[image_x_range[0]:image_x_range[1], image_y_range[0]:image_y_range[1]]
        # min_x =10
        # min_y = 700
        min_x = extent[0]  # Use updated extent for zoomed-in case
        min_y = extent[2]


    else:
        min_x = 1000
        min_y = 3000

    # Create figure and axis

    cmap = plt.get_cmap('hsv', len(ring_masks))
    cmap.set_under('white')
    background_image = np.where(binary_mask, 0.5, 1.0)
    flipped_image = np.flipud(np.where(binary_mask, 0.5, 1.0))
    ax.imshow(flipped_image, extent=extent, origin='lower', cmap='gray', aspect='auto', vmin=0,
              vmax=1)
    flipped_density_image = np.flipud(density_image)
    img = ax.imshow(flipped_density_image, extent=extent, origin='lower', cmap=cmap, aspect='auto', vmin=0.01, alpha=0.8)

    # Titles and labels
    ax.set_title('Transcripts Density Rings', fontsize=14)
    ax.set_xlabel('x', fontsize=12)
    ax.set_ylabel('y', fontsize=12)
    avg_ring_widths = pd.read_csv(os.path.join(output_path,"avg_ring_widths.csv"))['avg_ring_width'].values
    legend_patches = [
    mpatches.Patch(
        color=cmap(i / len(avg_ring_widths)),
        label=f'{int(avg_ring_widths[:i].sum())} \u03BCm'
    )
    for i, width in enumerate(avg_ring_widths) if width > 0]


    legend_patches.append(mpatches.Patch(color='grey', label='Monolayer'))

    plt.legend(handles=legend_patches, title="Ring Colors", loc='center left', bbox_to_anchor=(1, 0.5), borderaxespad=0, fontsize=12)

    # Add scale bar manually

    # scale_bar_length = 100  # Scale bar length in micrometers
    # pixel_size = 107.11 # Pixel size in nanometers
    # scale_bar_length_nm = scale_bar_length *(1000) #1000 for nanometer, but everything is ten times larger in each axis than a pixel already
    # scale_bar_length_pixels = scale_bar_length_nm / pixel_size  # Convert to pixels
    #
    # # Position the scale bar
    # scale_bar_x_start = 0.1  # Fraction of the width from the left
    # scale_bar_y_pos = 0.05  # Fraction of the height from the bottom
    # bar_start_x = extent[0] + scale_bar_x_start * (extent[1] - extent[0])
    # bar_end_x = bar_start_x + scale_bar_length_pixels * (extent[1] - extent[0]) / density_image.shape[1]
    # bar_y = extent[2] + scale_bar_y_pos * (extent[3] - extent[2])
    #
    # # Plot scale bar
    # ax.plot([bar_start_x, bar_end_x], [bar_y, bar_y], color='black', linewidth=3, solid_capstyle='butt')
    #
    #     # Add scale bar label
    # ax.text((bar_start_x + bar_end_x) / 2, bar_y - 0.02 * (extent[3] - extent[2]),
    #         f'{scale_bar_length} µm', color='black', fontsize=12, ha='center', va='top')
    scale_bar_length_um = 100
    pixel_size_nm = 107.11  # Pixel size in nanometersscale_bar_length_um = 100
    scale_bar_length_px = int(scale_bar_length_um * 1000.0 / pixel_size_nm)  # px
    # 8) Decide where "bottom-left" is visually.
    #    Because we inverted the y-axis:
    #       - left = min_x
    #       - bottom = max_y
    #    We'll offset by some margin from these edges.
    offset_x = 100
    offset_y = 100

    scale_bar_height = 200 #scale_bar_length_px / 4.665
    scale_bar_color = 'black'
    bar_left = min_x + offset_x
    # For the bottom-left visually, place the rectangle's "bottom" near max_y minus offset_y:
    bar_bottom = min_y - offset_y - scale_bar_height
    # 9) Place the scale bar rectangle
    rect = plt.Rectangle((extent[0] + 100, extent[2] + 100),  # Bottom-left corner
                         scale_bar_length_px,  # Width
                         200,  # Height
                         color='black', lw=0)
    ax.add_patch(rect)


    if image_x_range is not None:
        file_name = os.path.join(AUTONOMOUS_ZONATION_PLOTS_FOLDER_PATH,
                             'erosion_rings_zoom_in.pdf' if image_x_range else 'erosion_rings_zoom_in.pdf')
    else:
        file_name = os.path.join(AUTONOMOUS_ZONATION_PLOTS_FOLDER_PATH,
                                 'erosion_rings_zoom_in.pdf' if image_x_range else 'erosion_rings_full_monolayer.pdf')
    # Remove axis ticks and save
    ax.axis('off')
    plt.tight_layout()
    os.makedirs(AUTONOMOUS_ZONATION_PLOTS_FOLDER_PATH, exist_ok=True)
    plt.savefig(file_name, format='pdf', bbox_inches='tight')
    plt.show()
    plt.close()

def compute_transcript_density_in_rings(gene_name, binary_mask, erosion_step, num_iterations, data, xedges, yedges,
                                        plot_rings=True,output_path = AUTONOMOUS_ZONATION_PLOTS_FOLDER_PATH):
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
    avg_ring_widths = []
    for i in range(num_iterations):
        print(f"Iteration {i + 1} of {num_iterations}")
        # Erode the mask
        selem = disk(erosion_step)
        eroded_mask = erosion(previous_mask, selem)

        # Compute the ring region
        ring_region = previous_mask & (~eroded_mask)
        avg = calculate_iteration_width(previous_mask, eroded_mask)
        avg_ring_widths.append(avg)
        widths = np.asarray(avg_ring_widths, dtype=float)
        # Replace inf with NaN, then interpolate; finally fill any edge NaNs with 0
        s = pd.Series(widths).replace([np.inf, -np.inf], np.nan)
        s = s.interpolate(limit_direction="both").fillna(0.0)

        pd.DataFrame({'avg_ring_width': s.values}).to_csv(
            os.path.join(output_path, "avg_ring_widths.csv"),
            index=False)

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

    save_one_monolayer_masking_component_to_pickle(ring_masks,'monolayer_ring_masks')


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
        plt.close()
    print(avg_ring_widths)
    print(np.mean(avg_ring_widths))
    return np.array(avg_ring_widths), densities, areas, counts


def compute_transcript_density_in_rings_all_genes(binary_mask, erosion_step, num_iterations, data, xedges, yedges,
                                                  xy_spacing=XY_SPACING, save=False):
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
    if save:
        save_transcript_densities_unperturbed_monolayer(result_dict)
    return result_dict


def plot_density_profiles(result_dict, gene_names=None, normalize=False, spread_plots=False):
    """
    Plots the density profiles for multiple genes.

    Parameters:
    - result_dict: The dictionary returned by compute_transcript_density_in_rings_all_genes function.
    - gene_names: List of gene names to plot. If None, plots all genes.
    - normalize: Boolean flag to normalize densities by their maximum value (default False).
    - spread_plots: Boolean flag to create stacked plots for each gene (default False).

    Returns:
    - None
    """
    import numpy as np
    import os
    import matplotlib.pyplot as plt

    # Window in physical units (µm)
    LO, HI = 15, 70

    # If gene_names is None, plot all genes
    if gene_names is None:
        gene_names = list(result_dict.keys())

    if not gene_names:
        print("No genes to plot.")
        return

    # helper: distance axis in µm given number of rings
    avg_iteration_width = load_iteration_average_width()
    def _dist_axis(n):
        return np.round(np.arange(1, n + 1) * avg_iteration_width).astype(int)

    if spread_plots:
        # Create stacked subplots with shared X-axis
        num_genes = len(gene_names)
        fig, axs = plt.subplots(num_genes, 1, sharex=True, figsize=(9, 2.5*num_genes))
        if num_genes == 1:
            axs = [axs]  # Ensure axs is a list even if there is only one subplot

        yticks = [0, 0.5, 1]
        for idx, (ax, gene_name) in enumerate(zip(axs, gene_names)):
            if gene_name not in result_dict:
                print(f"Gene '{gene_name}' not found in the results.")
                continue

            # take profile and trim last 10 points as in your original code
            dens = result_dict[gene_name].get('densities', [])[:-10]
            if not dens:
                print(f"No density data available for gene '{gene_name}'.")
                continue

            y = np.asarray(dens, dtype=float)
            x = _dist_axis(y.size)

            # window to LO..HI µm
            mask = (x >= LO) & (x <= HI)
            xw, yw = x[mask], y[mask]
            if xw.size == 0:
                print(f"Gene '{gene_name}': no points in {LO}–{HI} µm.")
                continue

            # normalize AFTER slicing to the window (matches your intent)
            if normalize:
                yw = yw-min(yw)  # shift to zero min
                m = float(yw.max())
                if m > 0:
                    yw = yw / m
                else:
                    yw = yw  # leave as is (all zeros)

            ax.plot(xw, yw, marker='o', label=gene_name, linewidth=5)
            ax.set_ylabel('Normalized\n Density' if normalize else 'Density', fontsize=20)
            # Keep your original y-ticks choice
            ax.set_yticks(yticks, yticks)
            ax.tick_params(axis="y", labelsize=20)

            
            # Optional: nice x ticks every 4 µm (comment out if you prefer default)
            ax.set_xticks(np.arange(LO, HI + 1, 4), np.arange(LO, HI + 1, 4), rotation=45, fontsize=16)


            ax.grid(True)
            if idx == 0:
                ax.set_title('Density Profiles of Genes in Successive Rings', fontsize=20)
            ax.legend(loc='lower right', fontsize=18)

        axs[-1].set_xlabel('Distance to Monolayer Edge (μm)', fontsize=20)

        # === EXACT save/show logic & path as your original ===
        plt.tight_layout()
        os.makedirs(AUTONOMOUS_ZONATION_PLOTS_FOLDER_PATH, exist_ok=True)
        file_name = os.path.join(
            AUTONOMOUS_ZONATION_PLOTS_FOLDER_PATH,
            'monolayer_zonation_expression_profiles.pdf'
        )
        plt.savefig(file_name, format='pdf')
        plt.show()
        plt.close()

    else:
        # Plot all genes on the same plot (no saving in original)
        plt.figure(figsize=(12, 8))
        plotted_any = False

        for gene_name in gene_names:
            if gene_name not in result_dict:
                print(f"Gene '{gene_name}' not found in the results.")
                continue

            dens = result_dict[gene_name].get('densities', [])[:-10]
            if not dens:
                print(f"No density data available for gene '{gene_name}'.")
                continue

            y = np.asarray(dens, dtype=float)
            x = _dist_axis(y.size)

            # window
            mask = (x >= LO) & (x <= HI)
            xw, yw = x[mask], y[mask]
            if xw.size == 0:
                continue

            if normalize:
                m = float(yw.max())
                if m > 0:
                    yw = yw / m

            plt.plot(xw, yw, marker='o', label=gene_name)
            plotted_any = True

        plt.title('Density Profiles of Genes in Successive Rings')
        plt.xlabel('Distance to Monolayer Edge (μm)')
        plt.ylabel('Normalized Density' if normalize else 'Density (transcripts per unit area)', fontsize=20)
        plt.grid(True)
        if plotted_any:
            plt.xlim(LO, HI)
            plt.legend()
        plt.show()
        plt.close()

def plot_unperturbed_monolayer_gene_density_to_invivo_comparisons():
    """
    This function compares expression profiles measured in the unpderturbed monolayer by measuring in eroded rings the
    transcript densities in comparison to in vivo intestine expression measured in Moor. et al. 'Spatial Reconstruction
    of Single Enterocytes Uncovers Broad Zonation along the Intestinal Villus Axis' 2018

    The plot for each gene comparison is saved in the same directory
    """
    result_dict = load_unperturbed_monolayer_gene_densities()
    invivo_exp_raw = load_TPM_LCM_intestine_atlas()
    invivo_exp,genes_LCM =  get_LCM_atlas_gene_subset(invivo_exp_raw, ORGANOID_GENE_NAMES)
    gene_density_df = pd.DataFrame({gene: values['densities'] for gene, values in result_dict.items()})[5:20][::-1][genes_LCM]
    gene_density_smoothened = gene_density_df.apply(apply_savgol)
    gene_density_df_normalized = gene_density_smoothened.apply(lambda x: (x - x.min()) / (x.max() - x.min()))
    x_gene_density = np.linspace(0,1,gene_density_df.shape[0])
    x_invivo = np.linspace(0, 1, invivo_exp.shape[1])
    save_path = os.path.join(UNPERTURBED_DIR, 'invivo_to_monolayer_expression_profiles')

    for i, gene in enumerate(genes_LCM):
        print(f'{gene}, num {i}')
        plt.figure(figsize=(8, 8))
        plt.plot(x_gene_density, gene_density_df_normalized[gene], label='enteroid monolayer')
        plt.plot(x_invivo, invivo_exp.loc[gene], label='in vivo villus')
        plt.legend()
        plt.xlabel('bottom to top villus axis')
        plt.ylabel('expression')
        plt.title(f'{gene} invivo to enteroid monolayer expression profile comparison')
        plt.savefig(os.path.join(save_path, f'{gene}_invivo_monolayer_comparison.pdf'),bbox_inches='tight')
        #plt.show()
        plt.close()

def apply_savgol(column):
    return savgol_filter(column, window_length=15, polyorder=3)