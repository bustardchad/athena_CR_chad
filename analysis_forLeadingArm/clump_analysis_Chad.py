"""All to do with clump / surface finding.
"""
import sys, os
import numpy as np
from skimage import measure
import scipy.ndimage
from scipy.ndimage import gaussian_filter

from joblib import Parallel, delayed

def _march_cubes(args):
    if len(args) == 4:
        smoothing = True
    elif len(args) == 3:
        smoothing = False
    else:
        raise ValueError("Not recognized length of argumments: %d" %(len(args)))

    if smoothing:
        dat, threshold, dgrid, sigma = args
        if sigma > 0:
            dat = gaussian_filter(dat, sigma)
    else:
        dat, threshold, dgrid = args

    if threshold > np.max(dat):
        return 0.0
    if threshold < np.min(dat):
        return np.inf # np.product(np.shape(dat) * spacing)
    verts, faces, normals, values = measure.marching_cubes_lewiner(dat,
                                                                   level = threshold,
                                                                   spacing = dgrid)
    area = measure.mesh_surface_area(verts, faces)

    return area


def compute_surface_area(ds, scale_factors, threshold = 2.0, mode = 'density', ncpus = 1,
                         verbose = 0, sigma_smooth = None, **kwargs):
    """Computes surface of cold gas.

    Keyword Arguments:
    ds        -- DataSet
    scale_factors -- List of scale factors to multiply boxsize by (in x,y,z direction).
    mode      -- Mode to where to take the threshold.
                 Can be `temp` (default), `density` or any other field...
    threshold -- Threshold in units of *_cl (default 2.0), can be list
    sigma_smooth -- Width of Gaussian kernel (in cell lengths), can be list.

    **kwargs  -- passed to `measure.marching_cubes_lewiner`

    Returns:
    Surface area in units of r_cl**2
    """
    assert len(scale_factors) == 3, "`scale_factors` should be list of 3."
    scale_factors = np.asarray(scale_factors)

    #assert mode in ['temp', 'density', 't_cool'], "Unknown mode `%s`" %(mode)
    if not isinstance(threshold, (list, tuple, np.ndarray)):
        threshold = [threshold]

    #ad = ds.all_data()
    #T = np.reshape(ad['temp'], ds.domain_dimensions)
    g = ds.covering_grid(level=0, left_edge=ds.domain_left_edge,
                         dims=ds.domain_dimensions)
    if mode == 'temp':
        dat = g['temp'].to('T_cl').value
    elif mode == 'density':
        dat = g['density'].value
    else:
        dat = g[mode].value


    dgrid = tuple((scale_factors * ds.domain_width / ds.domain_dimensions).value)

    if sigma_smooth is not None:
        assert len(threshold) == 1, "Both sigma & different thresholds not supported yet."
        th = threshold[0]
        n = np.min([len(sigma_smooth), ncpus])
        if verbose > 0:
            print("Using %d CPUs" %(n))
        areas = Parallel(n_jobs=n, verbose=verbose)(delayed(_march_cubes)((dat, th, dgrid, sigma))
                                                    for sigma in sigma_smooth)
    else:
        n = np.min([len(threshold), ncpus])
        areas = Parallel(n_jobs=n,verbose=verbose)(delayed(_march_cubes)((dat, th, dgrid))
                                                   for th in threshold)
    if len(areas) == 1:
        return areas[0]
    print(areas)
    return areas


def test_compute_surface_area(fn):
    """Returns surface area.
    """
    import yt
    import yt_custom
    import athena_input as ai
    print("Computing surface area of %s" %(fn))
    ds = yt.load(fn)
    yt_custom.add_fields(ds)
    cfg = ai.read(os.path.dirname(fn) + "/../input")
    yt_custom.add_units(ds, cfg)

    A = compute_surface_area(ds)
    print("Area: %g r_cl^2" %(A))


def count_clumps(ds, threshold, mode = 'density'):
    """Simple way to count clumps. Returns number of clumps found.

    Keyword Arguments:
    ds             -- yt dataset
    threshold      -- Density threshold in code units. Can be list to evaluate several.
    """
    if not isinstance(threshold, (list, tuple, np.ndarray)):
        threshold = [threshold]


    # construct array of values
    g = ds.covering_grid(level=0, left_edge=ds.domain_left_edge,
                         dims=ds.domain_dimensions)
    dat = g[mode].value

    num_features = []
    for cth in threshold:
        if mode == 'temp':
            blobs = dat < cth
        else:
            blobs = dat > cth
        la, nf = scipy.ndimage.measurements.label(blobs)
        num_features.append(nf)

    if len(num_features) == 1:
        return num_features[0]
    return num_features

 
def test_clump_counting_2d():
    from skimage import measure
    from skimage import filters
    import matplotlib.pyplot as plt
    import numpy as np

    n = 12
    l = 256
    np.random.seed(1)
    im = np.zeros((l, l))
    points = l * np.random.random((2, n ** 2))
    im[(points[0]).astype(np.int), (points[1]).astype(np.int)] = 1
    im = filters.gaussian(im, sigma= l / (4. * n))
    blobs = im > 0.7 * im.mean()

    all_labels = measure.label(blobs)
    blobs_labels = measure.label(blobs, background=0)

    plt.figure(figsize=(9, 3.5))
    plt.subplot(131)
    plt.imshow(blobs, cmap='gray')
    plt.axis('off')
    plt.subplot(132)
    plt.imshow(all_labels, cmap='nipy_spectral')
    plt.axis('off')
    plt.subplot(133)
    plt.imshow(blobs_labels, cmap='nipy_spectral')
    plt.axis('off')

    plt.tight_layout()
    plt.show()



def test_clump_counting_3d():
    from skimage import measure
    from skimage import filters
    import matplotlib.pyplot as plt
    import numpy as np
    import scipy

    n = 12
    l = 256
    np.random.seed(1)
    im = np.zeros((l, l,l))
    points = np.random.random(size=(l,l,l))
    #im[(points[0]).astype(np.int), (points[1]).astype(np.int)] = 1
    im = points
    im = filters.gaussian(im, sigma= l / (50. * n))
    blobs = im > 0.7 * im.mean()

    labeled_array, num_features = scipy.ndimage.measurements.label(blobs)

    print("%d features found." %(num_features))

    plt.figure(figsize=(9, 3.5))
    plt.subplot(131)
    plt.imshow(blobs[:,:,0], cmap='gray')
    plt.axis('off')
    plt.subplot(132)
    plt.imshow(labeled_array[:,:,0], cmap='nipy_spectral')
    plt.axis('off')
    plt.subplot(133)
    plt.imshow(labeled_array[:,:,l/2], cmap='nipy_spectral')
    plt.axis('off')

    plt.tight_layout()
    plt.show()






if __name__ == '__main__':
    # test_compute_surface_area(sys.argv[1])

    test_clump_counting_3d()


