"""Manual clump analysis. Does not take yt clump finder as too slow.
"""
import scipy.ndimage.measurements
from skimage import measure
import scipy.ndimage

#from tqdm import tqdm
import numpy as np
import os
import pandas
import logging
import psutil

def _asel(a,ia,pad):
    ia = np.array(ia, dtype = int)
    # TODO: be smarter about boundaries...
    js = np.clip(ia[:,0]-pad,0,a.shape)
    je = np.clip(ia[:,1]+pad,0,a.shape)
    # return a[(ia[0][0]-pad):(ia[0][1]+pad),(ia[1][0]-pad):(ia[1][1]+pad),(ia[2][0]-pad):(ia[2][1]+pad)]
    return a[js[0]:je[0],js[1]:je[1],js[2]:je[2]]


def _analyze_clump(cid, la, g, comi, size, tempthresh_hot, mode, thresh_unit = None,
                   padding_hot = 2, mincells_hot = 16, compute_surface_area = True):
    """
    Keyword Arguments:
    cid -- Clump id
    la  -- array of clump IDs
    g   -- yt grid
    comi -- Center of mass in index units
    padding_hot  -- (default 2)
    mincells_hot -- (default 16)
    """
    assert 'temp' in mode, "mode `%s` not supported yet." %(mode)

    cdat = {'cid' : cid, 'isize' : size,
            'ix' : comi[0], 'iy' : comi[1], 'iz' : comi[2] }
    #ia = [int(i) for i in np.round(np.array(comi))]

    ia = []
    for i in range(3):
        ii = np.arange(la.shape[i])
        myax = tuple([ j for j in range(3) if i != j ])
        cm = np.any(la == cid, axis = myax)
        ia.append([ii[cm][0], ii[cm][-1]])

    #pad = int(np.round(0.5 * size**(1/3.))) + padding_hot
    pad = padding_hot
    while True:
        if thresh_unit is None:
            cgrid = g[mode].value
        else:
            cgrid = g[mode].to(thresh_unit).value
        mhot = _asel(cgrid,ia,pad) > tempthresh_hot
        if np.sum(mhot) > mincells_hot:
            break
        pad += 1

    mclump = _asel(la,ia,pad) == cid
    assert np.sum(mclump) == size, "%d versus %d (is = %s)" %(np.sum(mclump), size, str(ia))

    # Compute projected areas
    for iaxis in range(3):
        cdat['iarea%d' %(iaxis)] = np.sum(np.any(mclump,axis = iaxis))

    # Compute surface area
    if compute_surface_area is True:
        verts, faces, normals, values = measure.marching_cubes_lewiner(mclump,
                                                                       level = 0.5)
        cdat['isurface'] = measure.mesh_surface_area(verts, faces)


    # Compute quantities in and around clump
    cdat['ipad'] = pad

    for ck in [mode, 'density', 'velocity_x', 'velocity_y', 'velocity_z']:
        for cm, clbl in [(mhot, 'hot'), (mclump, 'clump')]:
            cdat[ck + '_' + clbl] = float(np.mean(_asel(g[ck],ia,pad)[cm]).value)
        if psutil.virtual_memory().percent > 93.0:
            logging.warning("Clearing memory.")
            g.clear_data() # to free some memory-->can't do that as will be used for other clumps

    return cdat


def _compute_clumps_to_particles(ds, la):
    ad = ds.all_data()
    pos = np.vstack([ad['particle_pos%s' %(i)].value for i in 'xyz']).T
    dx = ds.domain_width.value / ds.domain_dimensions
    ipos = np.array(np.floor(pos/dx),dtype=int)
    iblobs = la[ipos[:,0],ipos[:,1],ipos[:,2]]

    pids = np.array(ad['particle_tag'].value,dtype=int)

    mydat= np.array([pids,iblobs]).T
    df_particles = pandas.DataFrame(data=mydat,columns=['pid','clump']).set_index('pid')
    return df_particles


def clump_analysis(ds, threshold = None, thresh_hot = None, thresh_unit = None,
                   outdir = None, mode = 'temperature', minsize = 16,
                   clumps_to_particles = False,
                   **kwargs):
    """
    Keyword Arguments:
    ds           -- yt data set
    threshold    -- Threshold in code units.
    mode         -- (default 'temp')
    minsize      -- (default 16)
    clumps_to_particles -- If `True` will also return DF in which particle ID &
                     corresponding blob ID is listed
    """
    assert threshold is not None
    assert thresh_hot is not None
    if outdir is not None:
        ofn = os.path.join(outdir, str(ds)) + ".hdf5"
        if not os.path.exists(outdir):
            try:
                os.makedirs(outdir)
            except: # necessary for racetime
                pass
        if os.path.isfile(ofn):
            logging.info("%s exists already!", ofn)
            return
    logging.info("Analyze clumps for %s.", str(ds))

    g = ds.covering_grid(level=0, left_edge=ds.domain_left_edge,
                         dims=ds.domain_dimensions)
    if thresh_unit is None:
        dat = g[mode].value
    else:
        dat = g[mode].to(thresh_unit).value

    if mode == 'temp' or mode == 'temperature':
        blobs = dat < threshold
    else:
        blobs = dat > threshold

    la, nf = find_clumps(blobs)

    # Find sizes
    ids, sizes = np.unique(la, return_counts=True)
    logging.info("Found %d clumps above minimum size of %d cells.", np.sum(sizes > minsize) - 1, minsize)
    ids   = ids[sizes >= minsize][1:]
    sizes = sizes[sizes >= minsize][1:]
    if len(ids) == 0:
        logging.error("No clumps found!")
        return

    # Compute clump positions in code units
    logging.info("Computing clump positions.")
    comi = np.array(scipy.ndimage.measurements.center_of_mass(blobs,la,ids))

    # For particles
    if clumps_to_particles:
        logging.info("Computing clumps for particles...")
        df_particles = _compute_clumps_to_particles(ds, la)
        ofn_particles = os.path.join(outdir, str(ds)) + "_particles.hdf5"

    logging.info("Analyze clumps...")
    clumpdat = []
    #for i,cid in tqdm(enumerate(ids), total = len(ids)):
    for i,cid in enumerate(ids):
        clumpdat.append(_analyze_clump(cid, la, g, comi[i], sizes[i], thresh_hot, mode, thresh_unit = thresh_unit))
    g.clear_data() # Try to free some memory

    df = pandas.DataFrame(clumpdat).set_index('cid')

    # Convert some quantities to code units
    for i,k in enumerate('xyz'):
        df[k] = ((df['i' + k] / np.array(blobs.shape)[i]) * ds.domain_width[i] + ds.domain_left_edge[i])
    Vcell = np.product(ds.domain_width / ds.domain_dimensions)
    df['size'] = df['isize'] * Vcell

    if outdir is not None:
        logging.info("Writing to %s", ofn)
        df.to_hdf(ofn,'dat')
        if clumps_to_particles:
            df_particles.to_hdf(ofn_particles,'dat')
    else:
        if clumps_to_particles:
            return df, df_particles
        else:
            return df


def find_clumps(blobs):
    """
    Keyword Arguments:
    blobs -- Array of True / False entries to find clumps
    """

    la, nf = scipy.ndimage.measurements.label(blobs)
    logging.info("Found %d clumps.",  nf)

    return la, nf


def extra_info(df):
    """Fills in additional columns.
    """
    df['r_d'] = (df['size']/(4/3. * np.pi))**(1/3.)
    v_shear = np.sqrt((df.velocity_x_clump - df.velocity_x_hot)**2 + \
                      (df.velocity_y_clump - df.velocity_y_hot)**2 + \
                      (df.velocity_z_clump - df.velocity_z_hot)**2)
    df['v_shear'] = v_shear
    df['chi'] = df['density_clump'] / df['density_hot']
    df['P_clump'] = df['density_clump'] * df['temp_clump']
    df['P_hot'] = df['density_hot'] * df['temp_hot']
    df['t_cc'] = df['chi']**0.5 * df['r_d'] / df['v_shear']
    df['mass'] = df['size'] * df['density_clump']


if __name__ == '__main__':
    cfile = '/Users/maxbg/uni_aktuell/postdoc/hydro/athena++/tests/coag_fog/03_dropshape/out/cloud.out2.00010.athdf'
    #ds = test_clumping_factor(large_file)
    import yt
    import helpers
    helpers.logging_setup()

    ds = yt.load(cfile)

    Tfloor = 2.9094759921078023e-08
    df = clump_analysis(ds, 2 * Tfloor, 10 * Tfloor)


