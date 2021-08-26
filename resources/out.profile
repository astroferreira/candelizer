Timer unit: 1e-06 s

Total time: 4.88535 s
File: /home/ppxlf2/sources/candelizer/background.py
Function: search_input_position at line 160

Line #      Hits         Time  Per Hit   % Time  Line Contents
==============================================================
   160                                           @profile
   161                                           def search_input_position(cutout):
   162                                           
   163        20         58.0      2.9      0.0      sigma = 3.0 * gaussian_fwhm_to_sigma 
   164        20      94514.0   4725.7      1.9      kernel = Gaussian2DKernel(sigma, x_size=3, y_size=3)
   165        20        880.0     44.0      0.0      kernel.normalize()
   166        20    1868905.0  93445.2     38.3      threshold = detect_threshold(cutout.data, nsigma=2.)
   167        20     255815.0  12790.8      5.2      segm = detect_sources(cutout.data, threshold, npixels=5, filter_kernel=kernel)
   168                                           
   169        20        127.0      6.3      0.0      segm_deblend = deblend_sources(cutout.data, segm, npixels=5,
   170        20         18.0      0.9      0.0                                  filter_kernel=kernel, nlevels=5,
   171        20    2382572.0 119128.6     48.8                                  contrast=0.001)
   172                                           
   173        20      80268.0   4013.4      1.6      idx = np.random.choice(np.arange(len((np.where(segm_deblend.data==0)[0]))), 1)
   174        20      64733.0   3236.7      1.3      x = np.where(segm_deblend.data==0)[0][idx][0]
   175        20      56363.0   2818.2      1.2      y = np.where(segm_deblend.data==0)[0][idx][0]
   176                                           
   177        20      81098.0   4054.9      1.7      return cutout.wcs.pixel_to_world(x, y)

Total time: 15.1481 s
File: /home/ppxlf2/sources/candelizer/background.py
Function: find_CANDELS_sky_patch at line 179

Line #      Hits         Time  Per Hit   % Time  Line Contents
==============================================================
   179                                           @profile
   180                                           def find_CANDELS_sky_patch(bands, field):
   181                                           
   182        20         38.0      1.9      0.0      looking_for_input = True
   183        21         21.0      1.0      0.0      while looking_for_input:
   184                                           
   185        21         22.0      1.0      0.0          skypatches = []
   186        21         24.0      1.1      0.0          input_coord = None
   187                                           
   188        21         22.0      1.0      0.0          try:
   189        21       1850.0     88.1      0.0              print(field)
   190        21         53.0      2.5      0.0              sources = sources_dfs[field]
   191        21      20401.0    971.5      0.1              source = sources.iloc[np.random.randint(0, sources.shape[0])]
   192        21      40095.0   1909.3      0.3              coord = SkyCoord(ra=source.RA * u.deg, dec=source.DEC* u.deg)
   193                                                       
   194        81        347.0      4.3      0.0              for band, filter in zip(bands, FIELDS[field]):
   195                                           
   196                                                           #gal_img = fits.getdata(f'/data/captain/MERGERS/SKIRT/TNG50-1/PMxSF/TESTS/{filter}/30_100453_oct_1.fits')
   197                                           
   198        61        414.0      6.8      0.0                  wcs = WCS_MEMORY[field][filter]#header = fits.getheader(f'/home/ppxlf2/CANDELS/{field}/{FIELDS[field][filter]["fits_filename"]}.fits')
   199        61        164.0      2.7      0.0                  data     = FIELDS[field][filter]['data']#fits.getdata(f'/home/ppxlf2/CANDELS/{field}/{FIELDS[field][filter]["fits_filename"]}.fits', memmap=True)
   200                                                           
   201                                                           #wcs = WCS(header)
   202                                                           #wcs.sip = None
   203                                           
   204        61         75.0      1.2      0.0                  if input_coord is None:
   205                                           
   206        21     211969.0  10093.8      1.4                      cutout = Cutout2D(data, position=coord, size=512, wcs=wcs)
   207                                           
   208        21    9300649.0 442888.0     61.4                      if np.all(cutout.data == 0):
   209         1          9.0      9.0      0.0                          raise ValueError('Empty Patch of the Sky')
   210                                           
   211        20    4887260.0 244363.0     32.3                      input_coord = search_input_position(cutout)
   212                                           
   213                                           
   214        60     684325.0  11405.4      4.5                  cutout = Cutout2D(data, position=input_coord, size=256, wcs=wcs)
   215                                           
   216        60        238.0      4.0      0.0                  skypatches.append(cutout.data)
   217         1          2.0      2.0      0.0          except ValueError as e:
   218         1         56.0     56.0      0.0              print(e)
   219         1          2.0      2.0      0.0              continue
   220                                                   
   221        20         24.0      1.2      0.0          break
   222                                           
   223        20         25.0      1.2      0.0      return skypatches, field

Total time: 52.1698 s
File: candelizer.py
Function: run at line 23

Line #      Hits         Time  Per Hit   % Time  Line Contents
==============================================================
    23                                           @profile
    24                                           def run():
    25                                           
    26         1          4.0      4.0      0.0      if __name__ == '__main__':
    27                                           
    28                                                   """
    29                                                       TNG Data is stored as {snapshot}_{subfindID}. Example: 99_234510
    30                                                   """
    31         1          2.0      2.0      0.0          rootname = '60_303520'#sys.argv[1]
    32         1          3.0      3.0      0.0          output_headers = []
    33         2         13.0      6.5      0.0          for target_redshift in cfg.TARGET_REDSHIFTS:
    34         5         17.0      3.4      0.0              for orientation in cfg.ALLOWED_ORIENTATIONS:
    35        24        150.0      6.2      0.0                  for field in list(FIELDS.keys()):
    36        20         79.0      4.0      0.0                      filename = f'{rootname}_{orientation}'
    37                                           
    38        80        268.0      3.4      0.0                      for band in cfg.BANDS:
    39        60        769.0     12.8      0.0                          hdr = fits.Header()
    40        60        205.0      3.4      0.0                          output_headers.append(hdr)
    41                                           
    42                                           
    43        20         63.0      3.1      0.0                      base_folder = cfg.SOURCE_FOLDER
    44        20        838.0     41.9      0.0                      if not os.path.exists(f'{base_folder}/DATACUBE/{filename}_total.fits'):
    45                                                                   base_folder = cfg.SOURCE_FOLDER_FB
    46                                           
    47        20    3231943.0 161597.1      6.2                      cube = SkirtCube(filename, base_folder)
    48        20        238.0     11.9      0.0                      if target_redshift != 'full':
    49                                                                   cube.target_z = target_redshift
    50                                                           
    51                                           
    52        20      29546.0   1477.3      0.1                      dimming = (1/(1+cube.target_z))*((cosmo.luminosity_distance(SKIRT_Z0) / cosmo.luminosity_distance(cube.target_z))**2).value # Dimming calculated from the ratio of the luminosities of target and source frames
    53        20     180228.0   9011.4      0.3                      cube.data *= dimming 
    54        20   11334939.0 566746.9     21.7                      broadbands = [cube.integrate_filter(band)  for band in cfg.FILTERS]  
    55        20     575388.0  28769.4      1.1                      broadbands = [to_eps(band, header).value    for band, header in zip(broadbands, cfg.HEADERS)]  
    56                                           
    57        20        563.0     28.1      0.0                      if (cfg.fix_scaling) & (cube.z < REFERENCE_Z):
    58                                                                   correct_scale = COSMO.scale_factor(0.5)/COSMO.scale_factor(cube.z)
    59                                           
    60                                                                   for i, band in enumerate(broadbands):
    61                                                                       flux = band.sum()
    62                                                                       band = zoom(band, correct_scale)
    63                                                                       band = (band/band.sum())*flux
    64                                                                       broadbands[i] = band
    65                                           
    66                                           
    67        20         93.0      4.7      0.0                      if cfg.bg:
    68        20   15150678.0 757533.9     29.0                          bg_sections, _ = find_CANDELS_sky_patch(broadbands, field)
    69        80        554.0      6.9      0.0                          for i, band in enumerate(cfg.BANDS):
    70        60        441.0      7.3      0.0                              keyword = f'BG{band}'
    71        60    5444856.0  90747.6     10.4                              output_headers[i][keyword] = bg_sections[i].sum()
    72                                                               else:
    73                                                                   bg_sections = [None] * len(cfg.FILTERS)
    74                                                       
    75        20        381.0     19.1      0.0                      restframe = areia.ObservationFrame(REFERENCE_Z, pixelscale=0.05, exptime=1)
    76                                           
    77        20         90.0      4.5      0.0                      if cfg.restframe:
    78                                                                   obs_frames = [restframe for tps in cfg.TARGET_PS]
    79                                                               else:
    80        20        434.0     21.7      0.0                          obs_frames = [areia.ObservationFrame(cube.target_z, pixelscale=tps, exptime=exptime) for tps, exptime in zip(cfg.TARGET_PS, cfg.EXPTIME)]
    81                                           
    82        20   15179036.0 758951.8     29.1                      final = [areia.ArtificialRedshift(band, psf, bg, restframe, obs_frame, MAG=None, bg_position=None, config=cfg.areia_config).final_crop for band, psf, bg, obs_frame in zip(broadbands, cfg.PSFS, bg_sections, obs_frames)]
    83                                           
    84        20        172.0      8.6      0.0                      if cfg.shotnoise:
    85        80        588.0      7.3      0.0                          for i, filter in enumerate(FIELDS[field]):
    86        60      55482.0    924.7      0.1                              final[i] = final[i] + np.random.normal(0, FIELDS[field][filter]['skystd'], size=final[i].shape)
    87                                           
    88        80        776.0      9.7      0.0                      for i, (band, hdr) in enumerate(zip(cfg.BANDS, output_headers)):
    89        60        511.0      8.5      0.0                          output_path = f'{cfg.OUTPUT_FOLDER}/{band}/{rootname}_{orientation}_{target_redshift}_{field}.fits'
    90        60       1663.0     27.7      0.0                          dirname = os.path.dirname(output_path)
    91        60       2825.0     47.1      0.0                          if not os.path.exists(dirname):
    92                                                                       os.makedirs(dirname)
    93                                                               
    94        60     975956.0  16265.9      1.9                          fits.writeto(output_path, final[i], hdr, overwrite=True)

