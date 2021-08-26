Timer unit: 1e-06 s

Total time: 4.41261 s
File: /home/ppxlf2/sources/candelizer/background.py
Function: search_input_position at line 160

Line #      Hits         Time  Per Hit   % Time  Line Contents
==============================================================
   160                                           @profile
   161                                           def search_input_position(cutout):
   162                                           
   163        20         48.0      2.4      0.0      sigma = 3.0 * gaussian_fwhm_to_sigma 
   164        20      69624.0   3481.2      1.6      kernel = Gaussian2DKernel(sigma, x_size=3, y_size=3)
   165        20        736.0     36.8      0.0      kernel.normalize()
   166        20    1743683.0  87184.1     39.5      threshold = detect_threshold(cutout.data, nsigma=2.)
   167        20     295689.0  14784.5      6.7      segm = detect_sources(cutout.data, threshold, npixels=5, filter_kernel=kernel)
   168                                           
   169        20        143.0      7.2      0.0      segm_deblend = deblend_sources(cutout.data, segm, npixels=5,
   170        20         28.0      1.4      0.0                                  filter_kernel=kernel, nlevels=5,
   171        20    2068681.0 103434.1     46.9                                  contrast=0.001)
   172                                           
   173        20      61198.0   3059.9      1.4      idx = np.random.choice(np.arange(len((np.where(segm_deblend.data==0)[0]))), 1)
   174        20      51083.0   2554.2      1.2      x = np.where(segm_deblend.data==0)[0][idx][0]
   175        20      52515.0   2625.8      1.2      y = np.where(segm_deblend.data==0)[0][idx][0]
   176                                           
   177        20      69182.0   3459.1      1.6      return cutout.wcs.pixel_to_world(x, y)

Total time: 5.2647 s
File: /home/ppxlf2/sources/candelizer/background.py
Function: find_CANDELS_sky_patch at line 179

Line #      Hits         Time  Per Hit   % Time  Line Contents
==============================================================
   179                                           @profile
   180                                           def find_CANDELS_sky_patch(bands, field):
   181                                           
   182        20         26.0      1.3      0.0      looking_for_input = True
   183        21         21.0      1.0      0.0      while looking_for_input:
   184                                           
   185        21         28.0      1.3      0.0          skypatches = []
   186        21         20.0      1.0      0.0          input_coord = None
   187        21         19.0      0.9      0.0          try:
   188        21       1130.0     53.8      0.0              print(field)
   189        21         48.0      2.3      0.0              sources = sources_dfs[field]
   190        21      17703.0    843.0      0.3              source = sources.iloc[np.random.randint(0, sources.shape[0])]
   191        21      36002.0   1714.4      0.7              coord = SkyCoord(ra=source.RA * u.deg, dec=source.DEC* u.deg)
   192                                                       
   193        81        255.0      3.1      0.0              for band, filter in zip(bands, FIELDS[field]):
   194                                           
   195                                                           #gal_img = fits.getdata(f'/data/captain/MERGERS/SKIRT/TNG50-1/PMxSF/TESTS/{filter}/30_100453_oct_1.fits')
   196                                           
   197        61        254.0      4.2      0.0                  wcs  = WCS_MEMORY[field][filter]#header = fits.getheader(f'/home/ppxlf2/CANDELS/{field}/{FIELDS[field][filter]["fits_filename"]}.fits')
   198        61        128.0      2.1      0.0                  data = FIELDS[field][filter]['data']#fits.getdata(f'/home/ppxlf2/CANDELS/{field}/{FIELDS[field][filter]["fits_filename"]}.fits', memmap=True)
   199                                                           
   200                                                           #wcs = WCS(header)
   201                                                           #wcs.sip = None
   202                                           
   203        61         62.0      1.0      0.0                  if input_coord is None:
   204                                           
   205        21     204679.0   9746.6      3.9                      cutout = Cutout2D(data, position=coord, size=512, wcs=wcs)
   206                                           
   207        21      41872.0   1993.9      0.8                      if np.all(cutout.data == 0):
   208         1         10.0     10.0      0.0                          raise ValueError('Empty Patch of the Sky')
   209                                           
   210        20    4414265.0 220713.2     83.8                      input_coord = search_input_position(cutout)
   211                                           
   212                                           
   213        60     547848.0   9130.8     10.4                  cutout = Cutout2D(data, position=input_coord, size=256, wcs=wcs)
   214                                           
   215        60        215.0      3.6      0.0                  skypatches.append(cutout.data)
   216         1          3.0      3.0      0.0          except ValueError as e:
   217         1         68.0     68.0      0.0              print(e)
   218         1          2.0      2.0      0.0              continue
   219                                                   
   220        20         25.0      1.2      0.0          break
   221                                           
   222        20         21.0      1.1      0.0      return skypatches, field, input_coord

Total time: 35.6041 s
File: candelizer.py
Function: run at line 27

Line #      Hits         Time  Per Hit   % Time  Line Contents
==============================================================
    27                                           @profile
    28                                           def run():
    29                                               """
    30                                                   TNG Data is stored as {snapshot}_{subfindID}. Example: 99_234510
    31                                               """
    32         1          3.0      3.0      0.0      rootname = '60_303520'#sys.argv[1]
    33                                           
    34         1          7.0      7.0      0.0      seed = int(''.join(rootname.split('_')))
    35                                               
    36                                               
    37                                               """
    38                                                   Create a ParameterGrid to avoid 3 nested loops
    39                                                   The runs are independent on FIELD/REDSHIFT/ORIENTATION combination.
    40                                               """
    41         1          3.0      3.0      0.0      parameter_space = ParameterGrid({
    42         1          5.0      5.0      0.0                          'FIELDS': list(FIELDS.keys()),
    43         1          4.0      4.0      0.0                          'orientations' : cfg.ALLOWED_ORIENTATIONS,
    44         1         20.0     20.0      0.0                          'redshifts' : cfg.TARGET_REDSHIFTS
    45                                                                 })
    46                                           
    47        21        340.0     16.2      0.0      for i, params in enumerate(parameter_space):
    48                                           
    49        20         88.0      4.4      0.0          seed += i*787
    50        20        476.0     23.8      0.0          np.random.seed(seed)
    51                                                   
    52        20         68.0      3.4      0.0          field = params['FIELDS']
    53        20         57.0      2.9      0.0          orientation = params['orientations']
    54        20         60.0      3.0      0.0          target_redshift = params['redshifts']
    55                                                   
    56        20        812.0     40.6      0.0          primary_header = fits.Header()
    57                                               
    58        20         83.0      4.2      0.0          filename = f'{rootname}_{orientation}'
    59                                           
    60        20       5161.0    258.1      0.0          primary_header['ROOTNAME'] = rootname
    61        20       4259.0    212.9      0.0          primary_header['FILENAME'] = filename
    62        20       4179.0    208.9      0.0          primary_header['NP_SEED'] = seed
    63        20       4158.0    207.9      0.0          primary_header['SKIRTORI'] = orientation
    64        20       4238.0    211.9      0.0          primary_header['CANDELSF'] = field
    65        20       4342.0    217.1      0.0          primary_header['TARGET_Z'] = target_redshift
    66        20       4365.0    218.2      0.0          primary_header['BANDS'] = ','.join(cfg.BANDS)
    67                                           
    68        20        336.0     16.8      0.0          output_headers = {}
    69        80        288.0      3.6      0.0          for band in cfg.BANDS:
    70        60        750.0     12.5      0.0              hdr = fits.Header()
    71        60      12160.0    202.7      0.0              hdr['BAND'] = band
    72        60        233.0      3.9      0.0              output_headers[band] = hdr
    73                                           
    74        20    2914553.0 145727.6      8.2          cube = SkirtCube(filename, cfg.SOURCE_FOLDER)
    75                                           
    76        20        180.0      9.0      0.0          if target_redshift != 'full':
    77                                                       cube.target_z = target_redshift
    78                                           
    79        20      21928.0   1096.4      0.1          dimming = (1/(1+cube.target_z))*((cosmo.luminosity_distance(SKIRT_Z0) / cosmo.luminosity_distance(cube.target_z))**2).value # Dimming calculated from the ratio of the luminosities of target and source frames
    80        20     155220.0   7761.0      0.4          cube.data *= dimming 
    81                                           
    82        20       8936.0    446.8      0.0          primary_header['DIMMFACT'] = dimming
    83                                           
    84        20    9670350.0 483517.5     27.2          broadbands = [cube.integrate_filter(band)  for band in cfg.FILTERS]  
    85        20     505084.0  25254.2      1.4          broadbands = [to_eps(band, header).value    for band, header in zip(broadbands, cfg.HEADERS)]  
    86                                           
    87        20        454.0     22.7      0.0          if (cfg.fix_scaling) & (cube.z < REFERENCE_Z):
    88                                                       correct_scale = COSMO.scale_factor(0.5)/COSMO.scale_factor(cube.z)
    89                                           
    90                                                       for i, band in enumerate(broadbands):
    91                                                           flux = band.sum()
    92                                                           band = zoom(band, correct_scale)
    93                                                           band = (band/band.sum())*flux
    94                                                           broadbands[i] = band
    95                                           
    96                                           
    97        20         78.0      3.9      0.0          if cfg.bg:
    98        20    5267107.0 263355.3     14.8              bg_sections, _, coord = find_CANDELS_sky_patch(broadbands, field)
    99        20      32696.0   1634.8      0.1              primary_header['RA_BG'] = coord.ra.value
   100        20       8937.0    446.9      0.0              primary_header['DEC_BG'] = coord.dec.value
   101                                                   else:
   102                                                       bg_sections = [None] * len(cfg.FILTERS)
   103                                           
   104        20        277.0     13.8      0.0          restframe = areia.ObservationFrame(REFERENCE_Z, pixelscale=0.05, exptime=1)
   105                                           
   106        20         94.0      4.7      0.0          if cfg.restframe:
   107                                                       obs_frames = [restframe for tps in cfg.TARGET_PS]
   108                                                   else:
   109        20        389.0     19.4      0.0              obs_frames = [areia.ObservationFrame(cube.target_z, pixelscale=tps, exptime=exptime) for tps, exptime in zip(cfg.TARGET_PS, cfg.EXPTIME)]
   110                                           
   111        20        409.0     20.4      0.0          result = {}
   112        80        498.0      6.2      0.0          for band, psf, bg, obs_frame, bandname in zip(broadbands, cfg.PSFS, bg_sections, obs_frames, cfg.BANDS):
   113        60   13820170.0 230336.2     38.8              AR = areia.ArtificialRedshift(band, psf, bg, restframe, obs_frame, MAG=None, bg_position=None, config=cfg.areia_config)
   114        60        395.0      6.6      0.0              result[bandname] = {}
   115        60        220.0      3.7      0.0              result[bandname]['final'] = AR.final_crop
   116        60        192.0      3.2      0.0              result[bandname]['bg'] = AR.background
   117        60        197.0      3.3      0.0              result[bandname]['clean'] = AR.convolved
   118                                           
   119        80        425.0      5.3      0.0          for i, band in enumerate(cfg.BANDS):
   120        60        249.0      4.2      0.0              keyword = f'BG{band}'
   121                                                       
   122                                                       
   123        60     750985.0  12516.4      2.1              segmap, flux_per_pixel, flux, num_pix = segmentate(result[band]['bg'])
   124        60    1291450.0  21524.2      3.6              galsegmap, _, _, _ = segmentate(result[band]['clean'])
   125                                                       
   126        60       5574.0     92.9      0.0              overlap_index = np.where((segmap >= 1) & (galsegmap >= 1))
   127        60        278.0      4.6      0.0              overlapping = overlap_index[0].shape[0]
   128                                           
   129        60       2379.0     39.6      0.0              overlap_map = np.zeros_like(galsegmap)
   130        60        929.0     15.5      0.0              overlap_map[overlap_index] = 1
   131                                           
   132        60       4469.0     74.5      0.0              coverage = np.round(100 * overlapping / len(galsegmap[galsegmap >= 1]), 2)
   133                                           
   134        60      21883.0    364.7      0.1              output_headers[band]['FLUX_PP'] = flux_per_pixel
   135        60      13776.0    229.6      0.0              output_headers[band]['FLUX_INT'] = flux
   136        60      11952.0    199.2      0.0              output_headers[band]['NUM_PIX'] = num_pix
   137        60      11745.0    195.8      0.0              output_headers[band]['OVERLAP'] = overlapping
   138        60      14105.0    235.1      0.0              output_headers[band]['COVERAGE'] = coverage
   139                                           
   140        60        277.0      4.6      0.0              result[band]['segmap'] = segmap
   141        60        193.0      3.2      0.0              result[band]['clean_segmap'] = galsegmap
   142        60        206.0      3.4      0.0              result[band]['overlap_segmap'] = overlap_map
   143                                                   
   144        20         88.0      4.4      0.0          if cfg.shotnoise:
   145        80        354.0      4.4      0.0              for filter in FIELDS[field]:
   146        60      49134.0    818.9      0.1                  result[filter]['final'] = result[filter]['final'] + np.random.normal(0, FIELDS[field][filter]['skystd'], size=result[filter]['clean'].shape)
   147        60      45307.0    755.1      0.1                  result[filter]['clean'] = result[filter]['clean'] + np.random.normal(0, FIELDS[field][filter]['skystd'], size=result[filter]['clean'].shape)
   148                                           
   149        20      51728.0   2586.4      0.1          primary_header = config_to_header(primary_header)
   150                                           
   151        20     323278.0  16163.9      0.9          hdul = dict_to_hdu(fits.PrimaryHDU(header=primary_header), result, output_headers)
   152                                           
   153        20        160.0      8.0      0.0          output_path = f'{cfg.OUTPUT_FOLDER}/{rootname}_{orientation}_{target_redshift}_{field}.fits'
   154                                           
   155        20     548275.0  27413.8      1.5          hdul.writeto(output_path, overwrite=True)

