import sncosmo


#Change band names to ones used in P60 dataframe
band_g = sncosmo.get_bandpass('sdssg')
band_g.name = 'g'
sncosmo.register(band_g)

band_r = sncosmo.get_bandpass('sdssr')
band_r.name = 'r'
sncosmo.register(band_r)

band_i = sncosmo.get_bandpass('sdssi')
band_i.name = 'i'
sncosmo.register(band_i)

band_u = sncosmo.get_bandpass('sdssu')
band_u.name = 'u'
sncosmo.register(band_u)

