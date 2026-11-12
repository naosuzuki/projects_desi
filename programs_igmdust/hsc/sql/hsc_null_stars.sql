-- igmdust: HSC DR4 (S23B) Wide star sample for the zero-reddening null
-- test (point sources behind the same foreground galaxies must show no
-- color excess). Bright cut keeps extendedness reliable.
-- Expected: ~1e7 rows, ~1 GB.
SELECT f.object_id, f.ra, f.dec,
       f.g_cmodel_mag,    f.r_cmodel_mag,    f.i_cmodel_mag,
       f.z_cmodel_mag,    f.y_cmodel_mag,
       f.g_cmodel_magerr, f.r_cmodel_magerr, f.i_cmodel_magerr,
       f.z_cmodel_magerr, f.y_cmodel_magerr,
       f.a_g, f.a_r, f.a_i, f.a_z, f.a_y
FROM   s23b_wide.forced f
WHERE  f.isprimary
AND    f.i_cmodel_mag < 22.5
AND    f.i_extendedness_value = 0
AND    NOT f.i_pixelflags_saturatedcenter
AND    NOT f.i_pixelflags_edge
AND    NOT f.i_pixelflags_bad
;
