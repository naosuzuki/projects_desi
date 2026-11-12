-- igmdust: HSC DR4 (S23B) Wide background galaxy sample for the
-- galaxy-dust reddening measurement (crayon-style backgrounds).
-- Depth per Nao 2026-11-06 (second revision): i < 26 — the full Wide
-- depth; the analysis depth cut, the quiescent/red selection, and the
-- 25<i<26 blending-diagnostic split are all applied LOCALLY afterwards.
-- photo-z 0.4-1.7. Expected: ~4-6e8 rows, ~30-50 GB FITS.
-- Run via the CAS web UI (hscdata.mtk.nao.ac.jp/datasearch) or
-- 00_submit_hsc_queries.sh.
SELECT f.object_id, f.ra, f.dec,
       f.g_cmodel_mag,    f.r_cmodel_mag,    f.i_cmodel_mag,
       f.z_cmodel_mag,    f.y_cmodel_mag,
       f.g_cmodel_magerr, f.r_cmodel_magerr, f.i_cmodel_magerr,
       f.z_cmodel_magerr, f.y_cmodel_magerr,
       f.a_g, f.a_r, f.a_i, f.a_z, f.a_y,
       f.i_extendedness_value,
       p.photoz_best, p.photoz_risk_best, p.stellar_mass, p.sfr
FROM   s23b_wide.forced f
       JOIN s23b_wide.photoz_mizuki p USING (object_id)
WHERE  f.isprimary
AND    f.i_cmodel_mag < 26.0
AND    f.i_extendedness_value = 1
AND    NOT f.i_pixelflags_saturatedcenter
AND    NOT f.i_pixelflags_edge
AND    NOT f.i_pixelflags_bad
AND    p.photoz_best BETWEEN 0.4 AND 1.7
;
