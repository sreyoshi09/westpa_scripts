#!/bin/tcsh


foreach x ( 1 2 3 4 5 6 7 8 9 10 11 )
	echo "resid == $x"
	echo res_$x.h5
	python res_cont_h5.py /mnt/synuclein/8mer_run1/alpha.gro /mnt/synuclein/8mer_run1/bstates/model.psf "resid == $x" 'resid >= 12 && resid <=88' 8 /mnt/synuclein/8mer_run1 res_$x.h5
end
