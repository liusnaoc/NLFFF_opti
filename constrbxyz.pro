;;-
;Pro Optimization_fff, bx_in, by_in, bz_in, $
                       
 ;                     Bc_flag = bc_flag, $
cd,'F:\document-1\nlff_fortianh\nlfff2\pythonlfff\'
compile_opt idl2
restore,'F:\document-1\nlff_fortianh\nlfff2\test_data\CaseII_theory_B3_64x64_-1_1.sav'

;tb3x=congrid(tb3x,64,32,64)
;tb3y=congrid(tb3y,64,32,64)
;tb3z=congrid(tb3z,64,32,64)
bx=tb3x[*,*,0]
by=tb3y[*,*,0]
bz=tb3z[*,*,0]

writefits,'bx.fits',bx
writefits,'by.fits',by
writefits,'bz.fits',bz
end