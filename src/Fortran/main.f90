subroutine single_dmusic(x, y, len_t, rlim, nr, flim, nf, &
    N, J, K, psd, fs)
    !
    !		Pedro Pons-Villalonga (2022)
    !		Implementation of the one-dimensional DMUSIC algorithm.
    !
    !		INPUT:
    !			x:		time vector
    !			y:		time-series vector
    !           len_t:  Length of the time vector
    !			rlim: 	Limits on the damping constant vector (2-element array)
    !			nr:		Number of points in the damping constant vector
    !			flim:	Limits on the frequency vector (2-element array)
    !			nf:		Number of points in the frequency vector. Length of $psd
    !			N:		Length of $x, $y
    !			J:		Number of columns in $my
    !			K:		Number of ignored vectors in V
    !
    !       OUTPUT:
    !           psd:    Power spectrum
    !           freqs:  Frequency array
    
    use helper

    implicit none

    !   Input
    real(8),    intent(in)  :: x(len_t), y(len_t), rlim(2), flim(2)
    integer(8), intent(in)  :: nr, nf, J, K, N, len_t
    !   Output
    real(8),    intent(out)   :: psd(nf), fs(nf)
    !   Internal
    real(8)    :: DT, rs(nr), ws(nf), ym(N-J, J), P_real(nf), ymat(N-J, J)
    real(8)    :: fnyq, VT(J, J), vn(J-K, J), ri, ji, rnorm
    complex(8) :: P_comp(nf), tmps, rmat(J, nr*nf), tmpr
    integer(8) :: i, ik, jf, ir
    
    DT = (x(len_t) - x(1))/(len_t - 1) ! Sampling period
    fnyq = 1/(2*DT)
    
    rs    = real_linspace(rlim(1), rlim(2), nr)
    fs    = real_linspace(flim(1), flim(2), nf)

    ! TODO: Cambiar el paso flim, rlim a fs, rs
    call make_rmat(J, flim, rlim, DT, nf, nr, rmat) 
    ! END TODO:
    call make_ymat(y, N, J, ymat) ! Calculates Hankel matrix
    print *, "ymat shape\t", shape(ymat)
    call svdnov(ym, N-J, J, VT) ! SVD to get V^T
    vn = VT((K+1):(N-J), :) ! Discard the first K columns of V^T
    ! tmp = matmul(vn, rmat) 
    print *, shape(vn)

end subroutine single_dmusic



!end module dmusic
