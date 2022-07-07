module helper
    
    implicit none
    real(8), parameter :: pi = 4.D0*DATAN(1.D0)
    
contains
    
    subroutine make_ymat(sig, N, J, ymat)
!       Makes a [N-J, J] Hankel matrix ymat from the vector y

        !   Input
        integer(8), intent(in) :: N, J
        real(8), intent(in)    :: sig(N)
        !   Output
        real(8), intent(out)   :: ymat(N-J, J)
        !   Internal
        integer(8) :: i, nlj
        
        nlj = N - J

        do i=1, J
            ymat(:,i) = sig(i:(J+i))
        end do

    end subroutine make_ymat
    

    subroutine make_rmat(J, flim, alim, dt, len_f, len_a, rmat)
!       Makes the matrix of right vectors needed for dmusic analysis
        
        !   Input
        integer(8), intent(in)  :: J, len_f, len_a
        real(8),    intent(in)  :: alim(2), flim(2), dt
        !   Output
        complex(8), intent(out) :: rmat(J,len_f*len_a)
        !   Internal
        integer(8)              :: i, h1, h2
        complex(8)              :: tmp(J), sk
        real(8)                 :: w_arr(len_f), a_arr(len_a)
        
        w_arr = 2*pi*real_linspace(flim(1), flim(2), len_f)
        a_arr = real_linspace(alim(1), alim(2), len_a)

        do i=1, len_f*len_a
            h1 = mod((i-1), len_a)+1
            h2 = (i-1)/(len_a)+1
            sk = complex(a_arr(h1), w_arr(h2))
            call make_rvec(sk, J, dt, tmp)
            rmat(:, i) = tmp
        end do
        
    end subroutine make_rmat
    

    subroutine make_rvec(sk, J, dt, rvec)
!		Makes a normalized vector v_i = alpha*exp(i*sk)
!		where alpha is the normalization constant
!			alpha = sqrt(sum(abs(v_i)**2))

        !   Input
        integer(8), intent(in) :: J
        real(8),    intent(in) :: dt
        complex(8), intent(in) :: sk
        !   Output
        complex(8), intent(out) :: rvec(J)
        !   Internal
        complex(8) :: tvec(J)
        integer(8) :: i

        do i=1, J
            tvec(i) = exp((i-1)*dt*sk)
        end do
        
        call complex_normalize(tvec, size(tvec), rvec)
        
    end subroutine make_rvec
    

    subroutine complex_normalize(vec, N, res)
!		Given a complex vector as input, returns a normalized (also complex)
!		vector

        !   Input
        complex(8), intent(in)  :: vec(N)
        integer, intent(in)     :: N
        !   Output
        complex(8), intent(out) :: res(N)
        !   Internal
        real(8)                 :: sumv
        
        sumv = 1/sqrt(sum(abs(vec)**2))
        res = vec * sumv

    end subroutine complex_normalize
    

    function real_linspace(xi, xf, np) result(out)
!       Makes a linearly spaced vector between xi, xf
!       with np elements

        !   Input
        real(8)    :: xi, xf
        integer(8) :: np
        !   Output
        real(8)    :: out(np)
        !   Internal
        integer(8) :: i
        
        do i=1, np
            out(i) = xi + (xf-xi)/(np-1)*(i-1)
        end do

    end function real_linspace


    function trapezoid(x, y) result(out)
!       Calculates trapezoid integral for uniformly spaced vector x
!       and function y(x)

        !   Input
        real(8)    :: x(:), y(:)
        !   Output
        real(8)    :: out
        !   Internal
        real(8)    :: h
        integer(8) :: i

        h = (x(size(y)) - x(1))/(size(y)-1)
        
        out = (y(1) + y(size(y)))/2
        do i = 2, size(y)-1
            out = out + y(i)
        end do

        out = out*h
    
    end function trapezoid
    

    subroutine svdnov(M, d1, d2, VT)
!		Comentario: Da diferentes signos para los vectores a numpy.
!		No tengo muy claro cuánto importa esto
        
        integer(8), intent(in) :: d1, d2
        real(8),    intent(in) :: M(d1, d2)
        
        
        real(8), intent(out) :: VT(d2, d2)
        
        integer :: info
        real(8) :: U(d1, d1), work(5*min(d1,d2)), S(min(d1, d2))

        call DGESVD('N', 'A', d1, d2, M, d1, S, U, d1, VT, d2, work, 5*min(d1,d2), info)

    end subroutine svdnov


end module helper