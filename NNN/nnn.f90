program Ising_NNN
        use randomnumber
        implicit none

        integer :: pMC_total 
        integer, parameter :: pMC_medida = 150    
        integer :: num_medidas 

        real*8, parameter :: J1 = 1.0d0 
        real*8, parameter :: J2 = 2.0d0 

        integer, dimension(6) :: N_vals = (/16, 32, 64, 128, 256, 512/)
        real*8, dimension(25) :: T_vals
        
        integer :: N, in_val, it_val, pmc
        integer :: iseed, i, r
        real*8 :: T
        
        real*8 :: M_actual, E_actual
        real*8 :: sum_M, sum_M2
        real*8 :: sum_E, sum_E2, sum_E3, sum_E4
        real*8 :: avg_E, avg_E2, avg_E3, avg_E4, mu2, mu4
        
        real*8 :: m_N, err_M, e_N, err_E, c_N, err_C
        real*8, allocatable, dimension(:) :: sum_corr, sum_corr2
        real*8 :: avg_corr, avg_corr2
        
        integer*8, allocatable, dimension(:,:) :: S

        iseed = 1565867
        call dran_ini(iseed) !For the NN case it sufficees to set J2=0 and to center the interval of temperatures around 2.269 (Onsager)

        T_vals = (/ 6.00d0, 6.50d0, 7.00d0, 7.40d0, 7.60d0, 7.80d0, 7.90d0, 7.95d0, 8.02d0, 8.05d0, &
                    8.07d0, 8.09d0, 8.10d0, 8.25d0, 8.30d0, 8.35d0, 8.40d0, 8.50d0, 8.60d0, 8.80d0, &
                    9.00d0, 9.50d0, 10.00d0, 10.50d0, 11.00d0 /)       
        
        open(10, file="termodinamica-next-next.dat", status="replace")
        write(10, *) "# N      T         m_N            err_M          e_N            err_E          c_N            err_C"
        
        open(20, file="correlacion-next-next.dat", status="replace")
        write(20, *) "# T      Dist(r)   f(r)           err_f"

        do in_val = 1, 6
            N = N_vals(in_val)
            allocate(S(N,N))
            allocate(sum_corr(N/2))
            allocate(sum_corr2(N/2)) 
            
            if (N >= 256) then
                pMC_total = 250000
                write(*,*) ">> Iniciando red N =", N, "(Modo Rápido: 100k pMC)"
            else
                pMC_total = 1000000
                write(*,*) ">> Iniciando red N =", N, "(Modo Normal: 1M pMC)"
            end if
            num_medidas = pMC_total / pMC_medida

            do it_val = 1, 25
                T = T_vals(it_val)
                write(*, '(A, F5.2, A)', advance='no') "   ->T = ", T, " ... "

                S = 1

                do pmc = 1, 10000 
                    call paso_montecarlo(S, N, T)
                end do
                
                write(*, '(A)', advance='no') " Midiendo ... "

                sum_M = 0.0d0; sum_M2 = 0.0d0
                sum_E = 0.0d0; sum_E2 = 0.0d0; sum_E3 = 0.0d0; sum_E4 = 0.0d0
                sum_corr = 0.0d0; sum_corr2 = 0.0d0

                do pmc = 1, pMC_total
                    call paso_montecarlo(S, N, T)
                    
                    if (mod(pmc, pMC_medida) == 0) then
                        call medir_energia_y_mag(S, N, E_actual, M_actual)
                        
                        sum_M  = sum_M  + abs(M_actual)
                        sum_M2 = sum_M2 + M_actual**2
                        sum_E  = sum_E  + E_actual
                        sum_E2 = sum_E2 + E_actual**2
                        sum_E3 = sum_E3 + E_actual**3
                        sum_E4 = sum_E4 + E_actual**4
                        
                        if (N == 512) then
                            call acumular_correlacion(S, N, sum_corr, sum_corr2)
                        end if
                    end if
                end do

                m_N = (sum_M / dble(num_medidas)) / dble(N*N)
                err_M = sqrt(abs(sum_M2 / dble(num_medidas) - (sum_M / dble(num_medidas))**2) / dble(num_medidas)) / dble(N*N)
                
                avg_E  = sum_E / dble(num_medidas)
                avg_E2 = sum_E2 / dble(num_medidas)
                avg_E3 = sum_E3 / dble(num_medidas)
                avg_E4 = sum_E4 / dble(num_medidas)
                
                e_N = avg_E / dble(N*N) 
                err_E = sqrt(abs(avg_E2 - avg_E**2) / dble(num_medidas)) / dble(N*N)
                
                mu2 = avg_E2 - avg_E**2
                mu4 = avg_E4 - 4.0d0*avg_E3*avg_E + 6.0d0*avg_E2*(avg_E**2) - 3.0d0*(avg_E**4)
                
                c_N = mu2 / (dble(N*N) * T**2)
                err_C = sqrt(abs(mu4 - mu2**2) / dble(num_medidas)) / (dble(N*N) * T**2)

                write(10, '(I5, 7F15.6)') N, T, m_N, err_M, e_N, err_E, c_N, err_C
                
                if (N == 512) then
                    do r = 1, N/2
                        avg_corr = sum_corr(r) / dble(num_medidas)
                        avg_corr2 = sum_corr2(r) / dble(num_medidas)
                        write(20, '(F6.2, I6, 2F15.6)') T, r, avg_corr, sqrt(abs(avg_corr2 - avg_corr**2) / dble(num_medidas))
                    end do
                    write(20, *) "" 
                end if
                
                write(*,*) "Ok."
            end do
            
            write(10, *) "" 
            deallocate(S)
            deallocate(sum_corr)
            deallocate(sum_corr2)
        end do

        close(10)
        close(20)
        write(*,*) " COMPLETADO."

    contains

        subroutine paso_montecarlo(red, N_dim, Temp)
            integer*8, dimension(N_dim,N_dim), intent(inout) :: red
            integer, intent(in) :: N_dim
            real*8, intent(in) :: Temp
            
            integer*8 :: f, c, arriba, izq, dcha, abajo, intento
            real*8 :: prob, deltae, xi
            real*8 :: sum_NN, sum_NNN 

            do intento = 1, N_dim * N_dim
                f = int(1 + dran_u() * N_dim, 8)
                c = int(1 + dran_u() * N_dim, 8)

                arriba = f - 1; if (arriba < 1) arriba = N_dim
                abajo  = f + 1; if (abajo > N_dim) abajo = 1
                izq    = c - 1; if (izq < 1) izq = N_dim
                dcha   = c + 1; if (dcha > N_dim) dcha = 1

                sum_NN = dble(red(arriba, c) + red(abajo, c) + red(f, izq) + red(f, dcha))
                sum_NNN = dble(red(arriba, izq) + red(arriba, dcha) + red(abajo, izq) + red(abajo, dcha))

                deltae = 2.0d0 * dble(red(f, c)) * (J1 * sum_NN + J2 * sum_NNN)

                if (deltae <= 0.0d0) then
                    red(f, c) = -red(f, c)
                else
                    prob = exp(-deltae / Temp)
                    xi = dran_u()
                    if (xi < prob) then
                        red(f, c) = -red(f, c)
                    end if
                end if
            end do
        end subroutine paso_montecarlo

        subroutine medir_energia_y_mag(red, N_dim, E_tot, M_tot)
            integer*8, dimension(N_dim,N_dim), intent(in) :: red
            integer, intent(in) :: N_dim
            real*8, intent(out) :: E_tot, M_tot
            
            integer :: i_sub, j_sub, dcha_sub, abajo_sub, izq_sub
            real*8 :: E_NN, E_NNN
            
            M_tot = 0.0d0
            E_tot = 0.0d0
            
            do i_sub = 1, N_dim
                do j_sub = 1, N_dim
                    M_tot = M_tot + dble(red(i_sub, j_sub))
                    
                    dcha_sub = j_sub + 1; if (dcha_sub > N_dim) dcha_sub = 1
                    abajo_sub = i_sub + 1; if (abajo_sub > N_dim) abajo_sub = 1
                    izq_sub = j_sub - 1; if (izq_sub < 1) izq_sub = N_dim
                    
                    E_NN = dble(red(i_sub, dcha_sub) + red(abajo_sub, j_sub))
                    E_NNN = dble(red(abajo_sub, dcha_sub) + red(abajo_sub, izq_sub))
                    
                    E_tot = E_tot - dble(red(i_sub, j_sub)) * (J1 * E_NN + J2 * E_NNN)
                end do
            end do
        end subroutine medir_energia_y_mag

        subroutine acumular_correlacion(red, N_dim, sum_c, sum_c2)
            integer*8, dimension(N_dim,N_dim), intent(in) :: red
            integer, intent(in) :: N_dim
            real*8, dimension(N_dim/2), intent(inout) :: sum_c, sum_c2
            
            integer :: i_sub, j_sub, r_sub, vecino_x
            real*8 :: corr_actual
            
            do r_sub = 1, N_dim/2
                corr_actual = 0.0d0
                do i_sub = 1, N_dim
                    do j_sub = 1, N_dim
                        vecino_x = j_sub + r_sub
                        if (vecino_x > N_dim) vecino_x = vecino_x - N_dim
                        corr_actual = corr_actual + dble(red(i_sub, j_sub) * red(i_sub, vecino_x))
                    end do
                end do
                sum_c(r_sub)  = sum_c(r_sub)  + (corr_actual / dble(N_dim*N_dim))
                sum_c2(r_sub) = sum_c2(r_sub) + (corr_actual / dble(N_dim*N_dim))**2
            end do
        end subroutine acumular_correlacion

    end program Ising_NNN
