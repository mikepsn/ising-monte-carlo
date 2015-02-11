c------------------------------------------------------------------------------

      program main
            implicit none
            integer seed, Nx, Ny, i, j, k
		integer start, generate_spin, energy
            parameter(Nx=16)
            parameter(Ny=16)
            integer spins(Nx,Ny)
		integer ene, Ns
		real*8 beta

		open(unit=1, file='16x16.dat', status='unknown')

		write(*,*)'Enter Seed for Random Number Generator'
		read(*,*), seed
		call srand(seed)

            write(*,*)'Enter 1 for a Cold Start and 2 for a Hot Start'
            read(*,*) start

            do i = 1, Nx, +1 
                  do j = 1, Ny, +1 
                        spins(i,j) = generate_spin(start) 
                  end do
            end do

c            call print_lattice(spins, Nx, Ny)

		ene = energy(spins, Nx, Ny)
		write(*,*)'Initial Energy = ',ene

		Ns = 1000
		beta = 0.1d0
		do k = 1, 90, +1
			beta = beta + 0.01d0
			call metropolis(spins, beta, Ns, Nx, Ny)
		end do

		close(unit=1)
      end

c------------------------------------------------------------------------------

      integer function generate_spin(start)
            implicit none
            integer start, s
            real*8 rval, rand

            if (start .eq. 1) then
                  s = 1
            else
                  rval = rand()
                  if (rval .lt. 0.5d0) then
                        s = -1
                  else
                        s = 1
                  end if
            end if

            generate_spin = s
            return
      end

c------------------------------------------------------------------------------

      subroutine print_lattice(spins, Nx, Ny)
            implicit none
            integer i, j, Nx, Ny
            integer spins(Nx,Ny)

            write(*,'(a,$)')'50[;H[2J'
            write(*,'(a,$)')'[H'

            do i = 1, Nx, +1
                  do j = 1, Ny, +1
                        if (spins(i,j) .eq. 1) then
                              write(*,100) '+'
                        end if
                        if (spins(i,j) .eq. -1) then
                              write(*,100) '-'
                        end if
100                     format(a2, $)
                  end do
                  write(*,*)
            end do
      end

c------------------------------------------------------------------------------

	integer function energy(spins,Nx,Ny)
		implicit none
		integer i, j, Nx, Ny
		integer spins(Nx,Ny)
		integer left, right, up, down
		integer en, total_energy

		en = 0
		total_energy = 0
		left = 1
		right = 1
		up = 1
		down = 1

		do i = 1, Nx, +1
			do j = 1, Ny, +1

				left = i - 1
				right = i + 1
				up = j - 1 
				down = j + 1

				if (i .eq. 1) then 
					left = Nx
				end if

				if (i .eq. Nx) then
					right = 1
				end if

				if (j .eq. 1) then
					up = Ny
				end if

				if (j .eq. Ny) then
					down = 1
				end if

				en = -spins(i,j)*(spins(left,j)+
     > 			spins(i,up)+spins(right,j)+spins(i,down))
				total_energy = total_energy + en
			end do
		end do

		energy = total_energy
		return
	end

c------------------------------------------------------------------------------

	subroutine metropolis(spins, beta, Ns, Nx, Ny)
		implicit none	
		real*8 M, beta, eta, rand, error
		integer i_t, j_t, energy, currentEnergy, newEnergy
		integer magnetisation
		integer deltaEnergy
            integer sweep, i, j, Nx, Ny, Ns
            integer spins(Nx,Ny)
		integer thermalise, N, mag, msquared

		thermalise = 100
		N = 0 
		M = 0
		
		do sweep = 1, (Ns + thermalise), +1
			currentEnergy = energy(spins, Nx, Ny)

			i_t = idint(Nx*rand()) + 1
			j_t = idint(Ny*rand()) + 1
			spins(i_t, j_t) = -spins(i_t, j_t)

			newEnergy = energy(spins, Nx, Ny)
			deltaEnergy = newEnergy - currentEnergy
			eta = rand()

			if (deltaEnergy .ge. 0) then
				if (dexp(-beta*deltaEnergy) .le. eta) then
					spins(i_t, j_t) = -spins(i_t, j_t)
				end if
			end if

c			call print_lattice(spins, Nx, Ny)
c			pause

			if (sweep .gt. thermalise) then
				if (mod(sweep, 10) .eq. 0d0) then
					N = N + 1
					mag = magnetisation(spins, Nx, Ny)
					M = M + mag
					msquared = msquared + mag**2
				end if
			end if
		end do

		M = M/dfloat(N)
		msquared = msquared/dfloat(N)
		error = dsqrt((msquared - M**2)/dfloat(N))
		write(*,*)'beta, M, error = ', beta, dabs(M), error
		write(1,*) beta, dabs(M), error
	end
c------------------------------------------------------------------------------

	integer function magnetisation(spins, Nx, Ny)
		implicit none
		integer i, j, Nx, Ny
		integer spins(Nx,Ny)
		integer mag

		mag = 0

		do i = 1, Nx, +1
			do j = 1, Ny, +1
				mag = mag + spins(i,j)
			end do
		end do

		magnetisation = mag
		return
	end

