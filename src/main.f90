Program test2d
    Use Variabili
    Implicit none

	!call system("rm *.plt") !ATTENZIONE: cancella tutti i file plt presenti nella cartella! (su linux)
	call system("del *.plt") !ATTENZIONE: cancella tutti i file plt presenti nella cartella! (su windows)

    call input

    call leggi_gmsh

    call processa_elementi

    !Verifica delle caratteristiche geometria per elemento 10 e interfaccia 10
    write(*,*)'ele(10)%x0 = ',ele(10)%x0
    write(*,*)'ele(10)%area = ',ele(10)%area
    write(*,*)'interf(100)%length = ',interf(100)%length
    write(*,*)'interf(100)%normal = ',interf(100)%normal

	call init

	call write_tecplot_file(0,0.) !primo parametro determina il nome del file, secondo parametro determina la traslazione lungo y


    time=0.

	do k=1,kfinal

		call compute_dt
		call compute_fluxes
		call integ_time

		time=time+dt

		if(mod(k,kinf).eq.0)then
			write(*,*)' '
			write(*,*)'*************************************************'
			write(*,*)'k,time,dt = ',k,time,dt
			call compute_norm_residuals
			call compute_norm_entropy
			call save_wall_data

		end if


		if(mod(k,kout).eq.0)then
			call write_tecplot_file(k,0.)
		end if


	end do



End Program






