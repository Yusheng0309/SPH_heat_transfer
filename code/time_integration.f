      subroutine time_integration(x,vx, volume, density, p,Tem,u, c, s,
     &           e,  itype, hsml, ntotal, maxtimestep, dt )
     
c----------------------------------------------------------------------
c      x      :coordinates of particles                        [in/out]
c      vx     :velocities of particles                         [in/out]
c      volume :volume of particles                                 [in]
c      density:water density of particles                      [in/out]
c      p      :pressure  of particles                          [in/out]
c      u      :internal energy of particles                    [in/out]
c      Tem    :Temperature                                     [in/out]
c      c      :sound velocity of particles                        [out]
c      s      :entropy of particles, not used here                [out]
c      e      :total energy of particles                          [out]
c      itype  :types of particles                                  [in]
c               =1   ideal gas
c               =2   water
c               =3   tnt
c      hsml   :smoothing lengths of particles                  [in/out]
c      ntotal :total particle number                               [in]  
c      maxtimestep :maximum timesteps                              [in]
c      dt      :timestep                                           [in]
   

      implicit none     
      include 'param.inc.txt'
      
      integer itype(maxn), ntotal, maxtimestep
      double precision x(dim, maxn), vx(dim, maxn), volume(maxn), 
     &       density(maxn), p(maxn), u(maxn), c(maxn), s(maxn), e(maxn), 
     &       hsml(maxn), dt
      integer i, j, k, itimestep, d, current_ts, nstart ,ip       
      double precision  x_min(dim, maxn), v_min(dim, maxn), u_min(maxn),
     &       density_min(maxn), dx(dim,maxn), dvx(dim, maxn), du(maxn),  
     &       ddensity(maxn),  av(dim, maxn), ds(maxn), Tem_min(maxn),
     &       t(maxn), tdsdt(maxn), Tem(maxn), dTem(maxn),storageip(maxn)        
      double precision  time, temp_density, temp_u, temp_Tem
	integer ngrab
	character supp*4,name*40
               
      do i = 1, ntotal
        do d = 1, dim
          av(d, i) = 0.
        enddo
      enddo  
     
      do itimestep = nstart+1, nstart+maxtimestep   
	   
        current_ts=current_ts+1
        if (mod(itimestep,print_step).eq.0) then
         write(*,*)'______________________________________________'
         write(*,*)'  current number of time step =',
     &           itimestep,'     current time=', real(time+dt)
         write(*,*)'______________________________________________'
        endif      
      
c     If not first time step, then update thermal energy, density and 
c     velocity half a time step  

        if (itimestep .ne. 1) then

          do i = 1, ntotal                        
            Tem_min(i) = Tem(i)
            temp_Con=0.                      
	      if (dim.eq.1) temp_Tem=0. 
            Tem(i) = Tem(i) + (dt/2.)* (dTem(i)+temp_Tem)	       
            if(Tem(i).lt.0)  Tem(i) = 0.
		           
		  if (.not.summation_density) then    
              density_min(i) = density(i)
	      temp_density=0.
	      if (dim.eq.1) temp_density=-nsym*density(i)*vx(1,i)/x(1,i)
               density(i)=density(i)+(dt/2.)*(ddensityh(i)+temp_density)
            endif 
           
            do d = 1, dim
              v_min(d, i) = vx(d, i)
              vx(d, i) = 0.5 
            enddo
          enddo 

        endif
c-----       
c---  Definition of variables out of the function vector:    
      
        call single_step(itimestep, dt, ntotal, hsml, volume, x, vx, u,  
     &       s, density, p, C, tdsdt, dx, dvx, du, ds, ddensity, itype,
     &       av, Tem, dTem)  
                  
        if (itimestep .eq. 1) then
       
          do i=1,ntotal

             temp_Tem=0.                           
   	      if (dim.eq.1) temp_Tem=0.        
            Tem(i) = Tem(i) + (dt/2.)*(dTem(i) + temp_Tem)
            if(Tem(i).lt.0)  Tem(i) = 0. 
 		           
            if (.not.summation_density ) then
	      temp_density=0.
	      if (dim.eq.1) temp_density=-nsym*density(i)*vx(1,i)/x(1,i)
              density(i)=density(i)+(dt/2.)*(ddensity(i)+temp_density)
            endif
         
            do d = 1, dim        
              vx(d, i) = 0.5 
              x(d, i) = x(d, i) + dt * vx(d, i)
            enddo  	            
          enddo 
                  
        else   
                    
          do i=1,ntotal            
              temp_Tem=0.                           
	      if (dim.eq.1) temp_Tem=0.                       
            Tem(i) = Tem_min(i) + dt*(dTem(i)+temp_Tem)
            if(Tem(i).lt.0)  Tem(i) = 0.          
           
            
		  if (.not.summation_density ) then 
              temp_densityh=0.
	      if (dim.eq.1) temp_density=-nsym*density(i)*vx(1,i)/x(1,i)        	           
              density(i) = density_min(i) + dt*(ddensity(i)+temp_density)
            endif
                  
            do d = 1, dim                   
              vx(d, i) = 0.5 
              x(d, i) = x(d, i) + dt * vx(d, i)
            enddo
          enddo
        
        endif 
     
        time = time + dt

	if (mod(itimestep,save_step).eq.0) then
		ngrab=ngrab+1
          write(supp,'(i4.4)') ngrab
		name='HydrauJump10_'//supp
		write(*,*) name
		open(23,file=name)
		call output(x,vx,volume,density,p,u,c,itype,hsml,
     &  		        ntotal,22,Tem)       
          close(23)
	endif 



      if (mod(itimestep,print_step).eq.0) then
		write(*,*)
          write(*,101) 'x', 'velocity', 'dvx'    
          write(*,100) x(1,moni_particle), vx(1,moni_particle), 
     &                 dvx(1,moni_particle)    
      endif


101   format(1x,3(2x,a12))	 
100   format(1x,3(2x,e12.6))
	 
      enddo

      nstart=current_ts

      end