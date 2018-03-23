      subroutine single_step(itimestep, dt, ntotal, hsml, volume, x, vx,  
     &           u,s,density,p,t,tdsdt, dx, dvx, du, ds, ddensity,itype,
     &           av,Tem, dTem) 

c----------------------------------------------------------------------
c     Subroutine to determine the right hand side of a differential 
c     equation in a single step for performing time integration 
c     In this routine and its subroutines the SPH algorithms are performed.

c     itimestep: Current timestep number                            [in]
c     dt       : Timestep                                           [in]
c     ntotal   :  Number of particles                               [in]
c     hsml     :  Smoothing Length                                  [in]
c     volume   :  Particle volume                                   [in]
c     x        :  Particle position                                 [in]
c     vx       :  Particle velocity                                 [in]
c     u        :  Particle internal energy                          [in]
c     s        :  Particle entropy (not used here)                  [in]
c     density  :  water density                                 [in/out]
c     p        :  Pressure                                         [out]
c     Tem      :  Temperature                                   [in/out]
c     tdsdt    :  Production of viscous entropy t*ds/dt            [out]
c     dx       :  dx = vx = dx/dt                                  [out]
c     dvx      :  dvx = dvx/dt, force per unit volume              [out]
c     du       :  du  = du/dt                                      [out]
c     ds       :  ds  = ds/dt                                      [out]     
c     ddensity :  ddensity = ddensity/dt                           [out]
c     itype    :  Type of particle                                  [in]
c     av       :  Monaghan average velocity                        [out]

      implicit none
      include 'param.inc.txt'

      integer itimestep, ntotal, itype(maxn)    
      double precision dt, hsml(maxn), volume(maxn), x(dim,maxn),
     &       vx(dim,maxn), u(maxn), s(maxn), density(maxn), p(maxn), 
     &       t(maxn), tdsdt(maxn), dx(dim,maxn), dvx(dim,maxn), 
     &       du(maxn), ds(maxn), ddensity(maxn), av(dim, maxn),
     &       dTem(maxn), Tem(maxn), dTemi(maxn)            
      integer i, d, nvirt, niac, pair_i(max_interaction),
     &        pair_j(max_interaction), ns(maxn)
      double precision w(max_interaction), dwdx(dim,max_interaction),  
     &       indvxdt(dim,maxn),exdvxdt(dim,maxn),ardvxdt(dim,maxn),  
     &       avdudt(maxn), ahdudt(maxn), c(maxn),  eta(maxn),
     +       frdvxdt(dim,maxn)  
      double precision manning_n, Sfx_i, g

	manning_n = 0.0278 
	g = 9.81                             

      do  i=1,ntotal
        avdudt(i) = 0.
        ahdudt(i) = 0.
        do  d=1,dim
          indvxdt(d,i) = 0.
          ardvxdt(d,i) = 0.
          exdvxdt(d,i) = 0.
c	    frdvxdt(d,i) = 0.
        enddo
      enddo  
 
c---  Positions of virtual (boundary) particles: 

      nvirt = 0
      if (virtual_part) then 
        call virt_part(itimestep, ntotal,nvirt,hsml,volume,x,vx,
     &       density,u,p,itype,Tem)
      endif 
     
c---  Interaction parameters, calculating neighboring particles
c     and optimzing smoothing length

      if (nnps.eq.1) then 
        call direct_find(itimestep, ntotal+nvirt,hsml,x,niac,pair_i,
     &       pair_j,w,dwdx,ns)
      else if (nnps.eq.2) then
        call link_list(itimestep, ntotal+nvirt,hsml(1),x,niac,pair_i,
     &       pair_j,w,dwdx,ns)
      endif         
                        
c---  Density approximation or change rate
      
      if (summation_density) then      
        call sum_density(ntotal+nvirt,hsml,volume,niac,pair_i,pair_j,w,
     &       itype,density,x)          
      else             
        call con_density(ntotal+nvirt,volume,niac,pair_i,pair_j,
     &       dwdx,vx, itype,x,density, ddensity)         
      endif

c---  Dynamic viscosity:

      if (visc) call viscosity(ntotal+nvirt,itype,x,density,eta)
       
c---  Internal forces:

	call int_force(itimestep,dt,ntotal,hsml,volume,vx,niac,density,
     &     eta,pair_i,pair_j,dwdx,u,itype,x,t,c,p,indvxdt,tdsdt,du)
                  
c---  Artificial viscosity:

      if (visc_artificial) call art_visc(ntotal+nvirt,hsml,
     &   volume,x,vx,niac,density,c,pair_i,pair_j,w,dwdx,ardvxdt,avdudt)
      
c---  External forces:

      if (ex_force) call ext_force(ntotal+nvirt,volume,x,niac,
     &                   pair_i,pair_j,itype, hsml, exdvxdt)

c     Calculating the neighboring particles and undating HSML
      
      if (sle.ne.0) call h_upgrade(dt,ntotal, volume, vx, density, niac, 
     &                   pair_i, pair_j, dwdx, hsml)

      if (heat_artificial) call art_heat(ntotal+nvirt,hsml,
     &   volume,x,vx,niac,density,u, c,pair_i,pair_j,w,dwdx,ahdudt)

	if (trans_tem) call trans_temperature(x,vx, volume, density, niac, 
     &	    pair_i,pair_j, w, dwdx, hsml, ntotal, Tem, dTem,itimestep)
     
c     Calculating average velocity of each partile for avoiding penetration

      if (average_velocity) call av_vel(ntotal,volume,niac,pair_i,
     &                           pair_j, w, vx, density, av)

c---  Convert velocity and force to f and dfdt 
 	   
	do i=1,ntotal

c		if (x(1,i).gt.0)
c    +		call friction(ntotal,hsml,volume,x,vx,niac,density,c,
c    &             pair_i,pair_j,w,dwdx,frdvxdt) 
		if (x(1,i).gt.0)
     +		Sfx_i = manning_n**2 * vx(1,i)*sqrt( vx(1,i)**2 )  
     +		    / (density(i)**(4./3.))

 
		dvx(1,i) = 0 !indvxdt(1,i) + exdvxdt(1,i) + ardvxdt(1,i)
     +               !+ (- g * Sfx_i) 
c	    dTem(i) =  dTemi(i) 
	enddo	

      if (mod(itimestep,print_step).eq.0) then      
          write(*,*)
          write(*,*) '**** Information for particle ****', 
     &        		moni_particle         
          write(*,101)'internal a ','artifical a=',
     &         		'external a ','frictional a=',
     +                'total a '   
          write(*,100)indvxdt(1,moni_particle),ardvxdt(1,moni_particle),
     &                exdvxdt(1,moni_particle),
     +                dvx(1,moni_particle)-ardvxdt(1,moni_particle)-
     +                exdvxdt(1,moni_particle),
     +                dvx(1,moni_particle)          
      endif

101   format(1x,5(2x,a12))      
100   format(1x,5(2x,e12.6))      

      end