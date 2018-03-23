      subroutine int_force(itimestep,dt,ntotal,hsml,volume,vx,niac,
     & density,eta,pair_i,pair_j,dwdx,u,itype,x,t,c,p,dvxdt,tdsdt,dedt)

c----------------------------------------------------------------------

c     itimestep: Current timestep number                            [in]
c     dt     :   Time step                                          [in]
c     ntotal : Number of particles                                  [in]
c     hsml   : Smoothing Length                                     [in]
c     volume : Particle volume                                      [in]
c     vx     : Velocities of all particles                          [in]
c     niac   : Number of interaction pairs                          [in]
c     density: water density                                        [in]
c     eta    : Dynamic viscosity                                    [in]
c     pair_i : List of first partner of interaction pair            [in]
c     pair_j : List of second partner of interaction pair           [in]
c     dwdx   : Derivative of kernel with respect to x, y and z      [in]
c     itype  : Type of particle (material types)                    [in]
c     u      : Particle internal energy                             [in]
c     x      : Particle coordinates                                 [in]
c     itype  : Particle type                                        [in]
c     t      : Particle temperature                             [in/out]
c     c      : Particle sound speed                                [out]
c     p      : Particle pressure                                   [out]
c     dvxdt  : Acceleration with respect to x, y and z             [out] 
c     tdsdt  : Production of viscous entropy                       [out]
c     dedt   : Change of specific internal energy                  [out]

      implicit none
      include 'param.inc.txt'
      
      integer itimestep, ntotal,niac,pair_i(max_interaction),
     &        pair_j(max_interaction), itype(maxn) 
      double precision dt, hsml(maxn), volume(maxn), vx(dim,maxn),      
     &       density(maxn), eta(maxn), dwdx(dim,max_interaction),
     &       u(maxn), x(dim,maxn), t(maxn), c(maxn), p(maxn),          
     &       dvxdt(dim,maxn), tdsdt(maxn), dedt(maxn)
      integer i, j, k, d
      double precision  dvx(dim), txx(maxn), tyy(maxn),
     &       tzz(maxn), txy(maxn), txz(maxn), tyz(maxn), vcc(maxn),
     &       hxx, hyy, hzz, hxy, hxz, hyz, hvcc, he, densityij, h


c     Initialization of acceleration 

      do i=1,ntotal          
		do d=1,dim
			dvxdt(d,i) = 0.e0
		enddo 
      enddo

c      do i=1,ntotal            
c     Pressure from equation of state

c		if (abs(itype(i)).eq.1) then
c			call p_gas(density(i), u(i), p(i),c(i))  
c		else if (abs(itype(i)).eq.2) then	     
c			call p_art_water(density(i), p(i), c(i))
c		endif  
c      enddo

      do k=1,niac
		i = pair_i(k)
		j = pair_j(k)  
          
c     For SPH algorithm 2
          
		if (pa_sph.eq.2) then 		
			do d=1,dim             
				h = -9.81 * dwdx(1,k)			 
				dvxdt(d,i) = dvxdt(d,i) + volume(j)*h
				dvxdt(d,j) = dvxdt(d,j) - volume(i)*h 
			enddo    		       
		endif  		      
      enddo

      end