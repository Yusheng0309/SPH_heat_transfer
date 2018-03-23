      subroutine h_upgrade(dt,ntotal, volume, vx, density, niac, 
     &           pair_i, pair_j, dwdx, hsml)

c-----------------------------------------------------------------------
c     Subroutine to evolve smoothing length

c     dt     : time step                                            [in]
c     ntotal : Number of particles                                  [in]
c     volume : Particle volume                                      [in]
c     vx     : Velocities of all particles                          [in]
c     density: water density                                        [in]
c     niac   : Number of interaction pairs                          [in]
c     pair_i : List of first partner of interaction pair            [in]
c     pair_j : List of second partner of interaction pair           [in]
c     dwdx   : Derivative of kernel with respect to x, y and z      [in]
c     hsml   : Smoothing Length                                 [in/out]
   
      implicit none
      include 'param.inc.txt'
      
      integer ntotal, niac, pair_i(max_interaction), 
     &        pair_j(max_interaction)
      double precision volume(maxn), vx(dim, maxn), density(maxn),
     &       dwdx(dim, max_interaction), hsml(maxn)     
      integer i,j,k,d
      double precision dt, fac, dvx(dim), hvcc, vcc(maxn), dhsml(maxn)     

      if (sle.eq.0 ) then     

c---  Keep smoothing length unchanged. 
     
        return
      
      else if (sle.eq.2) then
      
c---  dh/dt = (-1/dim)*(h/density)*(ddensity/dt).

        do i=1,ntotal
          vcc(i) = 0.e0
        enddo
      
        do k=1,niac
          i = pair_i(k)
          j = pair_j(k)
          do d=1,dim
            dvx(d) = vx(d,j) - vx(d,i) 
          enddo
          hvcc = dvx(1)*dwdx(1,k)
          do d=2,dim
            hvcc = hvcc + dvx(d)*dwdx(d,k)
          enddo    
          vcc(i) = vcc(i) + volume(j)*hvcc/density(j)
          vcc(j) = vcc(j) + volume(i)*hvcc/density(i)         
        enddo  
        
        do i = 1, ntotal
          dhsml(i) = (hsml(i)/dim)*vcc(i)
          hsml(i) = hsml(i) + dt*dhsml(i)      
          if (hsml(i).le.0) hsml(i) = hsml(i) - dt*dhsml(i) 
        enddo
    
      else if(sle.eq.1) then
            
c        fac = 2.0
        do i = 1, ntotal          
c         hsml(i) = fac * (volume(i)/density(i))**(1./dim)
		hsml(i) = 1.2*0.5 * (0.031 / density(i))
        enddo
       
      endif 
       
      end 