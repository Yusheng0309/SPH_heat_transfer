      subroutine sum_density(ntotal,hsml,volume,niac,pair_i,pair_j,w,
     &           itype,density,x)

C----------------------------------------------------------------------
C     Subroutine to calculate the density with SPH summation algorithm.

C     ntotal : Number of particles                                  [in]
C     hsml   : Smoothing Length                                     [in]
C     volume : Particle volume                                      [in]
C     niac   : Number of interaction pairs                          [in]
C     pair_i : List of first partner of interaction pair            [in]
C     pair_j : List of second partner of interaction pair           [in]
C     w      : Kernel for all interaction pairs                     [in]
c     itype   : type of particles                                   [in]
c     x       : Coordinates of all particles                        [in]
c     density : water density                                      [out]
    
      implicit none
      include 'param.inc.txt'
      
      integer ntotal, niac, pair_i(max_interaction),
     &        pair_j(max_interaction), itype(maxn)  
      double precision hsml(maxn), volume(maxn), w(max_interaction),
     &       density(maxn) 
      integer i, j, k, d      
      double precision selfdens, hv(dim), r, wi(maxn), x(dim,maxn)     

c     wi(maxn)---integration of the kernel itself
        
      do d=1,dim
        hv(d) = 0.e0
      enddo

c     Self density of each particle: Wii (Kernel for distance 0)
c     and take contribution of particle itself:

      r=0.
      
c     Firstly calculate the integration of the kernel over the space

      do i=1,ntotal
        call kernel(r,hv,hsml(i),selfdens,hv)
        wi(i)=selfdens*volume(i)/density(i)
      enddo

      do k=1,niac
        i = pair_i(k)
        j = pair_j(k)
        wi(i) = wi(i) + volume(j)/density(j)*w(k)
        wi(j) = wi(j) + volume(i)/density(i)*w(k)
      enddo

c     Secondly calculate the density integration over the space

      do i=1,ntotal
        call kernel(r,hv,hsml(i),selfdens,hv)
        density(i) = selfdens*volume(i)
      enddo

c      Calculate SPH sum for density:
      do k=1,niac 	    
        i = pair_i(k)
        j = pair_j(k)
        density(i) = density(i) + volume(j)*w(k)
        density(j) = density(j) + volume(i)*w(k)
      enddo

c     Thirdly, calculate the normalized density, density=sum(density)/sum(w)
     
      if (nor_density) then 
        do i=1, ntotal
          density(i)=density(i)/wi(i)
        enddo
      endif 
      
      end
      
      subroutine con_density(ntotal,volume,niac,pair_i,pair_j,
     &           dwdx,vx,itype,x,density, ddensitydt)

c----------------------------------------------------------------------
c     Subroutine to calculate the density with SPH continuiity approach.

c     ntotal : Number of particles                                  [in]
c     volume : Particle volume                                      [in]
c     niac   : Number of interaction pairs                          [in]
c     pair_i : List of first partner of interaction pair            [in]
c     pair_j : List of second partner of interaction pair           [in]
c     dwdx   : derivation of Kernel for all interaction pairs       [in]
c     vx     : Velocities of all particles                          [in]
c     itype  : type of particles                                    [in]
c     x      : Coordinates of all particles                         [in]
c     density: water density                                        [in]
c     ddensitydt : water density change rate of each particle      [out]   

      implicit none
      include 'param.inc.txt'
      
      integer ntotal,niac,pair_i(max_interaction),
     &        pair_j(max_interaction), itype(maxn)    
      double precision volume(maxn), dwdx(dim, max_interaction),
     &       vx(dim,maxn), x(dim,maxn), density(maxn), ddensitydt(maxn)
      integer i,j,k,d    
      double precision    vcc, dvx(dim) 
      
      do i = 1, ntotal
        ddensitydt(i) = 0.
      enddo
     
      do k=1,niac      
        i = pair_i(k)
        j = pair_j(k)
        do d=1,dim
          dvx(d) = vx(d,i) - vx(d,j) 
        enddo        
        vcc = dvx(1)*dwdx(1,k)        
        do d=2,dim
          vcc = vcc + dvx(d)*dwdx(d,k)
        enddo    
        ddensitydt(i) = ddensitydt(i) + volume(j)*vcc
        ddensitydt(j) = ddensitydt(j) + volume(i)*vcc       
      enddo    
	 
      end