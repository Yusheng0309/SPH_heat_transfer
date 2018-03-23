      subroutine trans_temperature(x,vx,volume, density, niac, pair_i,
     &          pair_j,w,dwdx, hsml, ntotal, Tem, dTem,itimestep)
     
c----------------------------------------------------------------------
c      x      :coordinates of particles                       [in/out]
c      vx     :velocities of particles                        [in/out]
c      mass   :mass of particles                                  [in]
c      density:water density                                  [in/out]
c      u      :internal energy of particles                   [in/out]
c      Tem    :temperature                                    [in/out]
c      w      : Kernel for all interaction pairs                  [in]
c      dwdx   : Derivative of kernel with respect to x, y and z   [in]
c      hsml   :smoothing lengths of particles                 [in/out]
c      niac   : Number of interaction pairs                       [in]
c      ntotal :total particle number                              [in]  
  
      implicit none     
      include 'param.inc.txt'
      
      integer itype(maxn), ntotal, niac,pair_i(max_interaction),
     &        pair_j(max_interaction),itimestep
      double precision x(dim, maxn), vx(dim, maxn), volume(maxn), 
     &       density(maxn), w(max_interaction),
     &       hsml(maxn), dwdx(dim,max_interaction), Tem(maxn) 
      integer i, j, k, d       
      double precision dC(maxn),  av(dim, maxn)
      double precision  xl, Fij, rij, mhsml, rdwdx, dTTem(maxn), 
     &       dxiac(dim), mp, m, hv(dim), wi(maxn),driac,space_x, Di,
     &       dTem(maxn)

c     thermal conductivity
      parameter( Di =100.0) 
	

      space_x = 100.0 


c     wi(maxn)---integration of the kernel itself

	do d=1,dim
        hv(d) = 0.e0
      enddo 

      do i=1,ntotal
	  dTem(i) = 0.e0
	 dTTem(i) = 0.e0
      enddo


c     Secondly calculate the rho integration over the space
      do k=1,niac
        i = pair_i(k)
        j = pair_j(k)
        mhsml= (hsml(i)+hsml(j))/2.
	  driac = 0.e0
        rdwdx = 0.e0
        do d=1,dim
          dxiac(d) = x(d,i) -  x(d,j)
          driac = driac + dxiac(d)*dxiac(d)
          rdwdx  = rdwdx + dxiac(d)*dwdx(d,k)    
        enddo
	  Fij = rdwdx/(driac+0.01*mhsml**2)
        !note  the spatial derivative with respect to coordinates of particle i or j
	  !Fij is not the same when applying KGC                                 
	dTTem(i) = dTTem(i) + volume(j)/density(j)*2*Di*(Tem(i)-Tem(j))*Fij  
      dTTem(j) = dTTem(j) + volume(i)/density(i)*2*Di*(Tem(j)-Tem(i))*Fij 
	enddo
   
      do i=1, ntotal
	  dTem(i) = volume(i)/(space_x*density(i))*dTTem(i)
      enddo
   	 
      end

