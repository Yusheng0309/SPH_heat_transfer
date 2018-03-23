      program SPH

c----------------------------------------------------------------------
c     This is a three dimensional SPH code. the followings are the 
c     basic parameters needed in this codeor calculated by this code

c     volume :volume of particles                                 [in]
c     ntotal :total particle number ues                           [in]
c     dt     :time step used in the time integration              [in]
c     itype  :types of particles                                  [in]
c     x      :coordinates of particles                        [in/out]
c     vx     :velocities of particles                         [in/out]
c     density:water density of particles                      [in/out]
c     p      :pressure  of particles                          [in/out]
c     u      :internal energy of particles                    [in/out]
c     hsml   :smoothing lengths of particles                  [in/out]
c     c      :sound velocity of particles                        [out]
c     s      :entropy of particles                               [out]
c     e      :total energy of particles                          [out]

      implicit none     
      include 'param.inc.txt'

      integer ntotal, itype(maxn), maxtimestep, d, m, i, yesorno      
      double precision x(dim,maxn), vx(dim, maxn), volume(maxn),
     &                 density(maxn), p(maxn), u(maxn), c(maxn), 
     &                 s(maxn),e(maxn), hsml(maxn), dt, Tem(maxn)
      double precision s1, s2

      call time_print
      call time_elapsed(s1)      

      if (Hydrau_Jump) dt= 1.0 

      call input(x, vx, volume, density, p, u, itype, hsml, ntotal,Tem)   
	  
 1    write(*,*)'  ***************************************************' 
      write(*,*)'          Please input the maximal time steps '
      write(*,*)'  ***************************************************'
      read(*,*) maxtimestep 
	     
      call time_integration(x, vx, volume, density, p,Tem, u, c, s, e, 
     & itype, hsml, ntotal, maxtimestep, dt )
      call output(x, vx, volume, density, p, u, c, itype, hsml, ntotal,
     &            22, Tem) 
	     
      write(*,*)'  ***************************************************'
      write(*,*)' Are you going to run more time steps ? (0=No, 1=yes)'
      write(*,*)'  ***************************************************'
      read (*,*) yesorno
	     
      if (yesorno.ne.0) go to 1
      call time_print
      call time_elapsed(s2)      
      write (*,*)'        Elapsed CPU time = ', s2-s1
                           
      end