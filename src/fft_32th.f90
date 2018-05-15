subroutine fft_32th( t, pv )
! this is the analytical fourier expression for the 32-th orbit 
! keep four digits 

implicit none
  
! Input  and Output Declaration   
real(kind=dp), intent(in)        ::  t  
real(kind=dp), intent(out)        ::   pv(6) 
 
! Local Variable
real(kind=dp)  :: x, y, z, vx, vy, vz, feq , &
                  y0, z0, ampx, ampy, ampz, ampvx, ampvz 

feq =  1.802491547458d0 

! cosine
ampx = 0.235929626143d0

! sine 
y0 = 0.757312207266d0
ampy = 0.217866843104d0 

! sine 
z0 = 0.817078194813d0
ampz = 0.283714650034

! sine 
ampvx = 0.425261156918d0 

! cosine
ampvy = 0.392703143166d0

! cosine
ampvz = 0.511393258576d0

x   =  ampx * dcos( feq * t)
y   =  y0 + ampy * dsin( feq * t) 
z   =  z0 + ampz * dsin( feq * t) 

vx   =    ampvx * dsin( feq * t)  
vy   =    ampvy * dcos( feq * t)  
vz   =    ampvz * dcos( feq * t)     

pv = (/x, y, z, vx, vy, vz/)  

return 
end subroutine fft_32th
