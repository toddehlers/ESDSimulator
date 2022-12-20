
subroutine therm_output (time,temperature,ntime,nsurf,xdepth,ydepth,zdepth,run,nstep,istep,nrun)

  ! Subroutine for Pecube that writes the thermal history for every node
  ! at the specified time steps in Pecube.in

  integer nsurf,ntime,nstep,istep
  real*4 time(nsurf,ntime),temperature(nsurf,ntime)
  real*4 xdepth(nsurf,ntime),ydepth(nsurf,ntime),zdepth(nsurf,ntime)
  character run*100,count3*4,count4*4,nsurfc*5
  integer head1,head2,nrun

  ! Finds the length of time for the current output (head1) and the physical time
  ! at which this output is from (head2) 
  head1=nint(time(1,1)-time(1,ntime))
  head2=nint(time(1,ntime))

  write (count3,'(i4)') head1
  write (count4,'(i4)') head2
  write (nsurfc,'(i5)') nsurf
  if (head1.lt.1000) count3(1:1)='_'
  if (head1.lt.100) count3(1:2)='_0'
  if (head1.lt.10) count3(1:3)='_00'
  if (head2.lt.1000) count4(1:1)='_'
  if (head2.lt.100) count4(1:2)='_0'
  if (head2.lt.10) count4(1:3)='_00'
  if (nsurf.lt.10000) nsurfc(1:1)='0'
  if (nsurf.lt.1000) nsurfc(1:2)='00'
  if (nsurf.lt.100) nsurfc(1:3)='000'
  if (nsurf.lt.10) nsurfc(1:4)='0000'

  open (94,file=run(1:nrun)//'/X_hist'//count3//'My_surf'//count4//'Ma.dat',status='unknown')
  open (95,file=run(1:nrun)//'/Y_hist'//count3//'My_surf'//count4//'Ma.dat',status='unknown')
  open (96,file=run(1:nrun)//'/Z_hist'//count3//'My_surf'//count4//'Ma.dat',status='unknown')
  open (97,file=run(1:nrun)//'/Time_hist'//count3//'My_surf'//count4//'Ma.dat',status='unknown')
  open (98,file=run(1:nrun)//'/Temp_hist'//count3//'My_surf'//count4//'Ma.dat',status='unknown')

  ! Prints values to output files
  do k=1,ntime
    write (94,'('//nsurfc//'f12.3)') xdepth(:,k)
    write (95,'('//nsurfc//'f12.3)') ydepth(:,k)
    write (96,'('//nsurfc//'f12.3)') zdepth(:,k)
    write (97,'('//nsurfc//'f12.2)') time(:,k)
    write (98,'('//nsurfc//'f12.1)') temperature(:,k)
  enddo

  close(94)
  close(95)
  close(96)
  close(97)
  close(98)

  return
end