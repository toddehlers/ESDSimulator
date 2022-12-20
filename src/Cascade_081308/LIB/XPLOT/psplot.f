c
c subroutine equivalent to xplot to generate PostScript
c
c------------------------------------------------------------------

      subroutine Popen (window,window_length,fnme,nfnme)

      common /window_geometry/ iwind(4)
      common /factors/  ax,ay,bx,by

      integer*4       window(2),window_length
      character*256   fnme,fnmep

      if (nfnme.eq.-1) goto 1000

        do i=1,256
        fnmep(i:i)=' '
        enddo

        do i=1,nfnme
        fnmep(i:i)=fnme(i:i)
        enddo

      iwind(1)=(1200-window(1))/2
      iwind(2)=(1024-window(2))/2
      iwind(3)=window(1)
      iwind(4)=window(2)

      open (77,file=fnmep,status='unknown')

 1000 write (77,'(a)') '%!XPlot'

      xmin=0.
      xmax=iwind(3)
      ymin=0.
      ymax=iwind(4)

      xmi=((29.7-21.*1200./1024.)/2.+iwind(1)/1024.*21.)
      xma=xmi+iwind(3)/1024.*21.
      ymi=(1024.-iwind(2)-iwind(4))/1024.*21.
      yma=(1024.-iwind(2))/1024.*21.

      write (77,'(a2,4e12.6)') '%%',xmi,xma,ymi,yma
      write (77,*) 'gsave initgraphics'
      write (77,'(a)') '90 rotate'
      write (77,'(a)') '0 -595 translate'

c      asprat=xmax/ymax
c      asprat0=27./19.
c      if (abs(asprat).gt.asprat0) then
c      xmi=1.5
c      xma=28.5
c      ymi=1.+(19.-19.*asprat0/asprat)/2.
c      yma=ymi+19.*asprat0/asprat
c        else
c        xmi=1.5+(27.-27.*asprat/asprat0)/2.
c        xma=xmi+27.*asprat/asprat0
c        ymi=1.
c        yma=20.
c        endif

      if (nfnme.ne.-1) then
      ax=(xma-xmi)/(xmax-xmin)*72./2.54
      ay=(yma-ymi)/(ymax-ymin)*72./2.54
      bx=(xmi-xmin*(xma-xmi)/(xmax-xmin))*72./2.54
      by=(ymi-ymin*(yma-ymi)/(ymax-ymin))*72./2.54
      endif

      write (77,*) '/ax ',ax,' def'
      write (77,*) '/ay ',ay,' def'
      write (77,*)
     &  '/bx ',bx,' def'
      write (77,*)
     &  '/by ',by,' def'
      write (77,*)
     & '/m { ay mul by add exch ax mul bx add exch moveto} def'
      write (77,*)
     & '/l { ay mul by add exch ax mul bx add exch lineto} def'

      write (77,*) '/s {stroke} def'
      write (77,*) 
     & '/f1 {/Times-Roman findfont 5 scalefont setfont} def'
      write (77,*) 
     & '/f2 {/Times-Roman findfont 7 scalefont setfont} def'
      write (77,*) 
     & '/f3 {/Times-Roman findfont 9 scalefont setfont} def'
      write (77,*) 
     & '/f4 {/Times-Roman findfont 14 scalefont setfont} def'

      write (77,'(a)') 
     & '/sf1 {/Helvetica findfont 6 scalefont setfont} def'
      write (77,'(a)') 
     & '/sf2 {/Helvetica findfont 8 scalefont setfont} def'
      write (77,'(a)') 
     & '/sf3 {/Helvetica findfont 10 scalefont setfont} def'
      write (77,'(a)') 
     & '/sf4 {/Helvetica findfont 12 scalefont setfont} def'
      write (77,'(a)') 
     & '/sf5 {/Helvetica findfont 14 scalefont setfont} def'
      write (77,'(a)') 
     & '/sf6 {/Helvetica findfont 18 scalefont setfont} def'
      write (77,'(a)') 
     & '/sf7 {/Helvetica findfont 24 scalefont setfont} def'
      write (77,'(a)') 
     & '/sf8 {/Helvetica findfont 48 scalefont setfont} def'

      write (77,'(a)') 
     & '/ssf1 {/Symbol findfont 6 scalefont setfont} def'
      write (77,'(a)') 
     & '/ssf2 {/Symbol findfont 8 scalefont setfont} def'
      write (77,'(a)') 
     & '/ssf3 {/Symbol findfont 10 scalefont setfont} def'
      write (77,'(a)') 
     & '/ssf4 {/Symbol findfont 12 scalefont setfont} def'
      write (77,'(a)') 
     & '/ssf5 {/Symbol findfont 14 scalefont setfont} def'
      write (77,'(a)') 
     & '/ssf6 {/Symbol findfont 18 scalefont setfont} def'
      write (77,'(a)') 
     & '/ssf7 {/Symbol findfont 24 scalefont setfont} def'
      write (77,'(a)') 
     & '/ssf8 {/Symbol findfont 48 scalefont setfont} def'

      write (77,*) '/p {show} def '
      write (77,*) '/n {newpath} def'
      write (77,*) '/cp {closepath} def'
      write (77,*) '/sd {setdash} def'
      write (77,*) '/f {fill} def'
      write (77,*) '/sw {setlinewidth} def'
      write (77,*) '/sg {setrgbcolor} def'
      write (77,*) '/a {arc} def'
      write (77,*) '0 sw'
      write (77,*) 'initclip'
      write (77,*) 'n'

      if (nfnme.eq.-1) then
      aax=(xma-xmi)/(xmax-xmin)*72./2.54
      aay=(yma-ymi)/(ymax-ymin)*72./2.54
      bbx=(xmi-xmin*(xma-xmi)/(xmax-xmin))*72./2.54
      bby=(ymi-ymin*(yma-ymi)/(ymax-ymin))*72./2.54
      write (77,*) '/ax ',aax,' def'
      write (77,*) '/ay ',aay,' def'
      write (77,*)
     &  '/bx ',bbx,' def'
      write (77,*)
     &  '/by ',bby,' def'
      endif

      write (77,*) 1,1,1,' sg'
      write (77,*) 0,0,' m'
      write (77,*) iwind(3),0,' l'
      write (77,*) iwind(3),iwind(4),' l'
      write (77,*) 0,iwind(4),' l'
      write (77,*) 0,0,' l'
      write (77,*) 'cp s'
c from here
      write (77,*) 'n'
      write (77,*) 0,0,' m'
      write (77,*) iwind(3),0,' l'
      write (77,*) iwind(3),iwind(4),' l'
      write (77,*) 0,iwind(4),' l'
      write (77,*) 0,0,' l'
      write (77,*) 'cp clip s'
      write (77,*) 0,0,0,' sg'
c to here
      if (nfnme.eq.-1) then
      write (77,*) '/ax ',ax,' def'
      write (77,*) '/ay ',ay,' def'
      write (77,*)
     &  '/bx ',bx,' def'
      write (77,*)
     &  '/by ',by,' def'
      endif

      return
      end

c------------------------------------------------------------------

      subroutine Pclose

      write (77,'(a)') 's'
      write (77,'(a)')  'showpage'
      write (77,'(a)')  'grestore'
      write (77,'(a)') '%%EOF'

      close (77)

      return
      end

c------------------------------------------------------------------

      subroutine Pscale (ixmin,ixmax,iymin,iymax,
     &                   xmin,xmax,ymin,ymax,is)

      integer*4  window(4)

      common /window_geometry/ window
      common /factors/  ax,ay,bx,by

      if (is.eq.1) then

      win1=window(1)+ixmin
      win3=ixmax-ixmin
      win2=window(2)+iymin
      win4=iymax-iymin

      xxmin=xmin
      xxmax=xmax
      yymin=ymax
      yymax=ymin

c      xxmin=xmin+(xmax-xmin)*(0-ixmin)/(ixmax-ixmin)
c      xxmax=xmin+(xmax-xmin)*(window(3)-ixmin)/(ixmax-ixmin)
c      yymax=ymin+(ymax-ymin)*(0-iymin)/(iymax-iymin)
c      yymin=ymin+(ymax-ymin)*(window(4)-iymin)/(iymax-iymin)

c      asprat=abs((xxmax-xxmin)/(yymax-yymin))
c      asprat0=27./19.
c      if (asprat.gt.asprat0) then
c      xmi=1.5
c      xma=28.5
c      ymi=1.+(19.-19.*asprat0/asprat)/2.
c      yma=ymi+19.*asprat0/asprat
c        else
c        xmi=1.5+(27.-27.*asprat/asprat0)/2.
c        xma=xmi+27.*asprat/asprat0
c        ymi=1.
c        yma=20.
c        endif

      xmi=((29.7-21.*1200./1024.)/2.+win1/1024.*21.)
      xma=xmi+win3/1024.*21.
      ymi=(1024.-win2-win4)/1024.*21.
      yma=(1024.-win2)/1024.*21.

      write (77,*) '/ax ',(xma-xmi)/(xxmax-xxmin)*72./2.54,' def'
      write (77,*) '/ay ',(yma-ymi)/(yymax-yymin)*72./2.54,' def'
      write (77,*)
     &  '/bx ',(xmi-xxmin*(xma-xmi)/(xxmax-xxmin))*72./2.54,' def'
      write (77,*)
     &  '/by ',(ymi-yymin*(yma-ymi)/(yymax-yymin))*72./2.54,' def'

      ax=(xma-xmi)/(xxmax-xxmin)*72./2.54
      ay=(yma-ymi)/(yymax-yymin)*72./2.54
      bx=(xmi-xxmin*(xma-xmi)/(xxmax-xxmin))*72./2.54
      by=(ymi-yymin*(yma-ymi)/(yymax-yymin))*72./2.54

      else

      xmin=window(1)
      xmax=window(1)+window(3)
      ymin=window(2)
      ymax=window(2)+window(4)

c      asprat=(xmax-xmin)/(ymax-ymin)
c      asprat0=27./19.
c      if (abs(asprat).gt.asprat0) then
c      xmi=1.5
c      xma=28.5
c      ymi=1.+(19.-19.*asprat0/asprat)/2.
c      yma=ymi+19.*asprat0/asprat
c        else
c        xmi=1.5+(27.-27.*asprat/asprat0)/2.
c        xma=xmi+27.*asprat/asprat0
c        ymi=1.
c        yma=20.
c        endif

      xmi=((29.7-21.*1200./1024.)/2.+window(1)/1024.*21.)
      xma=xmi+window(3)/1024.*21.
      ymi=(1024.-window(2)-window(4))/1024.*21.
      yma=(1024.-window(2))/1024.*21.

      write (77,*) '/ax ',(xma-xmi)/(xmax-xmin)*72./2.54,' def'
      write (77,*) '/ay ',(yma-ymi)/(ymax-ymin)*72./2.54,' def'
      write (77,*)
     &  '/bx ',(xmi-xmin*(xma-xmi)/(xmax-xmin))*72./2.54,' def'
      write (77,*)
     &  '/by ',(ymi-ymin*(yma-ymi)/(ymax-ymin))*72./2.54,' def'
      write (77,*)
     & '/m { ay mul by add exch ax mul bx add exch moveto} def'
      write (77,*)
     & '/l { ay mul by add exch ax mul bx add exch lineto} def'

      ax=(xma-xmi)/(xmax-xmin)*72./2.54
      ay=(yma-ymi)/(ymax-ymin)*72./2.54
      bx=(xmi-xmin*(xma-xmi)/(xmax-xmin))*72./2.54
      by=(ymi-ymin*(yma-ymi)/(ymax-ymin))*72./2.54

      endif

      return
      end

c------------------------------------------------------------------

      subroutine Pplot (x,y,ip)

      if (ip.eq.3) then
      write (77,*) x,y,' m'
      else
      write (77,*) x,y,' l'
      endif

      return
      end

c--------------------------------------------------------------------

      subroutine Psave

      write (77,*) 's'

      return
      end

c--------------------------------------------------------------------

      subroutine Pflush

      write (77,*) 's'

      return
      end

c-------------------------------------------------------------------

      subroutine Pline (x,y,n,ip)
      
      real*4  x(n),y(n)

      write (77,*) 's'
      if (ip.eq.4) write (77,*) '[3 3] 0 sd'
      if (ip.eq.5) write (77,*) '[6 3] 0 sd'
      if (ip.eq.6) write (77,*) '[6 6] 0 sd'
      if (ip.eq.7) write (77,*) '[9 3] 0 sd'
      if (ip.eq.8) write (77,*) '[9 3 3 3] 0 sd'
      write (77,*) x(1),y(1),' m'
        do i=2,n
        write (77,*) x(i),y(i),' l'
        enddo
      write (77,*) 's'
      if (ip.gt.3) write (77,*) '[] 0 sd'

      return
      end

c-------------------------------------------------------------------

      subroutine Pfill (x,y,n)
      
      real*4  x(n),y(n)

      write (77,*) 's'
      write (77,*) x(1),y(1),' m'
        do i=2,n
        write (77,*) x(i),y(i),' l'
        enddo
      write (77,*) 'f s'

      return
      end

c--------------------------------------------------------------------

      subroutine Pthick (ithick)

      write (77,*) 's ',0.3*ithick,' sw'

      return
      end

c----------------------------------------------------------------------

      subroutine Pcmap (red,green,blue)

      integer*4      red(230),green(230),blue(230)

      common  /PS_color/  valr(230),valg(230),valb(230)

        do i=1,230
        valr(i)=float(red(i))/256.
        valg(i)=float(green(i))/256.
        valb(i)=float(blue(i))/256.
        enddo

      return
      end

c----------------------------------------------------------------------

      subroutine Ppen (ipen)

      common  /PS_color/  valr(230),valg(230),valb(230)

      write (77,*) 's ',valr(ipen),valg(ipen),valb(ipen),' sg'

      return
      end

c--------------------------------------------------------------------

      subroutine Psymbol (x,y,is,str,nstr)

      character*256 str,strp

      ip=0
        do i=1,nstr
        ip=ip+1
          if (str(i:i).eq.')' .or. str(i:i).eq.'(') then
          strp(ip:ip)=char(92)
          ip=ip+1
          endif
        strp(ip:ip)=str(i:i)
        enddo
      nstrp=ip

      write (77,*) x,y,' m'
      if (is.eq.1) write (77,*) 'f1'
      if (is.eq.2) write (77,*) 'f2'
      if (is.eq.3) write (77,*) 'f3'
      if (is.eq.4) write (77,*) 'f4'
      write (77,*) '('//strp(1:nstrp)//') p s'

      return
      end

c--------------------------------------------------------------------

      subroutine Psupersymbol (x,y,is,str,nstr)

      character*256 str,strp

      ip=0
        do i=1,nstr
        ip=ip+1
          if (str(i:i).eq.')' .or. str(i:i).eq.'(') then
          strp(ip:ip)=char(92)
          ip=ip+1
          endif
        strp(ip:ip)=str(i:i)
        enddo
      nstrp=ip

      write (77,*) x,y,' m'
      if (is.eq.1) write (77,*) 'sf1'
      if (is.eq.2) write (77,*) 'sf2'
      if (is.eq.3) write (77,*) 'sf3'
      if (is.eq.4) write (77,*) 'sf4'
      if (is.eq.5) write (77,*) 'sf5'
      if (is.eq.6) write (77,*) 'sf6'
      if (is.eq.7) write (77,*) 'sf7'
      if (is.eq.8) write (77,*) 'sf8'
      if (is.eq.-1) write (77,*) 'ssf1'
      if (is.eq.-2) write (77,*) 'ssf2'
      if (is.eq.-3) write (77,*) 'ssf3'
      if (is.eq.-4) write (77,*) 'ssf4'
      if (is.eq.-5) write (77,*) 'ssf5'
      if (is.eq.-6) write (77,*) 'ssf6'
      if (is.eq.-7) write (77,*) 'ssf7'
      if (is.eq.-8) write (77,*) 'ssf8'
      write (77,*) '('//strp(1:nstrp)//') p s'

      return
      end

c---------------------------------------------------------------------

      subroutine Pclear

      integer*4   window_dum(2)
      character*256 dumc

      rewind 77
      call Popen (window_dum,idum,dumc,-1)

      return
      end

c----------------------------------------------------------------------

      subroutine Pclip (xmin,xmax,ymin,ymax)

      write (77,*) 'gsave'
      write (77,*) 'newpath'
      write (77,*) xmin,ymin,' m'
      write (77,*) xmin,ymax,' l'
      write (77,*) xmax,ymax,' l'
      write (77,*) xmax,ymin,' l'
      write (77,*) xmin,ymin,' l'
      write (77,*) 'closepath clip'

      return
      end

c----------------------------------------------------------------------

      subroutine Pclipoff ()

      write (77,*) 'grestore'

      return
      end
      
c----------------------------------------------------------------------

      subroutine Pcircle (xc,yc,radius)

      common /factors/  ax,ay,bx,by

      radiusx=radius*ax
      radiusy=radius*ay
      write (77,*) xc*ax+bx,yc*ay+by,
     &             (radiusx+radiusy)/2.,0,360,' a s'

      return
      end

c----------------------------------------------------------------------

      subroutine Pcirclef (xc,yc,radius)

      common /factors/  ax,ay,bx,by

      radiusx=radius*ax
      radiusy=radius*ay
      write (77,*) xc*ax+bx,yc*ay+by,
     &             (radiusx+radiusy)/2.,0,360,' a f s'

      return
      end
